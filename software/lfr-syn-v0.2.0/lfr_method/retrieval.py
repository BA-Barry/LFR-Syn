"""Core extracted from retained independent_training256_candidates.py; see code_origins.json."""
from collections import Counter
from pathlib import Path
import time
from rdkit import Chem, DataStructs
from rdkit.Chem import QED
from . import common, roots as unified
from .permutation import permutation
ARMS = common.ARMS
OUT = Path(".")


def find_roots(source, inventory, forward, retros, saved_roots):
    previous = saved_roots.get(source['canonical_smiles'])
    if previous:
        metadata = {k: v for k, v in previous['costs'].items() if k in
                    ('exact_replay_root_counts', 'grouped_root_counts',
                     'selected_root_counts', 'all_exact_replay_roots')}
        return previous['selected'], dict(metadata, source_replay_forward_calls=0,
            retro_template_calls=0, root_seconds=0.,
            root_discovery_provenance=previous['provenance'])
    roots, costs = unified.find_roots(source, inventory, forward, retros)
    saved_roots[source['canonical_smiles']] = {
        'selected': roots, 'costs': costs, 'provenance': str(OUT / 'source_records')}
    return roots, costs


def generate_plan(source, model, inventory, forward, retros, cache, saved_roots):
    started = time.monotonic()
    roots, root_costs = find_roots(source, inventory, forward, retros, saved_roots)
    mol = Chem.MolFromSmiles(source['canonical_smiles'])
    _, environments = common.fingerprints(mol)
    start = time.monotonic()
    prediction, bias, contribution = common.local_contributions(model, common.features(mol))
    if abs(prediction-source['task_score']) > 1e-8:
        raise ValueError('Attribution and saved source oracle differ')
    folded = Counter(e % 2048 for e in environments)
    weights = {e: float(max(0., contribution[e % 2048])/folded[e % 2048]) for e in environments}
    attr_seconds = time.monotonic()-start
    costs = {a: {'logical_forward_attempts': 0, 'new_forward_attempts': 0, 'cached_forward_attempts': 0,
                  'logical_score_slots': 0, 'ranking_seconds': 0., 'new_forward_seconds': 0.,
                  'new_forward_errors': 0, 'logical_source_attributions': int(a != 'morgan')} for a in ARMS}
    root_records, matching_seconds = [], 0.
    for ri, root in enumerate(roots):
        sides = root['missing_sides']
        eligible = []
        precursor_fps = [common.fingerprints(Chem.MolFromSmiles(s)) for s in root['precursors']]
        for side in range(2):
            eligible.append(environments & precursor_fps[side][1] if side in sides else set())
        groups = [eligible[0] & eligible[1], eligible[0]-eligible[1], eligible[1]-eligible[0]]
        permuted, permutations = {}, []
        for gi, group in enumerate(groups):
            ordered = sorted(group)
            values = [weights[e] for e in ordered]
            shuffled = list(values)
            shuffled, seed = permutation(values, source, ri, gi)
            permuted.update(zip(ordered, shuffled))
            permutations.append({'environments': ordered, 'weights': values, 'permuted_weights': shuffled, 'seed': seed,
                'seed_coordinates': [source['task_index'], source['replicate'], source['position'], ri, gi]})
        side_records = []
        for side in range(2):
            if side not in sides:
                choices = {a: [{'smiles': root['precursors'][side], 'similarity': 1., 'retained_weight': 0.}] for a in ARMS}
                side_records.append({'side': side, 'preserved_stock': True, 'choices': choices})
                continue
            start = time.monotonic()
            query = forward[root['rxn_ref']].GetReactantTemplate(side)
            matches = [(s, float(DataStructs.TanimotoSimilarity(precursor_fps[side][0], fp)), envs)
                       for s,(m,fp,envs) in inventory.items() if m.HasSubstructMatch(query)]
            matching_seconds += time.monotonic()-start
            choices = {}
            for arm in ARMS:
                start = time.monotonic()
                ws = weights if arm == 'task_contribution' else permuted
                ranked = [(s,sim,0. if arm == 'morgan' else float(sum(ws[e] for e in eligible[side] if e in envs))) for s,sim,envs in matches]
                ranked.sort(key=lambda r: (-r[2],-r[1],r[0]))
                choices[arm] = [{'smiles': s, 'similarity': sim, 'retained_weight': w} for s,sim,w in ranked[:4 if len(sides)==1 else 2]]
                costs[arm]['ranking_seconds'] += time.monotonic()-start
            side_records.append({'side': side, 'preserved_stock': False, 'matching_stock_count': len(matches), 'choices': choices})
        root_records.append({**root, 'root_index': ri, 'sides': side_records, 'permutations': permutations,
            'B_C_same_lists': all([c['smiles'] for c in sr['choices']['task_contribution']] == [c['smiles'] for c in sr['choices']['permuted_contribution']] for sr in side_records)})
    arms = {}
    for arm in ARMS:
        planned, attempts, seen = [], [], set()
        for rank in range(4):
            for root in root_records:
                if len(planned) >= 8:
                    break
                sides = root['missing_sides']
                indices = [0,0]
                if len(sides) == 1:
                    indices[sides[0]] = rank
                else:
                    indices = [rank//2, rank%2]
                choices = [root['sides'][i]['choices'][arm] for i in range(2)]
                if any(indices[i] >= len(choices[i]) for i in range(2)):
                    continue
                reactants = [choices[i][indices[i]]['smiles'] for i in range(2)]
                costs[arm]['logical_forward_attempts'] += 1
                saved = cache.get(root['rxn_ref'], reactants, source['canonical_smiles'])
                fresh = saved is None
                if fresh:
                    start = time.monotonic()
                    products, error = set(), None
                    try:
                        for product_set in forward[root['rxn_ref']].RunReactants(tuple(Chem.MolFromSmiles(s) for s in reactants)):
                            for product in product_set:
                                if (s := common.canon(Chem.Mol(product))):
                                    products.add(s)
                    except Exception as exc:
                        error = repr(exc)
                        costs[arm]['new_forward_errors'] += 1
                    costs[arm]['new_forward_attempts'] += 1
                    costs[arm]['new_forward_seconds'] += time.monotonic()-start
                    cache.add(root['rxn_ref'], reactants, products, str(OUT), forward_error=error)
                    saved = cache.get(root['rxn_ref'], reactants, source['canonical_smiles'])
                else:
                    costs[arm]['cached_forward_attempts'] += 1
                evaluated = []
                for s in sorted(saved['products']):
                    if len(planned) >= 8:
                        break
                    if s == source['canonical_smiles'] or s in seen:
                        continue
                    seen.add(s)
                    pm = Chem.MolFromSmiles(s)
                    if Chem.MolToInchiKey(pm) in unified.NATIVE['stock_keys']:
                        continue
                    qed = float(QED.qed(pm))
                    candidate = {'smiles': s, 'qed': qed, 'root_index': root['root_index'], 'rank': rank,
                                 'rxn_ref': root['rxn_ref'], 'reactants': reactants, 'missing_side_count': len(sides),
                                 'quality_pass': qed >= source['qed']-.05}
                    evaluated.append(candidate)
                    if candidate['quality_pass']:
                        planned.append(candidate)
                attempts.append({'rxn_ref': root['rxn_ref'], 'root_index': root['root_index'], 'rank': rank,
                    'reactants': reactants, 'cache_hit': not fresh, 'forward_provenance': saved['provenance'],
                    'products': saved['products'], 'forward_error': saved.get('forward_error'), 'evaluated_products': evaluated})
            if len(planned) >= 8:
                break
        costs[arm]['logical_score_slots'] = len(planned)
        arms[arm] = {'planned_candidates': planned, 'forward_attempts': attempts, 'costs': costs[arm]}
    return {'source': source, 'source_qed': source['qed'], 'status': 'candidate_plan_saved_before_oracle_lookup',
        'roots': root_records, 'arms': arms, 'prediction_from_tree_paths': prediction, 'source_bias': bias,
        'costs': {**root_costs, 'source_RF_path_attributions': 1, 'attribution_seconds': attr_seconds,
                  'matching_seconds': matching_seconds, 'plan_generation_seconds': time.monotonic()-started}}
