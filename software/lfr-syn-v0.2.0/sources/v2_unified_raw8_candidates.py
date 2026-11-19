"""Unified one/two-missing-precursor retrieval with a disclosed official fix."""
from collections import Counter
import importlib.util
import random
import signal
import sys
import time
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import QED

import v2_common_slot_retrieval_candidates as common
from v2_official_syntwins_small_comparison import CountedReaction, ComputationBound, COMMIT
from v2_task_attribution_graph import EXPECTED

OUT = common.DEFAULT_DEV / 'unified_raw8_v1'
INTERNAL = common.ARMS
ARMS = INTERNAL + ('syntwins_compatible_depth1',)
_old_initialize = common.initialize_chemistry
_old_load_history = common.load_oracle_history
_old_load_forward = common.load_forward_history
NATIVE = {}


def load_history(dev, out):
    history, ledger = _old_load_history(dev, out)
    previous, _ = _old_load_history(dev, dev / 'common_slot_retrieval_v1')
    for task in history:
        history[task].update(previous[task])
    for name in ('dual_precursor_probe_v1', 'official_syntwins_small_comparison_v1'):
        path = dev / name / 'new_scores.json'
        for r in common.read(path)['records']:
            task = r.get('task', 'jnk3')
            old = history[task].get(r['smiles'])
            if old and abs(old['score']-r['score']) > 1e-10:
                raise ValueError('Saved oracle values disagree')
            history[task][r['smiles']] = {'score': r['score'], 'source': str(path)}
    return history, ledger


def load_forward(dev, out):
    cache, _ = _old_load_forward(dev, out)
    for name in ('common_slot_retrieval_v1', 'dual_precursor_probe_v1'):
        path = dev / name / 'forward_cache.json'
        for r in common.read(path)['records']:
            if cache.get(r['rxn_ref'], r['reactants'], r.get('excluded_source')) is None:
                cache.add(r['rxn_ref'], r['reactants'], r['products'], r['provenance'],
                          r.get('excluded_source'), r.get('forward_error'))
    # Older selected single-side root lists are not complete unified root lists.
    return cache, {}


def initialize(repo, stock_path):
    result = _old_initialize(repo, stock_path)
    inventory = result[0]
    NATIVE['stock_keys'] = set(Path(EXPECTED['stock']).read_text().splitlines())
    sys.path.insert(0, str(repo))
    from syntwins.data_utils import SynTwin_DataLoader
    import subprocess
    if subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip() != COMMIT:
        raise ValueError('SynTwins version changed')
    original = (repo / 'syntwins/syntree.py').read_text()
    old_slice = "bblocks = [bblocks[i] for i in ranks[:filter_k]]"
    new_slice = "bblocks = [bblocks[i] for i in (ranks if filter_k < 0 else ranks[:filter_k])]"
    old_role = """            for i, work_ref in rxns:
                ref_kdtree = self.kdtrees[work_ref][i]"""
    new_role = """            for i, work_ref in rxns:
                if i >= len(self.kdtrees[work_ref]) or i >= len(self.kd_bblocks[work_ref]):
                    continue
                ref_kdtree = self.kdtrees[work_ref][i]"""
    if original.count(old_slice) != 1 or original.count(old_role) != 1:
        raise ValueError('Expected native compatibility edit locations differ')
    compatible = original.replace(old_slice, new_slice).replace(old_role, new_role)
    copy_path = OUT / 'syntree_compatible.py'
    copy_path.write_text(compatible)
    spec = importlib.util.spec_from_file_location('syntwins.syntree_compatible', copy_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    data = common.DEFAULT_DEV / 'official_syntwins_small_comparison_v1/data'
    if data.joinpath('building_blocks.smi').read_text().splitlines() != list(inventory):
        raise ValueError('Native preprocessing inventory differs from the common inventory')
    for name in ('templates.csv', 'retro_templates.csv'):
        if data.joinpath(name).read_bytes() != repo.joinpath('data', name).read_bytes():
            raise ValueError('Native preprocessing templates differ')
    start = time.monotonic()
    loader = SynTwin_DataLoader(data_dir=str(data), _nBits=1024, verbose=False)
    if set(loader.work_rxns) != set(loader.template_refs):
        raise ValueError('Incomplete native reaction preprocessing')
    for reaction in loader.work_rxns.values():
        for ids in reaction.available_reactants_ids:
            if any(i < 0 or i >= len(inventory) for i in ids):
                raise ValueError('Native stock role index outside inventory')
    state = {'calls': Counter()}
    for reaction in loader.retro_rxns:
        reaction.rxn = CountedReaction(reaction.rxn, 'retro', state)
        reaction.work_rxn.rxn = CountedReaction(reaction.work_rxn.rxn, 'forward', state)
    NATIVE.update(loader=loader, tree_class=module.SynTree, state=state)
    common.save(OUT / 'native_initialization.json', {
        'seconds': time.monotonic()-start, 'fresh_stock_preprocessing': False,
        'prior_preprocessing': str(data), 'prior_first_initialization_seconds': common.read(data.parent / 'initialization.json')['seconds'],
        'official_commit': COMMIT, 'original_source_modified': False,
        'Graph_stock': EXPECTED['stock'], 'Graph_stock_key_count': len(NATIVE['stock_keys']),
        'compatibility_copy': str(copy_path), 'edits': [{'before': old_slice, 'after': new_slice}, {'before': old_role, 'after': new_role}],
        'zero_stock_reactions': [r for r,v in loader.work_rxns.items() if v.n_potential_reactions == 0],
    })
    return result


def find_roots(source, inventory, forward, retros):
    started = time.monotonic()
    molecule = Chem.MolFromSmiles(source['canonical_smiles'])
    seen, roots, costs = set(), [], Counter()
    for ref, retro in retros:
        if forward[ref].GetNumReactantTemplates() != 2:
            continue
        costs['retro_template_calls'] += 1
        try:
            pairs = retro.RunReactants((molecule,))
        except Exception:
            costs['retro_errors'] += 1
            continue
        for pair in pairs:
            precursors = [common.canon(Chem.Mol(m)) for m in pair]
            if len(precursors) != 2 or not all(precursors):
                continue
            key = (ref, tuple(precursors))
            if key in seen:
                continue
            seen.add(key)
            costs['source_replay_forward_calls'] += 1
            try:
                recovered = any(source['canonical_smiles'] in {common.canon(Chem.Mol(p)) for p in ps}
                    for ps in forward[ref].RunReactants(tuple(Chem.MolFromSmiles(s) for s in precursors)))
            except Exception:
                costs['source_replay_errors'] += 1
                recovered = False
            if recovered:
                missing = [i for i,s in enumerate(precursors) if Chem.MolToInchiKey(Chem.MolFromSmiles(s)) not in NATIVE['stock_keys']]
                roots.append({'rxn_ref': ref, 'precursors': precursors, 'missing_sides': missing})
    groups = {1: [], 2: []}
    for count in (1, 2):
        used_refs = set()
        for r in sorted(roots, key=lambda r: (r['rxn_ref'], tuple(r['precursors']))):
            if len(r['missing_sides']) == count and r['rxn_ref'] not in used_refs:
                groups[count].append(r)
                used_refs.add(r['rxn_ref'])
    selected = groups[1][:2] + groups[2][:2]
    extras = groups[1][2:] + groups[2][2:]
    selected.extend(extras[:4-len(selected)])
    selected.sort(key=lambda r: (len(r['missing_sides']), r['rxn_ref'], tuple(r['precursors'])))
    counts = Counter(len(r['missing_sides']) for r in roots)
    return selected, dict(costs, root_seconds=time.monotonic()-started,
        exact_replay_root_counts=dict(counts), grouped_root_counts={str(k):len(v) for k,v in groups.items()},
        selected_root_counts=dict(Counter(len(r['missing_sides']) for r in selected)),
        all_exact_replay_roots=roots)


def native_arm(source, inventory):
    state = NATIVE['state']
    state['calls'] = Counter()
    seed = 9298100+100*source['task_index']+source['position']
    import numpy as np
    np.random.seed(seed)
    random.seed(seed)
    tree = NATIVE['tree_class'](NATIVE['loader'], _nBits=1024, max_depth=1, k=5, verbose=False)
    started = time.monotonic()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(ComputationBound('60 second native-call request')))
    signal.alarm(60)
    try:
        outputs = tree.get_twins(source['canonical_smiles'], fast=True, return_as_list=False)
        status, error = 'complete', None
    except (Exception, ComputationBound) as exc:
        outputs, status, error = {}, 'execution_limitation', repr(exc)
    finally:
        signal.alarm(0)
    seconds = time.monotonic()-started
    planned, seen, excluded = [], set(), Counter()
    for native_smiles, record in outputs.items():
        mol = Chem.MolFromSmiles(native_smiles)
        smiles = common.canon(mol)
        if not smiles or smiles in seen:
            excluded['invalid_duplicate'] += 1
            continue
        seen.add(smiles)
        if smiles == source['canonical_smiles'] or Chem.MolToInchiKey(mol) in NATIVE['stock_keys']:
            excluded['source_or_stock'] += 1
            continue
        qed = float(QED.qed(mol))
        if qed < source['qed']-.05:
            excluded['qed'] += 1
            continue
        if len(planned) < 8:
            planned.append({'smiles': smiles, 'qed': qed, 'native_record': record})
        else:
            excluded['over_score_allowance'] += 1
    return {'planned_candidates': planned, 'native_outputs': outputs, 'native_status': status,
        'native_error': error, 'excluded': dict(excluded), 'seed': seed,
        'costs': {'logical_forward_attempts': state['calls']['forward'], 'new_forward_attempts': state['calls']['forward'],
            'logical_score_slots': len(planned), 'ranking_seconds': 0., 'new_forward_seconds': 0.,
            'native_total_seconds': seconds, 'native_retro_calls': state['calls']['retro'],
            'native_forward_includes_source_replay': True, 'logical_source_attributions': 0}}


def generate_plan(source, model, inventory, forward, retros, cache, saved_roots):
    started = time.monotonic()
    roots, root_costs = find_roots(source, inventory, forward, retros)
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
                  'new_forward_errors': 0, 'logical_source_attributions': int(a != 'morgan')} for a in INTERNAL}
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
            seed = 9297100+1000*source['task_index']+100*source['position']+10*ri+gi
            random.Random(seed).shuffle(shuffled)
            permuted.update(zip(ordered, shuffled))
            permutations.append({'environments': ordered, 'weights': values, 'permuted_weights': shuffled, 'seed': seed})
        side_records = []
        for side in range(2):
            if side not in sides:
                choices = {a: [{'smiles': root['precursors'][side], 'similarity': 1., 'retained_weight': 0.}] for a in INTERNAL}
                side_records.append({'side': side, 'preserved_stock': True, 'choices': choices})
                continue
            start = time.monotonic()
            query = forward[root['rxn_ref']].GetReactantTemplate(side)
            matches = [(s, float(DataStructs.TanimotoSimilarity(precursor_fps[side][0], fp)), envs)
                       for s,(m,fp,envs) in inventory.items() if m.HasSubstructMatch(query)]
            matching_seconds += time.monotonic()-start
            choices = {}
            for arm in INTERNAL:
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
    for arm in INTERNAL:
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
                    if Chem.MolToInchiKey(pm) in NATIVE['stock_keys']:
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
    arms[ARMS[-1]] = native_arm(source, inventory)
    return {'source': source, 'source_qed': source['qed'], 'status': 'candidate_plan_saved_before_oracle_lookup',
        'roots': root_records, 'arms': arms, 'prediction_from_tree_paths': prediction, 'source_bias': bias,
        'costs': {**root_costs, 'source_RF_path_attributions': 1, 'attribution_seconds': attr_seconds,
                  'matching_seconds': matching_seconds, 'plan_generation_seconds': time.monotonic()-started}}


def main():
    if len(sys.argv) != 1:
        raise SystemExit('This experiment reads unified_raw8_v1; run without arguments.')
    sys.argv.extend(['--output', str(OUT), '--sources', str(OUT / 'sources.json')])
    common.ARMS = ARMS
    common.initialize_chemistry = initialize
    common.load_oracle_history = load_history
    common.load_forward_history = load_forward
    common.generate_plan = generate_plan
    common.main()
    # The inherited main records candidate forward counts for all arms. Native
    # forward includes replay; keep it separate to avoid mixing cost definitions.
    path = OUT / 'candidate_stage.json'
    stage = common.read(path)
    rows = stage['sources']
    stage['new_work']['candidate_forward_attempts_internal'] = sum(
        row['arms'][arm]['costs']['new_forward_attempts'] for row in rows for arm in INTERNAL)
    stage['new_work']['official_total_forward_including_replay'] = sum(
        row['arms'][ARMS[-1]]['costs']['new_forward_attempts'] for row in rows)
    stage['new_work']['official_retro_template_calls'] = sum(
        row['arms'][ARMS[-1]]['costs'].get('native_retro_calls', 0) for row in rows)
    stage['new_work']['actual_total_forward_including_internal_source_replay'] = (
        stage['new_work']['candidate_forward_attempts_internal']+
        stage['new_work']['official_total_forward_including_replay']+
        stage['new_work']['source_replay_forward_calls'])
    stage['new_work']['actual_total_retro'] = (stage['new_work']['retro_template_calls']+
        stage['new_work']['official_retro_template_calls'])
    stage['new_work'].pop('candidate_forward_attempts', None)
    stage['cost_interpretation'] = common.read(OUT / 'design.json')['cost']
    stage['candidate_order'] = {'internal': common.read(OUT / 'design.json')['replacement_rule'],
                                'official': 'Native output order after the disclosed compatibility edits.'}
    stage['official_execution_limitations'] = [row['source']['sample_id'] for row in rows
        if row['arms'][ARMS[-1]].get('native_status') == 'execution_limitation']
    common.save(path, stage)


if __name__ == '__main__':
    main()
