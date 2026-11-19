"""Internal three-arm candidate stage for four independent training replicates.

CPU/RF only. Reads sources.json produced by common_slot runner.parents().
No training, sampling or Graph search is launched here. Candidate identities
are saved before score lookup. Existing chemistry implementations are reused;
the only retrieval change is a collision-free tuple-based permutation stream.
"""
import argparse
import json
import pickle
import random
import time
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import QED

import v2_unified_raw8_candidates as unified

common = unified.common
ARMS = ('morgan', 'task_contribution', 'permuted_contribution')
OUT = common.DEFAULT_DEV / 'independent_training256_v1'
HISTORY_NAMES = (
    'common_slot_retrieval_v1', 'dual_precursor_probe_v1',
    'official_syntwins_small_comparison_v1', 'unified_raw8_v1',
    'route_aware_raw8_v1', 'tap_vio_raw8_v1', 'reused_raw96_v1',
    'syntwins_depth3_raw96_v1', 'independent_training256_v1',
)


def identity(source):
    return {**common.source_identity(source),
            'replicate': source['replicate'], 'qed': source['qed']}


def validate_sources(rows):
    if len(rows) != 256 or len({s['sample_id'] for s in rows}) != 256:
        raise ValueError('Expected exactly 256 original raw slots with unique IDs')
    if Counter(s['task'] for s in rows) != Counter(gsk3b=128, jnk3=128):
        raise ValueError('Expected 128 original raw slots per task')
    for ti, task in enumerate(common.THRESHOLDS):
        selected = [s for s in rows if s['task'] == task]
        reps = sorted({s['replicate'] for s in selected})
        if reps != [0, 1]:
            raise ValueError('Expected training replicate IDs 0 and 1 per task')
        for rep in reps:
            cohort = [s for s in selected if s['replicate'] == rep]
            if sorted(s['position'] for s in cohort) != list(range(64)):
                raise ValueError(f'Expected positions 0..63 for {task} replicate {rep}')
        for s in selected:
            if s['task_index'] != ti:
                raise ValueError('Unexpected task index')
            if not s['valid']:
                if s['graph_verified']:
                    raise ValueError('Invalid source cannot be verified solved')
                continue
            if common.canon(Chem.MolFromSmiles(s['canonical_smiles'])) != s['canonical_smiles']:
                raise ValueError('Valid source must use canonical single-component SMILES')
            g = s['graph_result']
            if g['status'] not in ('solved', 'unsolved_within_budget'):
                raise ValueError('Execution error cannot be treated as an unsolved source')
            if g['seed'] != s['search_seed'] or g['canonical_smiles'] != s['canonical_smiles']:
                raise ValueError('Parent search chemical identity/seed differs from the source')
            if bool(g['solved']) != bool(s['graph_verified']):
                raise ValueError('Parent solved state differs from saved-route verification')
            if not Path(s['graph_result_path']).exists():
                raise ValueError('Missing saved parent result')


def permutation(values, source, root_index, group_index):
    seed = (9302000 + source['task_index'] * 100000 + source['replicate'] * 10000
            + source['position'] * 100 + root_index * 10 + group_index)
    shuffled = list(values)
    random.Random(seed).shuffle(shuffled)
    return shuffled, seed


def load_history(dev, out):
    lookup, ledger = common.load_oracle_history(dev, out)
    files = []
    cohort_dir = dev / 'independent_training256_v1'
    if (cohort_dir / 'generation_complete.json').exists():
        import independent_training256_prepare as preparation
        old_dev, old_out = preparation.DEV, preparation.OUT
        try:
            preparation.DEV, preparation.OUT = dev, cohort_dir
            training_lookup = preparation.history()
        finally:
            preparation.DEV, preparation.OUT = old_dev, old_out
        for task, values in training_lookup.items():
            for smiles, record in values.items():
                old = lookup[task].get(smiles)
                if old and abs(old['score'] - record['score']) > 1e-10:
                    raise ValueError('Training and saved task predictions disagree')
                if old is None:
                    lookup[task][smiles] = record
        files.extend(str(cohort_dir / f"{s['task']}-rep{s['replicate']}" /
                         'oracle-only/feedback.jsonl') for s in common.read(cohort_dir / 'design.json')['runs'])
    for name in HISTORY_NAMES:
        directory = dev / name
        for path in (directory / 'new_scores.json', directory / 'raw_scores.json',
                     directory / 'syntwins_depth3_k10/new_scores.json'):
            if not path.exists() or path == out / 'new_scores.json':
                continue
            data = common.read(path)
            if 'records' in data:
                records = data['records']
                if isinstance(records, dict):
                    records = [dict(r, task=t) for t, rs in records.items() for r in rs]
                else:
                    records = [dict(r, task=r.get('task', 'jnk3')) for r in records]
                entries = [(r['task'], r['smiles'], r['score']) for r in records]
            else:
                entries = [(r['task'], r['canonical_smiles'], r['task_score'])
                           for r in data['rows'] if r['valid']]
            for task, smiles, value in entries:
                if task not in lookup:
                    continue
                old = lookup[task].get(smiles)
                if old and abs(old['score'] - value) > 1e-10:
                    raise ValueError(f'Conflicting saved oracle value in {path}')
                if old is None:
                    lookup[task][smiles] = {'score': float(value), 'source': str(path)}
            files.append(str(path))
    return lookup, ledger, files


def load_forward(dev, out):
    cache, _ = common.load_forward_history(dev, out)
    roots, files = {}, []
    for name in HISTORY_NAMES:
        directory = dev / name
        path = directory / 'forward_cache.json'
        if path.exists():
            for r in common.read(path)['records']:
                if cache.get(r['rxn_ref'], r['reactants'], r.get('excluded_source')) is None:
                    cache.add(r['rxn_ref'], r['reactants'], r['products'], r['provenance'],
                              r.get('excluded_source'), r.get('forward_error'))
            files.append(str(path))
        # Complete unified root lists only. Earlier one-side lists are insufficient.
        if name not in ('unified_raw8_v1', 'route_aware_raw8_v1', 'reused_raw96_v1',
                        'independent_training256_v1'):
            continue
        for path in sorted((directory / 'source_records').glob('*.plan.json')):
            p = common.read(path)
            if 'all_exact_replay_roots' not in p.get('costs', {}):
                continue
            key = p['source']['canonical_smiles']
            roots.setdefault(key, {'selected': [{k: r[k] for k in
                ('rxn_ref', 'precursors', 'missing_sides')} for r in p['roots']],
                'costs': p['costs'], 'provenance': str(path)})
    return cache, roots, files


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


# generate_plan below follows v2_unified_raw8_candidates.generate_plan, with
# native enumeration separated and permutation seed including replicate.

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


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dev', type=Path, default=common.DEFAULT_DEV)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--sources', type=Path)
    parser.add_argument('--stock', type=Path, default=common.DEFAULT_STOCK)
    parser.add_argument('--models', type=Path, default=common.DEFAULT_MODEL)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    OUT = args.output
    source_path = args.sources or OUT / 'sources.json'
    RDLogger.DisableLog('rdApp.*')
    rows = sorted(common.read(source_path)['rows'], key=lambda s:
                  (s['task_index'], s['replicate'], s['position']))
    validate_sources(rows)
    if args.validate_only:
        print(json.dumps({'valid_source_interface': True, 'raw_denominator': len(rows)}))
        return
    design = common.read(OUT / 'design.json')
    common.ARMS = ARMS
    if (OUT / 'candidate_stage.json').exists():
        previous = common.read(OUT / 'candidate_stage.json')
        if [identity(r['source']) for r in previous['sources']] != [identity(s) for s in rows]:
            raise ValueError('Completed candidate stage has different sources')
        print(json.dumps({'reused_completed_candidate_stage': str(OUT / 'candidate_stage.json')}))
        return
    t0 = time.monotonic()
    inventory, forward, retros = common.initialize_chemistry(
        args.dev / 'external_baselines/SynTwins', args.stock)
    unified.NATIVE['stock_keys'] = set(Path(unified.EXPECTED['stock']).read_text().splitlines())
    cache, saved_roots, forward_files = load_forward(args.dev, OUT)
    oracle, ledger, score_files = load_history(args.dev, OUT)
    initialization_seconds = time.monotonic() - t0
    t0 = time.monotonic()
    models = {}
    for task in common.THRESHOLDS:
        if any(s['task'] == task and s['valid'] and not s['graph_verified'] for s in rows):
            with (args.models / f'{task}_current.pkl').open('rb') as stream:
                models[task] = pickle.load(stream)
    model_load_seconds = time.monotonic() - t0
    ip = OUT / 'candidate_invocations.json'
    invocations = common.read(ip)['invocations'] if ip.exists() else []
    invocations.append({'initialization_seconds': initialization_seconds,
        'model_load_seconds': model_load_seconds, 'inventory_size': len(inventory),
        'loaded_forward_cache_records': len(cache.records),
        'loaded_oracle_scores': {t: len(v) for t, v in oracle.items()},
        'score_history_paths': score_files, 'forward_history_paths': forward_files})
    common.save(ip, {'invocations': invocations})
    completed = []
    for source in rows:
        prefix = OUT / 'source_records' / (
            f'{source["task"]}-rep{source["replicate"]}-r{source["position"]:03d}')
        pp, cp = prefix.with_suffix('.plan.json'), prefix.with_suffix('.complete.json')
        if cp.exists():
            result = common.read(cp)
            if identity(result['source']) != identity(source):
                raise ValueError('Resumed candidate source differs')
        elif not source['valid'] or source['graph_verified']:
            result = common.retained_source(source)
            common.save(cp, result)
        else:
            if pp.exists():
                plan = common.read(pp)
                if identity(plan['source']) != identity(source):
                    raise ValueError('Resumed candidate plan differs')
            else:
                plan = generate_plan(source, models[source['task']], inventory,
                                     forward, retros, cache, saved_roots)
                common.save(OUT / 'forward_cache.json', {'records': cache.dump()})
                common.save(pp, plan)
            result = common.score_plan(plan, models[source['task']], oracle, ledger, OUT, args.models)
            common.save(cp, result)
        completed.append(result)
        common.save(OUT / 'candidate_progress.json', {'processed': len(completed),
            'raw_denominator': len(rows), 'last_source': source['sample_id']})
        print(json.dumps({'processed': len(completed), 'source': source['sample_id'],
            'status': result['status'], 'selected': {
                a: bool(result['arms'][a]['best_candidate']) for a in ARMS}}), flush=True)
    requests = []
    for r in completed:
        s = r['source']
        for arm in ARMS:
            candidate = r['arms'][arm]['best_candidate']
            if candidate:
                requests.append({'sample_id': s['sample_id'] + '-' + arm,
                    'source_id': s['sample_id'], 'arm': arm, 'task': s['task'],
                    'replicate': s['replicate'], 'position': s['position'],
                    'SMILES': candidate['smiles'], 'search_seed': s['search_seed'],
                    'candidate': candidate})
    summaries = {a: {'candidate_eligible_slots': sum(bool(r['arms'][a]['best_candidate']) for r in completed),
        **{k: sum(r['arms'][a]['costs'].get(k, 0) for r in completed) for k in
           ('logical_forward_attempts', 'logical_score_slots', 'logical_source_attributions', 'ranking_seconds')}}
        for a in ARMS}
    new_work = {k: sum(r['costs'].get(k, 0) for r in completed) for k in
        ('source_RF_path_attributions', 'retro_template_calls', 'source_replay_forward_calls',
         'attribution_seconds', 'root_seconds', 'matching_seconds')}
    new_work.update({
        'candidate_forward_attempts': sum(r['arms'][a]['costs']['new_forward_attempts'] for r in completed for a in ARMS),
        'candidate_forward_seconds': sum(r['arms'][a]['costs']['new_forward_seconds'] for r in completed for a in ARMS),
        'candidate_task_predictions': {t: len(rs) for t, rs in ledger['records'].items()},
        'candidate_score_seconds': sum(b['score_seconds'] for b in ledger['score_batches']),
        'initialization_seconds_all_invocations': sum(i['initialization_seconds'] for i in invocations),
        'model_load_seconds_all_invocations': sum(i['model_load_seconds'] for i in invocations),
        'training': 0, 'Graph': 0})
    stage = {'design': design, 'source_file': str(source_path), 'raw_denominator': len(rows),
        'sources': completed, 'arms': summaries, 'requests': requests, 'new_work': new_work,
        'candidate_order': 'Four exact-replay roots; rank outer/root inner/canonical product order; at most eight unique QED-qualified scored products per arm and source.',
        'permutation_seed': '9302000 + task_index*100000 + replicate*10000 + position*100 + root_index*10 + group_index',
        'graph_replacement_rule': 'Only a saved-graph-verified solve replaces the original slot; same per-query Graph12/45 ceiling, adopted historical seed recorded by runner.',
        'cost_interpretation': 'Actual shared execution and logical per-arm use are separate; this is not equal end-to-end wallclock.',
        'dependence': 'Four new training runs, two per task; 64 outputs within each run share a checkpoint. Duplicates retained in the denominator and share cached chemical measurements.'}
    common.save(OUT / 'candidate_stage.json', stage)
    print(json.dumps({'completed': str(OUT / 'candidate_stage.json'), 'arms': summaries,
                      'graph_requests': len(requests), 'new_work': new_work}), flush=True)


if __name__ == '__main__':
    main()

