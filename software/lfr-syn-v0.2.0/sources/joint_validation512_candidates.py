"""Bounded score-blind candidate policies for the new 512 fixed source slots.

Reuse the existing reaction/root/ranking implementation, then build Mix from
precursor proposals before any candidate task scores are read. This stage is
CPU/RF only. It does not train, sample, or make Graph2Edits predictions.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import pickle
import random
import time

ARMS = ('morgan', 'task_contribution', 'permuted_contribution', 'mix')
DEFAULT_DEV = Path('/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development')


def proposals(roots, arm):
    """The unchanged rank-outer/root-inner ordering; no task-score input."""
    for rank in range(4):
        for root in roots:
            missing = root['missing_sides']
            indices = [0, 0]
            if len(missing) == 1:
                indices[missing[0]] = rank
            elif len(missing) == 2:
                indices = [rank // 2, rank % 2]
            else:
                continue
            choices = [root['sides'][side]['choices'][arm] for side in range(2)]
            if any(indices[side] >= len(choices[side]) for side in range(2)):
                continue
            yield {'stream': arm, 'rank': rank, 'root_index': root['root_index'],
                   'rxn_ref': root['rxn_ref'], 'missing_side_count': len(missing),
                   'reactants': [choices[side][indices[side]]['smiles'] for side in range(2)]}


def alternating_proposals(roots):
    streams = [iter(proposals(roots, arm)) for arm in ('task_contribution', 'morgan')]
    active = [True, True]
    turn = 0
    while any(active):
        side = turn % 2
        turn += 1
        if not active[side]:
            continue
        try:
            yield next(streams[side])
        except StopIteration:
            active[side] = False


def execute_mix(roots, evaluate, max_attempts=16, max_products=8):
    """Apply chemistry/QED only; RF values never determine pool enrollment.

    evaluate(proposal) supplies evaluated_products in canonical order and a
    forward record. Duplicate reaction+ordered-reactants proposals advance
    their stream without spending a second logical combination allowance.
    """
    attempts, planned, duplicate_proposals = [], [], []
    seen_proposals, seen_products = set(), set()
    for proposal in alternating_proposals(roots):
        if len(attempts) >= max_attempts or len(planned) >= max_products:
            break
        key = (proposal['rxn_ref'], tuple(proposal['reactants']))
        if key in seen_proposals:
            duplicate_proposals.append(proposal)
            continue
        seen_proposals.add(key)
        record = evaluate(proposal)
        accepted_here = []
        for candidate in record['evaluated_products']:
            if len(planned) >= max_products:
                break
            if candidate['smiles'] in seen_products:
                continue
            seen_products.add(candidate['smiles'])
            if candidate['quality_pass']:
                planned.append(candidate)
                accepted_here.append(candidate['smiles'])
        attempts.append({**proposal, **record, 'enrolled_products': accepted_here})
    return {'planned_candidates': planned, 'forward_attempts': attempts,
            'skipped_duplicate_proposals': duplicate_proposals}


def permutation(values, source, root_index, group_index):
    seed = (10202000 + source['task_index'] * 1000000 + source['replicate'] * 10000
            + source['position'] * 100 + root_index * 10 + group_index)
    result = list(values)
    random.Random(seed).shuffle(result)
    return result, seed


def add_mix_plan(plan, forward, cache, stock_keys, out):
    from rdkit import Chem
    from rdkit.Chem import QED
    import independent_training256_candidates as old
    source = plan['source']
    costs = {'logical_forward_attempts': 0, 'new_forward_attempts': 0,
             'cached_forward_attempts': 0, 'new_forward_errors': 0,
             'logical_score_slots': 0, 'new_forward_seconds': 0.,
             'logical_source_attributions': 1,
             'ranking_seconds': sum(plan['arms'][arm]['costs']['ranking_seconds']
                                    for arm in ('morgan', 'task_contribution')),
             'ranking_reused_from_same_source_root_records': True}

    def evaluate(proposal):
        ref, reactants = proposal['rxn_ref'], proposal['reactants']
        costs['logical_forward_attempts'] += 1
        saved = cache.get(ref, reactants, source['canonical_smiles'])
        fresh = saved is None
        if fresh:
            started = time.monotonic()
            products, error = set(), None
            try:
                for product_set in forward[ref].RunReactants(
                        tuple(Chem.MolFromSmiles(smiles) for smiles in reactants)):
                    for product in product_set:
                        smiles = old.common.canon(Chem.Mol(product))
                        if smiles:
                            products.add(smiles)
            except Exception as exc:
                error = repr(exc)
                costs['new_forward_errors'] += 1
            costs['new_forward_seconds'] += time.monotonic() - started
            costs['new_forward_attempts'] += 1
            cache.add(ref, reactants, products, str(out), forward_error=error)
            saved = cache.get(ref, reactants, source['canonical_smiles'])
        else:
            costs['cached_forward_attempts'] += 1
        evaluated, excluded = [], Counter()
        for smiles in sorted(saved['products']):
            if smiles == source['canonical_smiles']:
                excluded['source'] += 1
                continue
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                excluded['invalid'] += 1
                continue
            if Chem.MolToInchiKey(molecule) in stock_keys:
                excluded['stock'] += 1
                continue
            qed = float(QED.qed(molecule))
            evaluated.append({'smiles': smiles, 'qed': qed, 'root_index': proposal['root_index'],
                              'rank': proposal['rank'], 'rxn_ref': ref, 'reactants': reactants,
                              'proposal_stream': proposal['stream'],
                              'missing_side_count': proposal['missing_side_count'],
                              'quality_pass': qed >= source['qed'] - .05})
        return {'cache_hit': not fresh, 'forward_provenance': saved['provenance'],
                'products': saved['products'], 'forward_error': saved.get('forward_error'),
                'evaluated_products': evaluated, 'excluded': dict(excluded)}

    started = time.monotonic()
    mix = execute_mix(plan['roots'], evaluate)
    costs['logical_score_slots'] = len(mix['planned_candidates'])
    mix['costs'] = costs
    mix['enrollment_rule'] = ('Before RF lookup, alternate contribution/Morgan precursor proposals; '
                              'deduplicate reaction+ordered-reactants; total16 combinations and8 '
                              'unique QED-qualified products, including cache hits.')
    plan['arms']['mix'] = mix
    plan['costs']['mix_plan_seconds'] = time.monotonic() - started
    plan['costs']['plan_generation_seconds'] += plan['costs']['mix_plan_seconds']
    return plan


def identity(source, common):
    return {**common.source_identity(source), 'replicate': source['replicate'],
            'cohort_id': source['cohort_id'], 'qed': source['qed']}


def validate_sources(rows, design, common):
    from rdkit import Chem
    expected = Counter({(spec['task'], spec['replicate']):
                        design['generation']['raw_draws_per_run'] for spec in design['runs']})
    if len(rows) != design['raw_denominator'] or len({s['sample_id'] for s in rows}) != len(rows):
        raise ValueError('Raw source coverage or slot identities differ from the planned new study')
    if Counter((s['task'], s['replicate']) for s in rows) != expected:
        raise ValueError('Incomplete source-training runs cannot replace the full planned denominator')
    for task, replicate in expected:
        selected = [s for s in rows if (s['task'], s['replicate']) == (task, replicate)]
        if sorted(s['position'] for s in selected) != list(range(expected[(task, replicate)])):
            raise ValueError(f'Unexpected raw source positions: {task} replicate {replicate}')
    for source in rows:
        if not source['valid']:
            if source['graph_verified']:
                raise ValueError('Invalid source cannot have an audited route')
            continue
        if common.canon(Chem.MolFromSmiles(source['canonical_smiles'])) != source['canonical_smiles']:
            raise ValueError('Source canonical identity differs')
        result = source['graph_result']
        if (result.get('status') not in ('solved', 'unsolved_within_budget') or
                result['seed'] != source['search_seed'] or
                result['canonical_smiles'] != source['canonical_smiles'] or
                bool(result['solved']) != bool(source['graph_verified'])):
            raise ValueError('Completed source route status/chemical identity/seed is inconsistent')
        if not Path(source['graph_result_path']).is_file():
            raise ValueError('Saved source route result is missing')


def joint_failed(source, thresholds):
    return bool(source['valid'] and (not source['graph_verified'] or
                                    source['task_score'] < thresholds[source['task']]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dev', type=Path, default=DEFAULT_DEV)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--sources', type=Path)
    parser.add_argument('--stock', type=Path)
    parser.add_argument('--models', type=Path)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    from rdkit import RDLogger
    import independent_training256_candidates as old
    common = old.common
    out = args.output or args.dev / 'joint_validation512_v1'
    source_path = args.sources or out / 'sources.json'
    stock = args.stock or common.DEFAULT_STOCK
    models_path = args.models or common.DEFAULT_MODEL
    RDLogger.DisableLog('rdApp.*')
    design = common.read(out / 'design.json')
    rows = sorted(common.read(source_path)['rows'], key=lambda s:
                  (s['task_index'], s['replicate'], s['position']))
    validate_sources(rows, design, common)
    if args.validate_only:
        print(json.dumps({'raw_denominator': len(rows), 'source_training_runs': len(design['runs']),
                          'valid_joint_failed': sum(joint_failed(s, common.THRESHOLDS) for s in rows)}))
        return
    final_path = out / 'candidate_stage.json'
    if final_path.exists():
        saved = common.read(final_path)
        if [identity(r['source'], common) for r in saved['sources']] != [identity(s, common) for s in rows]:
            raise ValueError('Completed candidate output belongs to different sources')
        print(json.dumps({'reuse_complete_candidate_stage': str(final_path)}), flush=True)
        return
    old.OUT, old.permutation = out, permutation
    common.ARMS = ARMS
    started = time.monotonic()
    inventory, forward, retros = common.initialize_chemistry(args.dev / 'external_baselines/SynTwins', stock)
    stock_keys = set(Path(old.unified.EXPECTED['stock']).read_text().splitlines())
    old.unified.NATIVE['stock_keys'] = stock_keys
    cache, saved_roots, forward_files = old.load_forward(args.dev, out)
    for path in sorted((out / 'source_records').glob('*.plan.json')):
        saved = common.read(path)
        if 'all_exact_replay_roots' in saved.get('costs', {}):
            saved_roots.setdefault(saved['source']['canonical_smiles'], {
                'selected': [{key: root[key] for key in ('rxn_ref', 'precursors', 'missing_sides')}
                             for root in saved['roots']],
                'costs': saved['costs'], 'provenance': str(path)})
    oracle, ledger, score_files = old.load_history(args.dev, out)
    initialization_seconds = time.monotonic() - started
    started = time.monotonic()
    models = {task: pickle.load((models_path / f'{task}_current.pkl').open('rb'))
              for task in common.THRESHOLDS
              if any(s['task'] == task and joint_failed(s, common.THRESHOLDS) for s in rows)}
    model_load_seconds = time.monotonic() - started
    invocation_path = out / 'candidate_invocations.json'
    invocations = common.read(invocation_path)['invocations'] if invocation_path.exists() else []
    invocations.append({'started_unix': time.time(), 'initialization_seconds': initialization_seconds,
                        'model_load_seconds': model_load_seconds, 'inventory_size': len(inventory),
                        'loaded_forward_cache_records': len(cache.records),
                        'score_history_paths': score_files, 'forward_history_paths': forward_files,
                        'loaded_oracle_scores': {task: len(records) for task, records in oracle.items()}})
    common.save(invocation_path, {'invocations': invocations})
    results = []
    for source in rows:
        prefix = out / 'source_records' / f"{source['cohort_id']}-r{source['position']:03d}"
        plan_path, complete_path = prefix.with_suffix('.plan.json'), prefix.with_suffix('.complete.json')
        if complete_path.exists():
            result = common.read(complete_path)
            if identity(result['source'], common) != identity(source, common):
                raise ValueError('Resumed candidate source differs')
        elif not joint_failed(source, common.THRESHOLDS):
            result = common.retained_source(source)
            result['status'] = 'retained_joint_success' if source['valid'] else 'retained_invalid'
            common.save(complete_path, result)
        else:
            if plan_path.exists():
                plan = common.read(plan_path)
                if identity(plan['source'], common) != identity(source, common):
                    raise ValueError('Resumed candidate plan source differs')
            else:
                plan = old.generate_plan(source, models[source['task']], inventory,
                                         forward, retros, cache, saved_roots)
                add_mix_plan(plan, forward, cache, stock_keys, out)
                common.save(out / 'forward_cache.json', {'records': cache.dump()})
                common.save(plan_path, plan)
            result = common.score_plan(plan, models[source['task']], oracle, ledger, out, models_path)
            for arm in ARMS:
                scored = result['arms'][arm]['scored_candidates']
                for candidate in scored:
                    candidate['task_pass'] = bool(candidate['task_score'] >=
                                                  max(source['task_score'], common.THRESHOLDS[source['task']]))
                result['arms'][arm]['eligible_candidates'] = sorted(
                    [c for c in scored if c['task_pass']],
                    key=lambda c: (-c['task_score'], -c['qed'], c['smiles']))
            common.save(complete_path, result)
        results.append(result)
        common.save(out / 'candidate_progress.json', {'processed': len(results), 'raw_denominator': len(rows),
                                                      'last_source': source['sample_id'], 'updated_unix': time.time()})
        print(json.dumps({'processed': len(results), 'source': source['sample_id'], 'status': result['status'],
                          'eligible': {arm: len(result['arms'][arm].get('eligible_candidates', []))
                                       for arm in ARMS}}), flush=True)
    summaries = {arm: {'candidate_eligible_slots': sum(bool(r['arms'][arm]['best_candidate']) for r in results),
                       'candidate_eligible_second_slots': sum(len(r['arms'][arm].get('eligible_candidates', [])) >= 2
                                                               for r in results),
                       **{key: sum(r['arms'][arm]['costs'].get(key, 0) for r in results)
                          for key in ('logical_forward_attempts', 'logical_score_slots', 'logical_source_attributions')}}
                 for arm in ARMS}
    new_work = {key: sum(r['costs'].get(key, 0) for r in results) for key in
                ('source_RF_path_attributions', 'retro_template_calls', 'source_replay_forward_calls',
                 'attribution_seconds', 'root_seconds', 'matching_seconds', 'mix_plan_seconds')}
    new_work.update({'candidate_forward_attempts': sum(r['arms'][arm]['costs'].get('new_forward_attempts', 0)
                                                      for r in results for arm in ARMS),
                     'candidate_forward_seconds': sum(r['arms'][arm]['costs'].get('new_forward_seconds', 0)
                                                     for r in results for arm in ARMS),
                     'candidate_task_predictions': {task: len(records) for task, records in ledger['records'].items()},
                     'candidate_score_seconds': sum(batch['score_seconds'] for batch in ledger['score_batches']),
                     'initializations': invocations, 'training': 0, 'Graph': 0})
    common.save(final_path, {'design': design, 'source_file': str(source_path), 'raw_denominator': len(rows),
                             'sources': results, 'arms': summaries, 'new_work': new_work,
                             'interpretation': 'All512 slots, eight source-training runs; candidates enrolled before RF '
                                               'lookup. Shared actual work and per-policy logical budgets are separate.'})
    print(json.dumps({'candidate_stage_complete': str(final_path), 'arms': summaries, 'new_work': new_work}), flush=True)


if __name__ == '__main__':
    main()
