"""Independently audit terminal saved five-arm endpoints, without science models.

Reads native saved vectors/provenance and the separate native route receipt.
Does not import the study evaluator, its statistics, Chemprop, RDKit or Torch.
SciPy is used only as a numerical reference for Student distribution functions.
The sole primary remains fixed-qualification B minus source-relative B.
"""
import argparse
import ast
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time
import traceback

POLICIES = ('Original', 'B_joint6plus6', 'C_joint6plus6', 'B_fixedtask6plus6', 'C_fixedtask6plus6')
TASKS = ('gsk3b', 'jnk3')
ENDPOINTS = ('external_joint', 'joint', 'solved')
FILTERS = ('all', 'no_pains', 'no_pains_brenk', 'mw_le500', 'source_qed_tolerance')
SOURCE_FIELDS = ('sample_id', 'task', 'task_index', 'replicate', 'cohort_id', 'position',
                 'canonical_smiles', 'valid', 'task_score', 'qed', 'search_seed', 'graph_verified', 'graph_result_path')
HISTORIES = ('graphga_transfer512_v1', 'graphga_fixedtask128_v1', 'graphga_fixedtask512_v1',
             'graphga_candidateallocation128_v1', 'graphga_retainedprecursor128_v1')
REFERENCE = {'repository': 'https://github.com/Tibogoss/CombiMOTS',
             'published_commit': '3f644b727982a95b4f87115ccb8d6b4c918a02ba',
             'software': 'Chemprop1.5.2', 'tasks': ['gsk3b_activity', 'jnk3_activity'], 'models': 10}
PACKAGES = {'chemprop': '1.5.2', 'torch': '2.4.1+cu121', 'rdkit': '2023.9.6', 'numpy': '1.26.4'}
TERMINAL = ('complete', 'complete_with_execution_limits')
BOOTSTRAPS = 100000
BOOTSTRAP_SEED = 2026101801


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def require(value, message):
    if not value:
        raise ValueError(message)


def near(a, b, tolerance=1e-10):
    return (type(a) in (int, float) and type(b) in (int, float)
            and math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tolerance)


def probability(value):
    return type(value) in (int, float) and math.isfinite(value) and 0 <= value <= 1


def compare(actual, reported, label='saved arithmetic'):
    """Compare independently recomputed fields, permitting unrelated metadata."""
    if isinstance(actual, dict):
        require(isinstance(reported, dict), label + ': dictionary required')
        for key, value in actual.items():
            require(key in reported, label + ': missing ' + key)
            compare(value, reported[key], label + '.' + key)
    elif isinstance(actual, list):
        require(isinstance(reported, list) and len(actual) == len(reported), label + ': list membership differs')
        for i, (a, b) in enumerate(zip(actual, reported)):
            compare(a, b, label + '[%d]' % i)
    elif type(actual) in (int, float):
        require(near(actual, reported), label + ': numeric mismatch')
    else:
        require(type(actual) is type(reported) and actual == reported, label + ': value/type mismatch')


def bounds(row, endpoint):
    value = row[endpoint]
    if type(value) is bool:
        return [int(value), int(value)]
    require(value is None, 'Explicit boolean/null endpoint required')
    lo, hi = row.get('endpoint_lower_' + endpoint, 0), row.get('endpoint_upper_' + endpoint, 1)
    require(lo in (0, 1) and hi in (0, 1) and lo <= hi, 'Invalid endpoint uncertainty bounds')
    return [int(lo), int(hi)]


def design_check(design):
    require(design['experiment'] == 'graphga_qualification_independent1024_v1', 'Wrong independent experiment')
    require(tuple(design['policies']) == POLICIES and design['raw_source_positions'] == 1024
            and design['source_positions_per_run'] == 64 and design['independent_optimizer_runs'] == 16,
            'Prespecified complete five-arm design differs')
    specs = {r['cohort_id']: r for r in design['runs']}
    require(len(specs) == len(design['runs']) == 16
            and {(r['task'], r['replicate']) for r in specs.values()} == {(t, r) for t in TASKS for r in range(8)},
            'Eight independent optimizations per fixed task required')
    for spec in specs.values():
        i = 2 * spec['replicate'] + TASKS.index(spec['task'])
        require(spec['task_index'] == TASKS.index(spec['task'])
                and spec['cohort_id'] == 'graphga-qualconfirm-%s-rep%d' % (spec['task'], spec['replicate'])
                and spec['generation_seed'] == 2026101501 + i
                and spec['initialization_seed'] == 2026101601 + i
                and spec['terminal_sampling_seed'] == 2026101701 + i
                and spec['search_seed_base'] == 12500000 + 10000 * i, 'Prospective run identity/seeds differ')
    require(design['method']['source_RF_cutoffs'] == {'gsk3b': .1, 'jnk3': .05}
            and design['method']['candidate_calls'] == [6, 6]
            and design['method']['candidate_shared_soft_seconds'] == 45., 'Retained task/search recipe differs')
    require(design['primary_contrast'].startswith('B_fixedtask6plus6 minus B_joint6plus6,')
            and 'no guiding RF conjunct' in design['primary_contrast'], 'Sole primary comparison changed')
    return specs


def selection_check(study, fixed, rows, route):
    design = read(study/'design.json')
    specs = design_check(design)
    require(fixed['design'] == design and fixed['experiment'] == design['experiment'], 'External fixed design differs')
    sources = read(study/'sources_with_parents.json')['rows']
    measured = read(study/'method_result.json')
    require(read(study/'method_execution.json')['status'] in TERMINAL, 'Method must be terminal')
    source_index = {s['sample_id']: s for s in sources}
    require(len(sources) == len(source_index) == 1024 and measured['selected_raw_n'] == 1024, 'Source/raw denominator differs')
    for source in sources:
        spec = specs.get(source['cohort_id'])
        require(spec is not None and tuple(source[k] for k in ('task', 'replicate', 'task_index'))
                == tuple(spec[k] for k in ('task', 'replicate', 'task_index')), 'Unexpected source optimizer identity')
        require(type(source['position']) is int and source['position'] in range(64)
                and source.get('requested_search_seed', source['search_seed']) == spec['search_seed_base'] + source['position'],
                'Raw draw position/requested parent seed differs')
        require(source['sample_id'] == source['cohort_id'] + '-r%03d' % source['position'], 'Sample primary key is not its original raw draw')
        require(type(source['valid']) is bool and source['valid'] == bool(source['canonical_smiles']), 'Source-valid identity differs')
    order = list(specs)
    ordered = sorted(sources, key=lambda s: (order.index(s['cohort_id']), s['position']))
    require(fixed['sources'] == ordered, 'Full fixed source objects changed')
    for cid in specs:
        group = [s for s in sources if s['cohort_id'] == cid]
        require(len(group) == 64 and {s['position'] for s in group} == set(range(64)), 'Missing/duplicated raw draw position')
    expected = {(sid, p) for sid in source_index for p in POLICIES}
    indexes = []
    for collection in (measured['rows'], fixed['rows'], rows):
        index = {(r['sample_id'], r['policy']): r for r in collection}
        require(len(collection) == len(index) == 5120 and set(index) == expected, 'Missing/duplicate/extra five-arm membership')
        indexes.append(index)
    native, selected, enriched = indexes
    require(route['status'] == 'passed' and Path(route['out']).resolve() == study.resolve()
            and route['summary']['raw_positions'] == 1024
            and route['summary']['independent_optimizer_runs'] == 16
            and route['summary']['five_policy_rows'] == 5120, 'Complete native route audit missing')
    native_checks = {(r['sample_id'], r['policy']): r for r in route['policy_row_checks']}
    require(len(route['policy_row_checks']) == len(native_checks) == 5120 and set(native_checks) == expected
            and all(r['passed'] is True for r in native_checks.values()), 'Five-arm native route/decision proof incomplete')
    for key, row in native.items():
        require(all(selected[key].get(k) == v for k, v in row.items()), 'Method selection changed before external evaluation')
        require(selected[key].get('source_valid') is source_index[key[0]]['valid'], 'Derived fixed source-valid flag differs')
        require(all(enriched[key].get(k) == v for k, v in selected[key].items()), 'Enrichment altered fixed method decisions')
        s = source_index[key[0]]
        for a, b in (('task', 'task'), ('replicate', 'replicate'), ('cohort_id', 'cohort_id'), ('position', 'position'),
                     ('source_smiles', 'canonical_smiles'), ('source_score', 'task_score'), ('source_qed', 'qed')):
            require(row[a] == s[b], 'Policy source identity differs: ' + a)
        declared = {k: s.get(k) for k in SOURCE_FIELDS}
        require(row['source_identity'] in (s, declared), 'Copied source identity differs')
        require(type(row['scientific_endpoint_available']) is bool and type(row['final_structure_known']) is bool
                and type(row['replaced']) is bool, 'Scientific endpoint availability must be explicit')
        require(probability(row['source_score']) and probability(row['final_score']), 'Saved RF probability invalid')
        require(row['status'] not in ('processing', 'pending'), 'Live policy cannot be enrolled as terminal')
        require(native_checks[key].get('unknown_preserved', False) == (not row['scientific_endpoint_available']),
                'Unknown endpoint differs from native route receipt')
    for cid in specs:
        folder = study/'method_outputs'/cid
        require(read(folder/'method_execution.json')['status'] in TERMINAL, 'Cohort scientific work is live')
        require(read(folder/'sources_with_parents.json')['rows'] == [s for s in sources if s['cohort_id'] == cid],
                'Aggregate/cohort source objects differ')
        crows = read(folder/'method_result.json')['rows']
        require(len(crows) == 320 and {(r['sample_id'], r['policy']): r for r in crows}
                == {k: r for k, r in native.items() if r['cohort_id'] == cid}, 'Aggregate/cohort complete decisions differ')
    return design, ordered, indexes


def vector_check(study, sources, fixed_rows, predictions, evaluation, analysis, expected_packages):
    ext, dev = study/'external_fixed_outputs_v1', study.parent
    require(predictions['tasks'] == REFERENCE['tasks'] and predictions['models'] == 10
            and predictions['model_reference'] == REFERENCE and evaluation['model_reference'] == REFERENCE,
            'Published asset declaration/task order differs')
    require(analysis['runtime']['packages'] == expected_packages, 'Actual external feature/runtime packages differ')
    declared_runtime = evaluation['named_history_feature_runtime_evidence']
    require({e['study'] for e in declared_runtime} == set(HISTORIES) and len(declared_runtime) == 5,
            'Named historical feature/runtime declarations incomplete')
    files, anchors = {}, {}

    def get(path):
        path = Path(path).resolve()
        require(path.is_relative_to(dev.resolve()), 'Vector provenance outside authorized development directory')
        if str(path) not in files:
            files[str(path)] = read(path)
        return path, files[str(path)]

    for evidence in declared_runtime:
        path, saved = get(evidence['runtime_path'])
        require(path.is_relative_to((dev/evidence['study']).resolve()), 'Historical runtime declared under wrong study')
        actual = saved['packages'] if 'packages' in saved else saved['runtime']['packages']
        require(actual == evidence['packages'] == expected_packages, 'Physical historical feature/runtime versions differ')

    required = {s['canonical_smiles'] for s in sources if s['valid']}
    required.update(r['final_smiles'] for r in fixed_rows if r['final_smiles'])
    require(set(predictions['predictions']) == required, 'Prediction identities differ from all fixed source/final structures')

    def native_vector(path, smiles, fold):
        path, record = get(path)
        if record.get('tasks') is not None:
            require(record['tasks'] == REFERENCE['tasks'], 'Provenance task-column order differs')
        if record.get('model_reference') is not None:
            require(record['model_reference'] == REFERENCE, 'Provenance model reference differs')
        else:
            anchor_path = path.parent/'evaluation_design.json'
            aggregate_path = path.parent/'predictions.json'
            key = str(anchor_path)
            if key not in anchors:
                _, header = get(anchor_path)
                _, aggregate = get(aggregate_path)
                require(header.get('repository') == REFERENCE['repository']
                        and header.get('published_commit') == REFERENCE['published_commit']
                        and header.get('tasks') == REFERENCE['tasks']
                        and 'Chemprop1.5.2' in str(header.get('models')) and '10-checkpoint' in str(header.get('models'))
                        and aggregate['tasks'] == REFERENCE['tasks'] and aggregate['models'] == 10,
                        'Legacy physical vector lacks actual compatible saved model declaration')
                anchors[key] = {'aggregate': aggregate, 'aggregate_path': str(aggregate_path), 'vectors_checked': 0}
            anchors[key]['vectors_checked'] += 1
        entries = record['predictions']
        if isinstance(entries, list):
            require(len(entries) == len(record['smiles']) and len(set(record['smiles'])) == len(record['smiles']),
                    'Legacy vector parallel membership differs')
            entries = dict(zip(record['smiles'], entries))
        value = entries[smiles]
        if isinstance(value, dict):
            require(record.get('models', 10) == 10 and len(value['folds']) == 10, 'Native aggregate checkpoint count differs')
            vector = value['folds'][fold]
            if 'mean' in value:
                require(all(near(value['mean'][i], statistics.fmean(v[i] for v in value['folds']), 1e-7)
                            for i in (0, 1)), 'Native saved mean conflicts with actual fold vectors')
        else:
            require(record['fold'] == fold, 'Native single-fold index differs')
            vector = value
        require(isinstance(vector, list) and len(vector) == 2 and all(map(probability, vector)), 'Invalid native two-task vector')
        if record.get('model_reference') is None:
            require(anchors[str(path.parent/'evaluation_design.json')]['aggregate']['predictions'][smiles]['folds'][fold] == vector,
                    'Legacy saved physical vector differs from its original aggregate')
        return vector

    # Conflicts are rejected for every required identity in every declared cache,
    # including an equivalent origin the merged analysis did not choose first.
    history_values = [{} for _ in range(10)]
    historical_paths = []
    for group in evaluation['history_groups']:
        for path in group['paths']:
            path, data = get(path)
            historical_paths.append(str(path))
            entries = data['predictions']
            if isinstance(entries, list):
                entries = dict(zip(data['smiles'], entries))
            for smiles in required.intersection(entries):
                folds = range(10) if isinstance(entries[smiles], dict) else (data['fold'],)
                for fold in folds:
                    vector = native_vector(path, smiles, fold)
                    previous = history_values[fold].get(smiles)
                    require(previous is None or all(near(a, b, 1e-7) for a, b in zip(previous, vector)),
                            'Required same-model historical vector conflict exceeds 1e-7')
                    history_values[fold].setdefault(smiles, vector)

    counts = Counter()
    physical_refs = set()
    for smiles, record in predictions['predictions'].items():
        require(len(record['folds']) == len(record['fold_provenance']) == 10, 'Complete ten-fold vector/provenance required')
        require(all(len(v) == 2 and all(map(probability, v)) for v in record['folds']), 'Invalid saved prediction vector')
        for i in (0, 1):
            require(near(record['mean'][i], statistics.fmean(v[i] for v in record['folds']), 1e-7)
                    and near(record['fold_sd'][i], statistics.pstdev(v[i] for v in record['folds']), 1e-7),
                    'Saved ten-fold mean/SD arithmetic differs')
        for fold, (vector, origin) in enumerate(zip(record['folds'], record['fold_provenance'])):
            require(origin['recorded_smiles'] == smiles, 'Saved canonical vector provenance identity differs')
            path = Path(origin['selected_path']).resolve()
            require(native_vector(path, smiles, fold) == vector, 'Analysis vector differs from actual physical origin')
            physical_refs.add(str(path))
            require(origin['kind'] in ('historical', 'current_missing_only'), 'Unknown fold provenance kind')
            if origin['kind'] == 'historical':
                require(str(path) in historical_paths and smiles in history_values[fold]
                        and all(near(a, b, 1e-7) for a, b in zip(vector, history_values[fold][smiles])),
                        'Historical origin not in actual declared compatible cache')
            else:
                require(path == (ext/('new_fold_%02d_predictions.json' % fold)).resolve()
                        and smiles not in history_values[fold], 'Missing-only inference duplicated compatible historical measurement')
            for equivalent in origin.get('equivalent_cache_paths', []):
                other = native_vector(equivalent, smiles, fold)
                require(all(near(a, b, 1e-7) for a, b in zip(vector, other)), 'Equivalent physical vector conflict exceeds 1e-7')
                physical_refs.add(str(Path(equivalent).resolve()))
            counts[origin['kind']] += 1

    cost = analysis['execution_and_reuse']['physical_new_inference_cumulative']
    fold_costs, total, seconds = [], 0, 0.
    for fold in range(10):
        path, record = get(ext/('new_fold_%02d_predictions.json' % fold))
        require(record['model_reference'] == REFERENCE and record['tasks'] == REFERENCE['tasks'] and record['fold'] == fold,
                'New native prediction file identity differs')
        batches = record['batches']
        membership = [s for b in batches for s in b['smiles']]
        require(len(membership) == len(set(membership)) == len(record['predictions'])
                and set(membership) == set(record['predictions']) and set(membership) <= required
                and sum(b['structures'] for b in batches) == len(membership), 'Physical inference batch membership differs')
        require(all(type(b['seconds']) in (int, float) and math.isfinite(b['seconds']) and b['seconds'] >= 0 for b in batches),
                'Invalid physical inference duration')
        require(all(s not in history_values[fold] for s in membership), 'Saved new work duplicated declared old exact predictions')
        duration = sum(b['seconds'] for b in batches)
        fold_costs.append({'fold': fold, 'new_structure_vectors': len(membership), 'saved_batch_seconds': duration, 'path': str(path)})
        total += len(membership)
        seconds += duration
    compare({'new_fold_structure_vectors_cumulative': total, 'saved_batch_seconds_cumulative': seconds,
             'fold_files': fold_costs}, cost, 'Physical inference work')
    require(predictions['logical_checkpoint_structure_predictions'] == len(required) * 10
            and predictions['reused_historical_fold_structure_predictions'] == counts['historical']
            and predictions['model_structure_predictions_this_invocation'] + predictions['reused_current_fold_structure_predictions']
                == counts['current_missing_only']
            and predictions['model_structure_predictions_this_invocation'] <= total
            and predictions['historical_recomputed'] == 0, 'Logical/physical saved prediction provenance accounting differs')
    require(analysis['execution_and_reuse']['training_generation_RF_forward_parent_work_this_external_stage'] == 0,
            'Unexpected non-external scientific work')
    return {'prediction_structures': len(required), 'saved_native_fold_vectors_checked': len(required) * 10,
            'physical_vector_provenance_files_read': len(physical_refs), 'provenance_kind_counts': dict(counts),
            'required_vector_history_files_read': len(set(historical_paths)), 'physical_new_vectors_cumulative': total,
            'physical_new_saved_batch_seconds': seconds,
            'legacy_model_declaration_anchors': [{'evaluation_design': k, 'aggregate_path': v['aggregate_path'],
                                                  'vectors_checked': v['vectors_checked']} for k, v in anchors.items()]}


def row_arithmetic(rows, predictions, profiles, descriptors):
    require(profiles['rdkit_version'] == descriptors['rdkit_version']
            and profiles['profiles'] == descriptors['profiles'], 'Saved descriptor outputs differ from physical cache')
    similarity = {(r['source'], r['final']): r['similarity'] for r in descriptors['source_final_similarities']}
    values = predictions['predictions']
    checks = []
    for row in rows:
        index = TASKS.index(row['task'])
        source = values.get(row['source_smiles'])
        final = values.get(row['final_smiles'])
        source_mean = statistics.fmean(v[index] for v in source['folds']) if source else 0.
        final_mean = statistics.fmean(v[index] for v in final['folds']) if final else 0.
        source_profile = profiles['profiles'].get(row['source_smiles'])
        final_profile = profiles['profiles'].get(row['final_smiles'])
        require(near(row['source_external_probability'], source_mean)
                and row['source_descriptors'] == source_profile
                and row['source_external_score_provenance'] == (source['fold_provenance'] if source else []),
                'Source external arithmetic/profile/provenance differs')
        if source_profile is not None:
            require(source_profile.get('canonical_smiles', row['source_smiles']) == row['source_smiles']
                    and near(source_profile['qed'], row['source_qed']), 'Source saved profile identity/QED differs')
        available = row['scientific_endpoint_available']
        for endpoint in ('solved', 'joint'):
            require(available == (type(row[endpoint]) is bool), 'Unknown scientific route/RF result was upgraded or lost')
            bounds(row, endpoint)
        if available:
            require(row['final_structure_known'] == bool(row['final_smiles']), 'Scientific final structure availability differs')
            require(row['descriptors'] == final_profile and near(row['external_probability'], final_mean)
                    and near(row['external_change'], final_mean - source_mean)
                    and row['external_fold_probabilities'] == ([v[index] for v in final['folds']] if final else [])
                    and row['external_score_provenance'] == (final['fold_provenance'] if final else []),
                    'Final external arithmetic/profile/provenance differs')
            require(row['external_score_available'] is True and row['external_valid'] == bool(row['final_smiles']),
                    'Known external score availability differs')
            if final:
                require(near(row['external_fold_sd'], statistics.pstdev(v[index] for v in final['folds'])), 'Final checkpoint SD differs')
            for name in ('qed', 'molecular_weight'):
                expected = final_profile[name] - source_profile[name] if final_profile and source_profile else None
                compare(expected, row[name + '_delta'], 'Row property delta')
            if final_profile is not None:
                require(final_profile.get('canonical_smiles', row['final_smiles']) == row['final_smiles']
                        and near(final_profile['qed'], row['final_qed'])
                        and row['final_scaffold'] == final_profile.get('murcko_scaffold'), 'Final profile/QED/scaffold differs')
            require(row['source_morgan_similarity'] == similarity[(row['source_smiles'], row['final_smiles'])], 'Saved source/final similarity differs')
            external = bool(row['solved'] and final_mean >= .5)
            require(row['external_joint'] is external, 'Primary external endpoint must omit guiding RF conjunct')
            eligible = row['source_valid'] and row['final_score'] >= (.1 if row['task'] == 'gsk3b' else .05)
            require(row['joint'] is bool(row['solved'] and eligible), 'RF plus route endpoint arithmetic differs')
            if row['replaced']:
                require(row['solved'] and row['final_smiles'] != row['source_smiles'] and eligible
                        and row['final_qed'] >= row['source_qed'] - .05, 'Replacement task/QED qualification differs')
                if row['policy'] in POLICIES[1:3]:
                    require(row['final_score'] >= row['source_score'], 'Strict arm RF guarantee violated')
        else:
            require(not row['replaced'] and not row['final_structure_known'] and row['external_joint'] is None
                    and row['external_valid'] is None and row['external_score_available'] is False,
                    'Unknown final policy was promoted to a known output')
            for field in ('external_probability', 'external_change', 'external_fold_sd', 'qed_delta',
                          'molecular_weight_delta', 'descriptors', 'final_scaffold', 'source_morgan_similarity'):
                require(row[field] is None, 'Operational fallback used as known scientific property: ' + field)
            require(row['external_fold_probabilities'] == row['external_score_provenance'] == []
                    and bounds(row, 'external_joint') == [0, bounds(row, 'solved')[1]], 'Unknown final vectors/bounds differ')
            require(row['operational_fallback_smiles'] == row['final_smiles']
                    and near(row['operational_fallback_external_probability'], final_mean)
                    and row['operational_fallback_descriptors'] == final_profile,
                    'Separate operational fallback provenance differs')
        checks.append({'sample_id': row['sample_id'], 'policy': row['policy'], 'passed': True,
                       'scientific_endpoint_available': available, 'external_joint_bounds': bounds(row, 'external_joint')})
    return checks


def distribution(values):
    values = sorted(float(v) for v in values if v is not None)
    if not values:
        return {'n': 0, 'mean': None, 'median': None}
    return {'n': len(values), 'mean': sum(values) / len(values), 'median': statistics.median(values),
            'min': values[0], 'max': values[-1]}


def summary(rows):
    n = len(rows)
    result = {'raw_n': n, 'replaced_slots': sum(r['replaced'] for r in rows),
              'route_execution_unknown_slots': sum(not r['scientific_endpoint_available'] for r in rows),
              'distinct_source_task_structures': len({(r['task'], r['source_smiles']) for r in rows if r['source_smiles']}),
              'distinct_final_task_structures_observed': len({(r['task'], r['final_smiles']) for r in rows if r['final_structure_known']}),
              'known_route_failures_within_saved_recipe': sum(r['solved'] is False for r in rows),
              'candidate_query_status_counts': dict(Counter(a['query']['status'] for r in rows for a in r.get('attempts', []) if 'query' in a)),
              'candidate_budget_unclosed_attempts': sum(a.get('query', {}).get('status') == 'unsolved_within_budget' for r in rows for a in r.get('attempts', []))}
    for endpoint in ENDPOINTS:
        success = [r for r in rows if r[endpoint] is True]
        identities = {(r['task'], r['final_smiles']) for r in success}
        raw = [sum(bounds(r, endpoint)[i] for r in rows) for i in (0, 1)]
        possible = sum(bounds(r, endpoint)[1] for r in rows if r[endpoint] is None)
        result[endpoint] = {'known_successes': len(success), 'unknown_slots': sum(r[endpoint] is None for r in rows),
            'raw_count_bounds': raw, 'raw_rate_bounds': [v / n for v in raw],
            'distinct_task_structure_known_successes': len(identities),
            'distinct_task_structure_count_bounds': [len(identities), len(identities) + possible],
            'distinct_nonempty_successful_task_scaffolds': len({(r['task'], r['final_scaffold']) for r in success if r.get('final_scaffold')}),
            'acyclic_successful_task_structures': len({(r['task'], r['final_smiles']) for r in success if not r.get('final_scaffold')})}
    profiles = [r['descriptors'] for r in rows if r.get('descriptors') is not None]
    result['profile_diagnostics'] = {'raw_denominator': n, 'observed_final_profile_slots': len(profiles),
        'pains_flagged_slots': sum(bool(p['pains']) for p in profiles), 'brenk_flagged_slots': sum(bool(p['brenk']) for p in profiles),
        'mw_above500_slots': sum(p['molecular_weight'] > 500 for p in profiles)}
    return result


def properties(rows):
    result = {}
    for label, chosen in (('all_positions', rows), ('replacement_positions', [r for r in rows if r['replaced']])):
        n = len(chosen)
        observed = [r for r in chosen if r['scientific_endpoint_available']]
        unknown = [r for r in chosen if not r['scientific_endpoint_available']]
        rf = [r['final_score'] - r['source_score'] for r in observed]
        decreased = sum(v < 0 for v in rf)
        result[label] = {'raw_denominator': n, 'unknown_final_positions': len(unknown),
            'RF_change_observed': distribution(rf),
            'external_probability_change_observed': distribution(r.get('external_change') for r in observed),
            'QED_change_observed': distribution(r.get('qed_delta') for r in observed),
            'molecular_weight_change_observed': distribution(r.get('molecular_weight_delta') for r in observed),
            'RF_decreased_count_bounds': [decreased, decreased + len(unknown)],
            'RF_decreased_fraction_bounds': [decreased / n, (decreased + len(unknown)) / n] if n else None,
            'mean_RF_change_identification_bounds': [(sum(rf) + sum(i - r['source_score'] for r in unknown)) / n for i in (0, 1)] if n else None,
            'mean_external_change_identification_bounds': [(sum(r.get('external_change') or 0 for r in observed)
                + sum(i - r['source_external_probability'] for r in unknown)) / n for i in (0, 1)] if n else None}
    result['audited_routes_below_external_cutoff'] = sum(r['solved'] is True and r['external_probability'] < .5 for r in rows)
    return result


def concentration(rows):
    successes = [r for r in rows if r['external_joint'] is True]
    products = Counter((r['task'], r['final_smiles']) for r in successes)
    replacements = Counter((r['task'], r['final_smiles']) for r in rows if r['replaced'])
    cohorts = defaultdict(set)
    for row in successes:
        cohorts[(row['task'], row['final_smiles'])].add(row['cohort_id'])
    top = [v for k, v in products.most_common(2)]
    return {'raw_denominator': len(rows), 'known_successful_positions': len(successes),
        'distinct_successful_task_structures': len(products), 'top_one_successful_positions': top[0] if top else 0,
        'top_two_successful_positions': sum(top), 'top_one_fraction_of_raw_positions': top[0] / len(rows) if top else 0.,
        'top_two_fraction_of_known_successes': sum(top) / len(successes) if successes else None,
        'distinct_replacement_task_structures': len(replacements),
        'successful_products': [{'task': k[0], 'canonical_smiles': k[1], 'positions': v, 'optimizer_runs': len(cohorts[k])}
                                for k, v in products.most_common()],
        'products_successful_in_multiple_runs': sum(len(cids) > 1 for cids in cohorts.values())}


def filtered(rows):
    output = {}
    for name in FILTERS:
        good = []
        for row in rows:
            p, source = row.get('descriptors'), row.get('source_descriptors')
            if row['external_joint'] is not True or p is None:
                continue
            passes = (name == 'all' or (name == 'no_pains' and not p['pains'])
                or (name == 'no_pains_brenk' and not p['pains'] and not p['brenk'])
                or (name == 'mw_le500' and p['molecular_weight'] <= 500)
                or (name == 'source_qed_tolerance' and source is not None and p['qed'] >= source['qed'] - .05))
            if passes:
                good.append(row)
        possible = sum(bounds(r, 'external_joint')[1] for r in rows if r['external_joint'] is None)
        output[name] = {'raw_denominator': len(rows), 'external_joint_count_bounds': [len(good), len(good) + possible],
            'distinct_successful_task_structures': len({(r['task'], r['final_smiles']) for r in good}),
            'distinct_nonempty_successful_task_scaffolds': len({(r['task'], r['final_scaffold']) for r in good if r.get('final_scaffold')}),
            'acyclic_successful_task_structures': len({(r['task'], r['final_smiles']) for r in good if not r.get('final_scaffold')})}
    return output


def paired_runs(rows, design, terms, endpoint):
    index = {(r['sample_id'], r['policy']): r for r in rows}
    sources = [r for r in rows if r['policy'] == 'Original']
    runs = []
    for spec in design['runs']:
        chosen = [r for r in sources if r['cohort_id'] == spec['cohort_id']]
        lo = hi = unknown = positive = negative = 0
        for s in chosen:
            selected = [(weight, index[(s['sample_id'], policy)]) for policy, weight in terms]
            for weight, row in selected:
                a, b = bounds(row, endpoint)
                lo += weight * (a if weight > 0 else b)
                hi += weight * (b if weight > 0 else a)
            missing = any(row[endpoint] is None for weight, row in selected)
            unknown += missing
            if not missing:
                delta = sum(weight * int(row[endpoint]) for weight, row in selected)
                positive += delta > 0
                negative += delta < 0
        require(len(chosen) == 64, 'Paired contrast lost raw draw positions')
        runs.append({'cohort_id': spec['cohort_id'], 'task': spec['task'], 'replicate': spec['replicate'], 'raw_n': 64,
            'net_count_bounds': [lo, hi], 'effect_bounds': [lo / 64, hi / 64],
            'known_positive_positions': positive, 'known_negative_positions': negative, 'unknown_pairs': unknown})
    return runs


def scipy_inference(runs, *, primary, bootstrap_cache=None):
    grouped = {task: [r for r in runs if r['task'] == task] for task in TASKS}
    require(len(runs) == 16 and all(len(grouped[t]) == 8 for t in TASKS), 'Independent-run statistical membership differs')
    intervals = {t: [100 * statistics.fmean(r['effect_bounds'][i] for r in grouped[t]) for i in (0, 1)] for t in TASKS}
    unknown = sum(r['unknown_pairs'] for r in runs)
    result = {'task_equal_identification_bounds_pp': [statistics.fmean(intervals[t][i] for t in TASKS) for i in (0, 1)],
        'technical_unknown_pairs': unknown, 'task_effect_bounds_pp': intervals, 'welch': None, 'bootstrap_sensitivity': None}
    if unknown:
        result['interval_status'] = 'technical_unknowns_identification_bounds_only'
        return result
    # Reference implementation independent of the evaluator's incomplete-beta code.
    from scipy.stats import t as student
    effects = {t: [r['effect_bounds'][0] for r in grouped[t]] for t in TASKS}
    means = {t: statistics.fmean(effects[t]) for t in TASKS}
    variances = {t: statistics.variance(effects[t]) for t in TASKS}
    task_summaries = {}
    for task in TASKS:
        info = {'optimizer_runs': 8, 'effects_pp': [100 * x for x in effects[task]],
            'mean_effect_pp': 100 * means[task], 'run_variance': variances[task], 'interval_95_pp': None}
        if primary:
            info['p_two_sided'] = None
        if variances[task] == 0:
            info['interval_status'] = 'not_estimable_zero_observed_run_variance'
        else:
            se = math.sqrt(variances[task] / 8)
            critical = float(student.ppf(.975, 7))
            info.update(interval_status='complete_data_paired_run_t',
                        interval_95_pp=[100 * (means[task] - critical * se), 100 * (means[task] + critical * se)])
            if primary:
                info['p_two_sided'] = float(min(1., 2 * student.sf(abs(means[task] / se), 7)))
        task_summaries[task] = info
    mean = statistics.fmean(means.values())
    variance_terms = [variances[t] / 8 for t in TASKS]
    sum_variance = sum(variance_terms)
    welch = {'task_equal_mean_effect_pp': 100 * mean, 'task_summaries': task_summaries,
        'interval_95_pp': None, 'alpha': .05, 'standard_error_pp': 50 * math.sqrt(sum_variance)}
    if primary:
        welch['p_two_sided'] = None
    if sum_variance == 0:
        welch.update(interval_status='not_estimable_zero_observed_run_variance', degrees_of_freedom=None)
    else:
        se = math.sqrt(sum_variance) / 2
        df = sum_variance ** 2 / sum(v ** 2 / 7 for v in variance_terms)
        critical = float(student.ppf(.975, df))
        welch.update(interval_status='complete_data_welch', degrees_of_freedom=df,
            interval_95_pp=[100 * (mean - critical * se), 100 * (mean + critical * se)], t_statistic=mean / se)
        if primary:
            welch['p_two_sided'] = float(min(1., 2 * student.sf(abs(mean / se), df)))
        if primary:
            request = {'effects_by_task': effects, 'replicates': BOOTSTRAPS, 'seed': BOOTSTRAP_SEED,
                       'resampling_unit': 'complete paired optimizer run; never raw positions, structures or checkpoints'}
            cached = read(bootstrap_cache) if bootstrap_cache and Path(bootstrap_cache).exists() else None
            if cached is not None:
                require(cached['request'] == request and cached['status'] == 'complete', 'Existing audit bootstrap does not match paired-run effects')
                boot = cached['result']
            else:
                generator = random.Random(BOOTSTRAP_SEED)
                sampled = []
                for _ in range(BOOTSTRAPS):
                    first = sum(effects[TASKS[0]][generator.randrange(8)] for _ in range(8)) / 8
                    second = sum(effects[TASKS[1]][generator.randrange(8)] for _ in range(8)) / 8
                    sampled.append(50 * (first + second))
                sampled.sort()
                def quantile(q):
                    pos = (len(sampled) - 1) * q
                    a, b = math.floor(pos), math.ceil(pos)
                    return sampled[a] + (sampled[b] - sampled[a]) * (pos - a)
                boot = {'replicates': BOOTSTRAPS, 'seed': BOOTSTRAP_SEED,
                        'interval_95_pp': [quantile(.025), quantile(.975)],
                        'task_equal_mean_effect_pp': 100 * mean, 'resampling_unit': request['resampling_unit']}
                if bootstrap_cache:
                    save(bootstrap_cache, {'status': 'complete', 'request': request, 'result': boot,
                        'scientific_model_calls': 0, 'finished_unix': time.time()})
            result['bootstrap_sensitivity'] = boot
    result['welch'] = welch
    result['interval_status'] = welch['interval_status']
    return result


def arithmetic_check(rows, design, analysis, bootstrap_cache):
    arms = {p: [r for r in rows if r['policy'] == p] for p in POLICIES}
    compare({'experiment': design['experiment'], 'raw_source_positions': 1024, 'policy_rows': 5120,
        'independent_optimizer_runs': 16, 'independent_runs_per_fixed_task': 8,
        'primary_contrast': design['primary_contrast']}, analysis, 'Prespecified analysis')
    diagnostics = {'policies': {p: summary(arms[p]) for p in POLICIES},
        'by_cohort': {spec['cohort_id']: {p: summary([r for r in arms[p] if r['cohort_id'] == spec['cohort_id']]) for p in POLICIES} for spec in design['runs']},
        'by_task': {t: {p: summary([r for r in arms[p] if r['task'] == t]) for p in POLICIES} for t in TASKS},
        'property_changes': {p: properties(arms[p]) for p in POLICIES},
        'concentration': {p: concentration(arms[p]) for p in POLICIES},
        'structural_sensitivity': {p: filtered(arms[p]) for p in POLICIES}}
    source_counts = Counter((r['task'], r['source_smiles']) for r in arms['Original'] if r['source_smiles'])
    diagnostics['input_duplicates'] = {'raw_positions': 1024, 'valid_source_positions': sum(r['source_valid'] for r in arms['Original']),
        'distinct_task_structures': len(source_counts), 'repeated_valid_structure_positions': sum(source_counts.values()) - len(source_counts)}
    compare(diagnostics, analysis, 'Direct five-arm diagnostics')
    comparisons = [(((POLICIES[3], 1), (POLICIES[1], -1)), 'external_joint', True)]
    pairs = ((POLICIES[4], POLICIES[2]), (POLICIES[1], POLICIES[2]), (POLICIES[3], POLICIES[4]),
             (POLICIES[3], 'Original'), (POLICIES[4], 'Original'))
    comparisons.extend((((a, 1), (b, -1)), e, False) for a, b in pairs for e in ENDPOINTS)
    comparisons.extend((((POLICIES[3], 1), (POLICIES[1], -1)), e, False) for e in ('joint', 'solved'))
    interaction = ((POLICIES[3], 1), (POLICIES[4], -1), (POLICIES[1], -1), (POLICIES[2], 1))
    comparisons.extend((interaction, e, False) for e in ENDPOINTS)
    reported = [analysis['primary']] + analysis['exploratory_secondary_contrasts']
    require(len(reported) == len(comparisons) == 21, 'Primary/secondary contrast membership differs')
    recomputed = []
    for (terms, endpoint, primary), saved in zip(comparisons, reported):
        require(saved['terms'] == [{'policy': p, 'weight': w} for p, w in terms] and saved['endpoint'] == endpoint,
                'Prespecified primary/exploratory contrast order changed')
        require(saved['scope'].startswith('Prespecified sole primary' if primary else 'Prespecified exploratory secondary'),
                'Primary and exploratory labels were exchanged')
        runs = paired_runs(rows, design, terms, endpoint)
        infer = scipy_inference(runs, primary=primary, bootstrap_cache=bootstrap_cache if primary else None)
        compare({'raw_n': 1024, 'runs': runs, **infer}, saved, 'Independent paired-run/SciPy inference')
        if not primary and saved['welch'] is not None:
            require('p_two_sided' not in saved['welch'] and all('p_two_sided' not in s for s in saved['welch']['task_summaries'].values()),
                    'Exploratory contrast acquired undeclared significance test')
        recomputed.append({'terms': saved['terms'], 'endpoint': endpoint, 'primary': primary, 'runs': runs, **infer})
    return {'paired_run_contrasts_checked': len(recomputed), 'sole_primary': recomputed[0],
            'exploratory_secondary': recomputed[1:], 'direct_policy_diagnostics': diagnostics}


def ordering_check(study, evaluation, evaluation_code):
    ext = study/'external_fixed_outputs_v1'
    require(evaluation['raw_positions'] == 1024 and evaluation['optimizer_runs'] == 16
            and tuple(evaluation['policies']) == POLICIES and evaluation['primary_threshold'] == .5
            and evaluation['external_scores_used_in_selection'] is False and evaluation['fixed_outputs_no_reselection'] is True,
            'External design denominator/selection rule differs')
    module = ast.parse(Path(evaluation_code).read_text(encoding='utf-8-sig'))
    execute = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'execute')
    def name(node):
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Attribute):
            return name(node.value) + '.' + node.attr
        return ''
    calls = defaultdict(list)
    for node in ast.walk(execute):
        if isinstance(node, ast.Call):
            calls[name(node.func)].append(node.lineno)
    selection = min(calls['write_fixed'])
    later = ('validate_named_history_runtimes', 'adapter.history_groups', 'adapter.merge_history', 'original.obtain_predictions')
    require(all(selection < min(calls[k]) for k in later), 'Executed code reads scores/predictors before saving complete choices')
    stage = read(study/'stage_execution/external.json')
    require(stage['status'] == 'complete' and stage['exit_code'] == 0, 'External process did not terminate successfully')
    method = read(study/'method_execution.json')
    require(method['ended_unix'] <= stage['started_unix'], 'External stage preceded saved method decisions')
    fixed_time = (ext/'fixed_outputs.json').stat().st_mtime
    folds = [ext/('new_fold_%02d_predictions.json' % f) for f in range(10)]
    measured = [p for p in folds if read(p)['predictions']]
    require(stage['started_unix'] <= fixed_time <= stage['ended_unix'], 'Fixed selection not saved during measured external invocation')
    # On a resumed external invocation, fixed outputs may be rewritten identically.
    # Selection remains bound to the complete method objects and executed ordering.
    return {'executed_evaluation_code': str(evaluation_code), 'write_fixed_line': selection,
            'later_call_lines': {k: calls[k] for k in later}, 'method_ended_unix': method['ended_unix'],
            'external_started_unix': stage['started_unix'], 'external_ended_unix': stage['ended_unix'],
            'fixed_outputs_mtime': fixed_time, 'saved_nonempty_new_fold_mtimes': [p.stat().st_mtime for p in measured],
            'limitation': 'Actual code order, unchanged full objects and execution timestamps; no read-system-call instrumentation.'}


def audit(study, route_path, bootstrap_cache, evaluation_code=None, *, expected_packages=PACKAGES):
    started = time.time()
    study = Path(study).resolve()
    ext = study/'external_fixed_outputs_v1'
    fixed, predictions, analysis = (read(ext/name) for name in ('fixed_outputs.json', 'predictions.json', 'analysis.json'))
    rows = [json.loads(line) for line in (ext/'enriched_rows.jsonl').read_text(encoding='utf-8').splitlines()]
    evaluation = read(ext/'evaluation_design.json')
    require(analysis['evaluation_design'] == evaluation, 'Analysis/copied evaluation design differs')
    route = read(route_path)
    design, sources, indexes = selection_check(study, fixed, rows, route)
    vector_receipt = vector_check(study, sources, fixed['rows'], predictions, evaluation, analysis, expected_packages)
    row_checks = row_arithmetic(rows, predictions, read(ext/'profiles.json'), read(ext/'structure_descriptors.json'))
    arithmetic = arithmetic_check(rows, design, analysis, bootstrap_cache)
    ordering = ordering_check(study, evaluation, evaluation_code or study/'evaluate_independent.py')
    unknowns = read(ext/'execution_unknowns.json')
    require(unknowns['rows'] == [r for r in rows if not r['scientific_endpoint_available']]
            and unknowns['by_status'] == dict(Counter(r['status'] for r in unknowns['rows'])), 'Explicit unknown output membership differs')
    from scipy import __version__ as scipy_version
    return {'status': 'passed', 'scope': 'Complete saved 16-run/1024-position/5120-row endpoint/provenance/diagnostic/run-statistics audit',
        'study': str(study), 'route_audit': str(route_path), 'started_unix': started, 'ended_unix': time.time(),
        'python_executable': sys.executable, 'SciPy_reference_version': scipy_version,
        'scientific_work': {'RF_predictions': 0, 'forward_reactions': 0, 'Graph_searches': 0,
                            'Chemprop_predictions': 0, 'model_constructors_called': 0, 'descriptor_calculations': 0},
        'source_or_method_files_modified': False, 'unknown_evidence_promoted_to_success': False,
        'summary': {'raw_source_positions': 1024, 'independent_optimizer_runs': 16, 'five_policy_rows': 5120,
            'endpoint_row_checks': len(row_checks), 'paired_run_contrasts_checked': arithmetic['paired_run_contrasts_checked'],
            'technical_unknown_rows': len(unknowns['rows']), **vector_receipt},
        'ordering_evidence': ordering, 'independent_arithmetic': arithmetic, 'row_endpoint_checks': row_checks,
        'limitations': ['Released checkpoint declaration, exact saved origin vectors, compatible feature/runtime versions and conflicts are checked; no weights were repredicted or independently authenticated.',
            'Older fold files without model-reference headers require actual saved study metadata plus an exact aggregate binding; absent evidence fails the audit rather than inventing metadata.',
            'Saved chemical profiles are bound to their physical descriptor cache without recalculating descriptors. Canonical provenance keys must match saved canonical identities.',
            'The separate native route audit proves saved route/stock/decision semantics; this script does not establish laboratory synthesis or biological activity.',
            'Sixteen optimizer runs, conditional on two fixed tasks/assets, are the inference units. Repeated positions, repeated products and ten checkpoints are not independent experiments.',
            'The sole primary is B fixed minus B strict. Secondary comparisons are exploratory and cannot replace an unsupported primary.',
            'Technical unknowns retain paired identification bounds only; operational source fallbacks do not become observed final-policy properties.',
            'A within-task 100000-run bootstrap is calculated once and saved by its actual paired effects/seed; it is a sensitivity interval, not independent new scientific evidence.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--route-audit', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--evaluation-code', type=Path)
    parser.add_argument('--bootstrap-cache', type=Path)
    args = parser.parse_args()
    scope = args.study.resolve()/'witness_reaudit'
    cache = args.bootstrap_cache or scope/'independent_primary_bootstrap_20261005.json'
    require(args.report.resolve().is_relative_to(scope) and cache.resolve().is_relative_to(scope), 'Audit outputs outside scoped folder')
    try:
        receipt = audit(args.study, args.route_audit, cache, args.evaluation_code)
        code = 0
    except Exception:
        receipt = {'status': 'verification_failed_evidence_not_upgraded', 'scope': 'Complete terminal saved endpoint audit',
                   'error': traceback.format_exc(), 'scientific_model_calls': 0, 'source_or_method_files_modified': False}
        code = 1
    save(args.report, receipt)
    print(json.dumps(receipt.get('summary', receipt), indent=2, allow_nan=False), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
