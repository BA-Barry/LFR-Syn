"""Missing-only external evaluation of the new sixteen-run qualification study.

All five method outputs are saved before external histories are read or any
predictor is constructed. Retained chemistry/predictor implementations and
exact fold readers are used unchanged. No generation, guiding RF, forward
reaction or route model is called here.
"""
import argparse
from collections import Counter, defaultdict
import copy
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path
import statistics
import sys
import time

from independent_statistics import infer_from_runs

POLICIES = ('Original', 'B_joint6plus6', 'C_joint6plus6', 'B_fixedtask6plus6', 'C_fixedtask6plus6')
TASKS = ('gsk3b', 'jnk3')
ENDPOINTS = ('external_joint', 'joint', 'solved')
TERMINAL = ('complete', 'complete_with_execution_limits')
SOURCE_FIELDS = ('sample_id', 'task', 'task_index', 'replicate', 'cohort_id', 'position',
                 'canonical_smiles', 'valid', 'task_score', 'qed', 'search_seed', 'graph_verified', 'graph_result_path')
FILTERS = ('all', 'no_pains', 'no_pains_brenk', 'mw_le500', 'source_qed_tolerance')
HISTORICAL_STUDIES = ('graphga_transfer512_v1', 'graphga_fixedtask128_v1',
                      'graphga_fixedtask512_v1', 'graphga_candidateallocation128_v1',
                      'graphga_retainedprecursor128_v1')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def probability(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and 0 <= value <= 1


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dependencies(code_dir, old_study):
    sys.path.insert(0, str(Path(code_dir).resolve()))
    original = load_module('_qualification_external', Path(code_dir) / 'joint_validation512_external.py')
    adapter = load_module('_qualification_history', Path(old_study) / 'external_adapter/graphga_external.py')
    saved_analysis = load_module('_qualification_saved_arithmetic', Path(old_study) / 'analysis_adapter/analyze_rows.py')
    return original, adapter, saved_analysis


def validate_design(design):
    require(design['experiment'] == 'graphga_qualification_independent1024_v1', 'Wrong independent study')
    require(tuple(design['policies']) == POLICIES, 'Five prespecified policies required')
    require(design['raw_source_positions'] == 1024 and design['source_positions_per_run'] == 64,
            'All 1024 raw positions required')
    runs = design['runs']
    require(len(runs) == 16 and len({r['cohort_id'] for r in runs}) == 16, 'Sixteen distinct complete optimizer runs required')
    require({(r['task'], r['replicate']) for r in runs} == {(t, r) for t in TASKS for r in range(8)}, 'Eight runs per fixed task required')
    for run in runs:
        i = 2 * run['replicate'] + TASKS.index(run['task'])
        require(run['task_index'] == TASKS.index(run['task'])
                and run['cohort_id'] == 'graphga-qualconfirm-%s-rep%d' % (run['task'], run['replicate'])
                and run['generation_seed'] == 2026101501 + i
                and run['initialization_seed'] == 2026101601 + i
                and run['terminal_sampling_seed'] == 2026101701 + i
                and run['search_seed_base'] == 12500000 + 10000 * i, 'Prespecified run identity/seeds changed')
    require(design['method']['source_RF_cutoffs'] == {'gsk3b': .1, 'jnk3': .05}, 'Task RF thresholds changed')
    require(design['method']['candidate_calls'] == [6, 6]
            and design['method']['candidate_shared_soft_seconds'] == 45., 'Candidate budgets changed')
    return {run['cohort_id']: run for run in runs}


def endpoint_bounds(row, endpoint):
    value = row[endpoint]
    if type(value) is bool:
        return (int(value), int(value))
    require(value is None, 'Explicit endpoint boolean or null required')
    lo, hi = row.get('endpoint_lower_' + endpoint, 0), row.get('endpoint_upper_' + endpoint, 1)
    require(lo in (0, 1) and hi in (0, 1) and lo <= hi, 'Invalid endpoint bounds')
    return int(lo), int(hi)


def validate_row(row, source, design):
    for field, source_field in (('sample_id', 'sample_id'), ('task', 'task'), ('cohort_id', 'cohort_id'),
            ('replicate', 'replicate'), ('position', 'position'), ('source_smiles', 'canonical_smiles'),
            ('source_score', 'task_score'), ('source_qed', 'qed')):
        require(row[field] == source[source_field], 'Policy source identity changed: ' + field)
    require('source_valid' not in row or row['source_valid'] is source['valid'], 'Copied source-valid identity changed')
    retained = {key: source.get(key) for key in SOURCE_FIELDS}
    require(row['source_identity'] in (source, retained), 'Full or declared source identity required')
    require(type(source['valid']) is bool and source['valid'] == bool(source['canonical_smiles']), 'Invalid source-valid identity')
    require(probability(row['source_score']) and probability(row['final_score']), 'Invalid retained RF probability')
    require(type(row['scientific_endpoint_available']) is bool and type(row['replaced']) is bool, 'Explicit method availability required')
    require(type(row['final_structure_known']) is bool, 'Explicit scientific final-structure availability required')
    require(row['status'] not in ('processing', 'pending'), 'A live policy cannot be scored as terminal')
    for key in ('source_qed', 'final_qed'):
        require(row[key] is None or probability(row[key]), 'Invalid QED')
    for endpoint in ('solved', 'joint'):
        endpoint_bounds(row, endpoint)
        require(row['scientific_endpoint_available'] == (type(row[endpoint]) is bool), 'Unknown endpoint manufactured or lost')
    if row['scientific_endpoint_available']:
        cutoff = design['method']['source_RF_cutoffs'][row['task']]
        require(row['joint'] == bool(row['solved'] and source['valid'] and row['final_score'] >= cutoff), 'Guiding RF endpoint changed')
        require(row['final_structure_known'] == bool(row['final_smiles']), 'Known output structure identity changed')
        if row['replaced']:
            require(row['solved'] and row['final_smiles'] != row['source_smiles']
                    and row['final_score'] >= cutoff and row['final_qed'] is not None and row['source_qed'] is not None
                    and row['final_qed'] >= row['source_qed'] - .05, 'Fixed task/QED qualification violated')
            if row['policy'] in POLICIES[1:3]:
                require(row['final_score'] >= row['source_score'], 'Strict source RF guarantee violated')
    else:
        require(not row['replaced'] and not row['final_structure_known'], 'Unknown result cannot become a known replacement')


def fixed_inputs(study):
    study = Path(study)
    design = read(study / 'design.json')
    specs = validate_design(design)
    execution = read(study / 'method_execution.json')
    require(execution['status'] in TERMINAL, 'Aggregate method execution has not terminated')
    measured = read(study / 'method_result.json')
    source_record = read(study / 'sources_with_parents.json')
    sources, rows = copy.deepcopy(source_record['rows']), copy.deepcopy(measured['rows'])
    require(measured['selected_raw_n'] == 1024 and len(sources) == 1024 and len(rows) == 5120, 'Complete five-policy denominator required')
    index = {s['sample_id']: s for s in sources}
    require(len(index) == 1024, 'Duplicate raw sample primary key')
    for source in sources:
        spec = specs.get(source['cohort_id'])
        require(spec is not None and (source['task'], source['replicate'], source['task_index'])
                == (spec['task'], spec['replicate'], spec['task_index']), 'Source outside specified optimizer runs')
        require(type(source['position']) is int and source['position'] in range(64), 'Source raw draw position invalid')
        # Search seed may be adopted by exact reuse; its original request is distinct.
        requested = source.get('requested_search_seed', source['search_seed'])
        require(requested == spec['search_seed_base'] + source['position'], 'Source requested search seed changed')
    for cohort in specs:
        selected = [s for s in sources if s['cohort_id'] == cohort]
        require(len(selected) == 64 and {s['position'] for s in selected} == set(range(64)), 'Missing or duplicate raw draw position')
    keys = [(r['sample_id'], r['policy']) for r in rows]
    expected = {(sid, policy) for sid in index for policy in POLICIES}
    require(len(set(keys)) == len(keys) and set(keys) == expected, 'Missing/duplicate/extra policy membership')
    for row in rows:
        validate_row(row, index[row['sample_id']], design)
    origins = []
    for cohort in specs:
        directory = study / 'method_outputs' / cohort
        state = read(directory / 'method_execution.json')
        require(state['status'] in TERMINAL, 'Cohort execution is not terminal')
        saved_sources = read(directory / 'sources_with_parents.json')['rows']
        saved_rows = read(directory / 'method_result.json')['rows']
        actual_sources = [s for s in sources if s['cohort_id'] == cohort]
        actual_rows = [r for r in rows if r['cohort_id'] == cohort]
        require({s['sample_id']: s for s in saved_sources} == {s['sample_id']: s for s in actual_sources}
                and len(saved_sources) == 64, 'Aggregate/per-cohort source objects differ')
        require({(r['sample_id'], r['policy']): r for r in saved_rows}
                == {(r['sample_id'], r['policy']): r for r in actual_rows} and len(saved_rows) == 320,
                'Aggregate/per-cohort complete policy objects differ')
        origins.append({'cohort_id': cohort, 'directory': str(directory), 'status': state['status']})
    # Native Graph kernels retain valid inside source_identity, not a flat field.
    # Derive this redundant analysis field after comparing unchanged cohort files.
    for row in rows:
        row.setdefault('source_valid', index[row['sample_id']]['valid'])
    order = list(specs)
    sources.sort(key=lambda s: (order.index(s['cohort_id']), s['position']))
    rows.sort(key=lambda r: (order.index(r['cohort_id']), r['position'], POLICIES.index(r['policy'])))
    return design, sources, rows, measured, origins


def write_fixed(out, design, sources, rows, origins):
    data = {'experiment': design['experiment'], 'design': design, 'sources': sources, 'rows': rows,
            'row_provenance': origins, 'scope': 'All fresh raw positions/five policies saved before external history reads and predictor construction.'}
    path = Path(out) / 'fixed_outputs.json'
    require(not path.exists() or read(path) == data, 'Existing fixed outputs differ; no reselection is permitted')
    save(path, data)


def profile_sources(dev):
    paths = []
    for name in HISTORICAL_STUDIES:
        for relative in ('external_fixed_outputs_v1/profiles.json', 'candidate_profiles.json'):
            path = Path(dev) / name / relative
            if path.is_file():
                paths.append(path)
    return paths


def current_packages():
    packages = {}
    for name in ('chemprop', 'torch', 'rdkit', 'numpy'):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return packages


def validate_named_history_runtimes(dev, packages):
    """The completed main iteration caches retain actual feature runtimes."""
    evidence = []
    for name in HISTORICAL_STUDIES:
        directory = Path(dev) / name / 'external_fixed_outputs_v1'
        runtime_path = directory / 'external_runtime.json'
        if runtime_path.is_file():
            runtime = read(runtime_path)
        else:
            runtime_path = directory / 'analysis.json'
            runtime = read(runtime_path)['runtime']
        require(runtime['packages'] == packages, 'Named full-vector feature/runtime versions differ: ' + name)
        evidence.append({'study': name, 'runtime_path': str(runtime_path), 'packages': runtime['packages']})
    return evidence


def seed_descriptor_cache(out, dev, chemistry):
    path = Path(out) / 'structure_descriptors.json'
    cached = read(path) if path.exists() else {'rdkit_version': chemistry.version,
        'fingerprint': 'Morgan radius2/2048,chiral=False', 'scaffold': 'Murcko includeChirality=False; empty acyclic counted separately',
        'profiles': {}, 'source_final_similarities': []}
    require(cached['rdkit_version'] == chemistry.version, 'Current descriptor runtime differs')
    prior, origins = {}, {}
    for old_path in profile_sources(dev):
        old = read(old_path)
        require(old['rdkit_version'] == chemistry.version, 'Historical descriptor runtime differs')
        for smiles, value in old['profiles'].items():
            require(smiles not in prior or prior[smiles] == value, 'Historical descriptor conflict')
            require(smiles not in cached['profiles'] or cached['profiles'][smiles] == value, 'Current/historical descriptor conflict')
            prior.setdefault(smiles, value)
            origins.setdefault(smiles, str(old_path))
            cached['profiles'].setdefault(smiles, value)
    similarities = {(r['source'], r['final']): r['similarity'] for r in cached['source_final_similarities']}
    for name in HISTORICAL_STUDIES:
        old_path = Path(dev) / name / 'external_fixed_outputs_v1/structure_descriptors.json'
        if old_path.is_file():
            old = read(old_path)
            require(old['rdkit_version'] == chemistry.version, 'Historical similarity runtime differs')
            for row in old['source_final_similarities']:
                key = row['source'], row['final']
                require(key not in similarities or similarities[key] == row['similarity'], 'Historical similarity conflict')
                similarities.setdefault(key, row['similarity'])
    cached['source_final_similarities'] = [{'source': a, 'final': b, 'similarity': value}
        for (a, b), value in sorted(similarities.items(), key=lambda item: str(item[0]))]
    save(path, cached)
    return {'retained_profile_sources': [str(p) for p in profile_sources(dev)], 'seeded_unique_profiles': len(prior),
            'seeded_unique_similarity_pairs': len(similarities), 'profile_origins': origins}


def distribution(values):
    values = sorted(float(v) for v in values if v is not None)
    if not values:
        return {'n': 0, 'mean': None, 'median': None}
    return {'n': len(values), 'mean': statistics.fmean(values), 'median': statistics.median(values),
            'min': values[0], 'max': values[-1]}


def mask_unknown_finals(rows):
    """Keep fixed selection fields unchanged; mark main scientific fields unknown.

    The old kernel reports an operational source fallback in an unknown policy.
    Its properties are never used as known final-policy properties or success.
    """
    for row in rows:
        if row['scientific_endpoint_available']:
            continue
        row['operational_fallback_external_probability'] = row.get('external_probability')
        row['operational_fallback_descriptors'] = row.get('descriptors')
        row['operational_fallback_external_fold_probabilities'] = row.get('external_fold_probabilities', [])
        row['operational_fallback_external_score_provenance'] = row.get('external_score_provenance', [])
        row['operational_fallback_smiles'] = row.get('operational_fallback_smiles', row['final_smiles'])
        for key in ('external_probability', 'external_change', 'external_fold_sd', 'qed_delta',
                    'molecular_weight_delta', 'descriptors', 'final_scaffold', 'source_morgan_similarity'):
            row[key] = None
        row['external_fold_probabilities'] = []
        row['external_score_provenance'] = []
        row['external_valid'] = None
        row['external_score_available'] = False
        row['external_joint'] = None
        row['endpoint_lower_external_joint'] = 0
        row['endpoint_upper_external_joint'] = endpoint_bounds(row, 'solved')[1]
        row['output_structure_status'] = 'scientific_final_unknown_operational_fallback_separate'


def summarize(rows):
    n = len(rows)
    result = {'raw_n': n, 'replaced_slots': sum(r['replaced'] for r in rows),
              'route_execution_unknown_slots': sum(not r['scientific_endpoint_available'] for r in rows),
              'distinct_source_task_structures': len({(r['task'], r['source_smiles']) for r in rows if r['source_smiles']}),
              'distinct_final_task_structures_observed': len({(r['task'], r['final_smiles']) for r in rows if r['final_structure_known']})}
    result['known_route_failures_within_saved_recipe'] = sum(r['solved'] is False for r in rows)
    result['candidate_query_status_counts'] = dict(Counter(a['query']['status']
        for r in rows for a in r.get('attempts', []) if 'query' in a))
    result['candidate_budget_unclosed_attempts'] = sum(a.get('query', {}).get('status') == 'unsolved_within_budget'
        for r in rows for a in r.get('attempts', []))
    for endpoint in ENDPOINTS:
        bounds = [endpoint_bounds(r, endpoint) for r in rows]
        successes = [r for r in rows if r[endpoint] is True]
        structures = {(r['task'], r['final_smiles']) for r in successes}
        unknown_upper = sum(r[endpoint] is None and endpoint_bounds(r, endpoint)[1] for r in rows)
        result[endpoint] = {'known_successes': len(successes), 'unknown_slots': sum(r[endpoint] is None for r in rows),
            'raw_count_bounds': [sum(b[i] for b in bounds) for i in (0, 1)],
            'raw_rate_bounds': [sum(b[i] for b in bounds) / n for i in (0, 1)],
            'distinct_task_structure_known_successes': len(structures),
            'distinct_task_structure_count_bounds': [len(structures), len(structures) + unknown_upper],
            'distinct_nonempty_successful_task_scaffolds': len({(r['task'], r['final_scaffold']) for r in successes if r.get('final_scaffold')}),
            'acyclic_successful_task_structures': len({(r['task'], r['final_smiles']) for r in successes if not r.get('final_scaffold')})}
    profiles = [r['descriptors'] for r in rows if r.get('descriptors') is not None]
    result['profile_diagnostics'] = {'raw_denominator': n, 'observed_final_profile_slots': len(profiles),
        'pains_flagged_slots': sum(bool(p['pains']) for p in profiles), 'brenk_flagged_slots': sum(bool(p['brenk']) for p in profiles),
        'mw_above500_slots': sum(p['molecular_weight'] > 500 for p in profiles),
        'note': 'Alerts remain in the primary denominator; observed-profile diagnostics do not claim activity or safety.'}
    return result


def contrast(rows, design, terms, endpoint, *, primary=False, bootstrap_replicates=100000):
    indexes = {policy: {r['sample_id']: r for r in rows if r['policy'] == policy} for policy, _ in terms}
    sources = [r for r in rows if r['policy'] == 'Original']
    runs = []
    for spec in design['runs']:
        selected = [s for s in sources if s['cohort_id'] == spec['cohort_id']]
        lo = hi = unknown = wins = losses = 0
        for source in selected:
            values = [(weight, indexes[policy][source['sample_id']]) for policy, weight in terms]
            intervals = [(weight, endpoint_bounds(row, endpoint)) for weight, row in values]
            lo += sum(weight * b[0 if weight > 0 else 1] for weight, b in intervals)
            hi += sum(weight * b[1 if weight > 0 else 0] for weight, b in intervals)
            if any(row[endpoint] is None for weight, row in values):
                unknown += 1
            else:
                value = sum(weight * int(row[endpoint]) for weight, row in values)
                wins += int(value > 0)
                losses += int(value < 0)
        require(len(selected) == 64, 'Every run retains all paired raw positions')
        runs.append({'cohort_id': spec['cohort_id'], 'task': spec['task'], 'replicate': spec['replicate'],
            'raw_n': 64, 'net_count_bounds': [lo, hi], 'effect_bounds': [lo / 64, hi / 64],
            'known_positive_positions': wins, 'known_negative_positions': losses, 'unknown_pairs': unknown})
    result = {'terms': [{'policy': p, 'weight': w} for p, w in terms], 'endpoint': endpoint,
              'raw_n': len(sources), 'runs': runs,
              'scope': 'Prespecified sole primary' if primary else 'Prespecified exploratory secondary; unadjusted interval, not a substitute for the primary',
              **infer_from_runs(runs, bootstrap=primary, replicates=bootstrap_replicates)}
    if not primary and result['welch'] is not None:
        result['welch'].pop('p_two_sided', None)
        for task in result['welch']['task_summaries'].values():
            task.pop('p_two_sided', None)
    return result


def property_changes(rows):
    result = {}
    for label, selected in (('all_positions', rows), ('replacement_positions', [r for r in rows if r['replaced']])):
        n = len(selected)
        known = [r for r in selected if r['scientific_endpoint_available']]
        unknown = len(selected) - len(known)
        rf_changes = [r['final_score'] - r['source_score'] for r in known]
        decreased = sum(v < 0 for v in rf_changes)
        result[label] = {'raw_denominator': n, 'unknown_final_positions': unknown,
            'RF_change_observed': distribution(rf_changes),
            'external_probability_change_observed': distribution(r.get('external_change') for r in known),
            'QED_change_observed': distribution(r.get('qed_delta') for r in known),
            'molecular_weight_change_observed': distribution(r.get('molecular_weight_delta') for r in known),
            'RF_decreased_count_bounds': [decreased, decreased + unknown],
            'RF_decreased_fraction_bounds': [decreased / n, (decreased + unknown) / n] if n else None,
            'mean_RF_change_identification_bounds': [(sum(rf_changes) + sum(i - r['source_score'] for r in selected if not r['scientific_endpoint_available'])) / n for i in (0, 1)] if n else None,
            'mean_external_change_identification_bounds': [(sum(r.get('external_change') or 0 for r in known) + sum(i - r['source_external_probability'] for r in selected if not r['scientific_endpoint_available'])) / n for i in (0, 1)] if n else None,
            'observed_mean_note': 'Observed-property n is reported beside the unchanged raw denominator. Replacement subsets are selected outcomes, not full-position effects.'}
    result['audited_routes_below_external_cutoff'] = sum(r['solved'] is True and r['external_probability'] < .5 for r in rows)
    return result


def concentration(rows):
    successes = [r for r in rows if r['external_joint'] is True]
    products = Counter((r['task'], r['final_smiles']) for r in successes)
    replacements = Counter((r['task'], r['final_smiles']) for r in rows if r['replaced'])
    cohorts = defaultdict(set)
    for row in successes:
        cohorts[(row['task'], row['final_smiles'])].add(row['cohort_id'])
    highest = [v for key, v in products.most_common(2)]
    return {'raw_denominator': len(rows), 'known_successful_positions': len(successes),
        'distinct_successful_task_structures': len(products), 'top_one_successful_positions': highest[0] if highest else 0,
        'top_two_successful_positions': sum(highest),
        'top_one_fraction_of_raw_positions': highest[0] / len(rows) if highest else 0.,
        'top_two_fraction_of_known_successes': sum(highest) / len(successes) if successes else None,
        'distinct_replacement_task_structures': len(replacements),
        'successful_products': [{'task': k[0], 'canonical_smiles': k[1], 'positions': v, 'optimizer_runs': len(cohorts[k])} for k, v in products.most_common()],
        'products_successful_in_multiple_runs': sum(len(v) > 1 for v in cohorts.values()),
        'global_union_note': 'Descriptive structure union and shared products, not additional independent optimizer-run evidence.'}


def filter_pass(row, name):
    profile = row.get('descriptors')
    if profile is None:
        return False
    return {'all': True, 'no_pains': not profile['pains'],
        'no_pains_brenk': not profile['pains'] and not profile['brenk'], 'mw_le500': profile['molecular_weight'] <= 500,
        'source_qed_tolerance': row['source_descriptors'] is not None and profile['qed'] >= row['source_descriptors']['qed'] - .05}[name]


def structural_sensitivity(rows):
    result = {}
    for name in FILTERS:
        successes = [r for r in rows if r['external_joint'] is True and filter_pass(r, name)]
        possible = sum(endpoint_bounds(r, 'external_joint')[1] for r in rows if r['external_joint'] is None)
        result[name] = {'raw_denominator': len(rows), 'external_joint_count_bounds': [len(successes), len(successes) + possible],
            'distinct_successful_task_structures': len({(r['task'], r['final_smiles']) for r in successes}),
            'distinct_nonempty_successful_task_scaffolds': len({(r['task'], r['final_scaffold']) for r in successes if r.get('final_scaffold')}),
            'acyclic_successful_task_structures': len({(r['task'], r['final_smiles']) for r in successes if not r.get('final_scaffold')})}
    return result


def mechanism_counts(study, design):
    path = Path(study) / 'candidate_stage.json'
    require(path.is_file(), 'Complete shared candidate-stage record required')
    stage = read(path)
    require(stage['selected_raw_n'] == 1024 and len(stage['sources']) == 1024, 'Mechanism stage denominator differs')
    saved_sources = read(Path(study) / 'sources_with_parents.json')['rows']
    enrolled = {source['sample_id']: source for source in saved_sources}
    observed = {entry['source']['sample_id']: entry['source'] for entry in stage['sources']}
    require(len(enrolled) == len(observed) == 1024 and observed == enrolled,
            'Mechanism pools do not belong to the complete saved source identities')
    by_task = {}
    for task in TASKS:
        entries = [s for s in stage['sources'] if s['source']['task'] == task]
        require(len(entries) == 512, 'Mechanism task membership differs')
        current = {'raw_positions': len(entries), 'source_RF_exact_one_slots': sum(s['source']['task_score'] == 1 for s in entries),
            'source_RF_distribution': distribution(s['source']['task_score'] for s in entries),
            'triggered_positions': sum(bool(s['source']['valid']) and not (s['source'].get('graph_verified') and s['source']['task_score'] >= design['method']['source_RF_cutoffs'][task]) for s in entries),
            'arms': {}}
        for arm, labels in (('B', ('task_contribution', 'B')), ('C', ('permuted_contribution', 'C'))):
            totals = Counter()
            for entry in entries:
                matching = [entry.get('arms', {}).get(label) for label in labels if label in entry.get('arms', {})]
                require(len(matching) <= 1, 'Duplicate chemical arm labels')
                if not matching:
                    continue
                record = matching[0]
                incomplete = record.get('native_status') not in (None, 'complete')
                if incomplete:
                    totals['technical_unknown_pool_positions'] += 1
                    continue
                candidates = record.get('scored_candidates', [])
                totals['observed_pool_positions'] += 1
                totals['formed_candidate_positions'] += bool(candidates)
                totals['logical_score_slots'] += len(record.get('planned_candidates', []))
                fixed = [c for c in candidates if c['task_score'] >= design['method']['source_RF_cutoffs'][task]]
                strict = [c for c in fixed if c['task_score'] >= entry['source']['task_score']]
                totals['fixed_RF_qualified_positions'] += bool(fixed)
                totals['source_relative_RF_qualified_positions'] += bool(strict)
                totals['fixed_top_two_slots'] += min(2, len(fixed))
                totals['strict_top_two_slots'] += min(2, len(strict))
                totals['fixed_two_candidate_positions'] += len(fixed) >= 2
                totals['strict_two_candidate_positions'] += len(strict) >= 2
            current['arms'][arm] = dict(totals)
        by_task[task] = current
    return {'by_task': by_task, 'selection_use': 'None: diagnostics of the shared, already saved candidate formation and scoring pool.',
            'unknown_pool_note': 'Incomplete pool observations are counted separately, not classified as chemical coverage failures.'}


def analyze(design, rows, study, *, bootstrap_replicates=100000):
    arms = {p: [r for r in rows if r['policy'] == p] for p in POLICIES}
    primary = contrast(rows, design, ((POLICIES[3], 1), (POLICIES[1], -1)), 'external_joint', primary=True,
                       bootstrap_replicates=bootstrap_replicates)
    pairs = ((POLICIES[4], POLICIES[2]), (POLICIES[1], POLICIES[2]), (POLICIES[3], POLICIES[4]),
             (POLICIES[3], 'Original'), (POLICIES[4], 'Original'))
    secondary = [contrast(rows, design, ((a, 1), (b, -1)), endpoint) for a, b in pairs for endpoint in ENDPOINTS]
    secondary.extend(contrast(rows, design, ((POLICIES[3], 1), (POLICIES[1], -1)), e) for e in ('joint', 'solved'))
    interaction = ((POLICIES[3], 1), (POLICIES[4], -1), (POLICIES[1], -1), (POLICIES[2], 1))
    secondary.extend(contrast(rows, design, interaction, e) for e in ENDPOINTS)
    inputs = arms['Original']
    union = {(r['task'], r['source_smiles']) for r in inputs if r['source_smiles']}
    counts = Counter((r['task'], r['source_smiles']) for r in inputs if r['source_smiles'])
    return {'experiment': design['experiment'], 'raw_source_positions': 1024, 'policy_rows': 5120,
        'independent_optimizer_runs': 16, 'independent_runs_per_fixed_task': 8,
        'primary_contrast': design['primary_contrast'], 'primary': primary, 'exploratory_secondary_contrasts': secondary,
        'policies': {p: summarize(arms[p]) for p in POLICIES},
        'by_cohort': {r['cohort_id']: {p: summarize([x for x in arms[p] if x['cohort_id'] == r['cohort_id']]) for p in POLICIES} for r in design['runs']},
        'by_task': {t: {p: summarize([x for x in arms[p] if x['task'] == t]) for p in POLICIES} for t in TASKS},
        'property_changes': {p: property_changes(arms[p]) for p in POLICIES},
        'concentration': {p: concentration(arms[p]) for p in POLICIES},
        'structural_sensitivity': {p: structural_sensitivity(arms[p]) for p in POLICIES},
        'shared_candidate_mechanisms': mechanism_counts(study, design),
        'input_duplicates': {'raw_positions': 1024, 'valid_source_positions': sum(r['source_valid'] for r in inputs),
            'distinct_task_structures': len(union), 'repeated_valid_structure_positions': sum(counts.values()) - len(union)},
        'limits': ['The sixteen new optimizer runs are the statistical units. Positions, duplicate molecules and ten checkpoints are not independent repetitions.',
            'Old development optimizer runs are excluded from all new primary run-level inference.',
            'Released Chemprop checkpoints and optimized RF share the broad ExCAPE source; incomplete original training composition prevents a claim of independent biological validation.',
            'Native-route audit is inherited from saved method records and independently rechecked after scientific execution.',
            'Alerts and low properties remain in the primary denominator. Secondary contrasts are exploratory; no positive secondary replaces the sole primary.',
            'Technical unknowns retain identifiable paired bounds. Source fallback properties do not become known final-policy properties.']}


def execute(args, *, dependencies_override=None, chemistry=None, predictor_factory=None, bootstrap_replicates=100000, runtime_override=None):
    started = time.monotonic()
    original, adapter, saved_analysis = dependencies_override or dependencies(args.code_dir, args.old_study)
    original.configure_cpu()
    chemistry = chemistry or original.Chemistry()
    design, sources, rows, method, origins = fixed_inputs(args.study)
    for row in rows:
        for key in ('source_smiles', 'final_smiles'):
            require(not row[key] or chemistry.canonical(row[key]) == row[key], 'Saved canonical structure identity changed')
    out = Path(args.study) / 'external_fixed_outputs_v1'
    out.mkdir(parents=True, exist_ok=True)
    # The entire fresh study's choices precede ALL external-cache reads below.
    write_fixed(out, design, sources, rows, origins)
    packages = runtime_override if runtime_override is not None else current_packages()
    runtime_evidence = validate_named_history_runtimes(args.dev, packages)
    required, pairs = original.required_structures(sources, rows)
    groups, reuse_note = adapter.history_groups(args.dev, original)
    groups[0:0] = [(name, [Path(args.dev) / name / 'external_fixed_outputs_v1/predictions.json']) for name in HISTORICAL_STUDIES]
    for label, paths in groups:
        require(all(path.is_file() for path in paths), 'Recover named historical vectors instead of repeating inference: ' + label)
    history, history_origins, history_audit = adapter.merge_history(groups, original, chemistry.canonical, required)
    coverage = {}
    for fold in range(10):
        _, current = original.new_fold(out, fold)
        coverage[str(fold)] = {'historical': sum(s in history[fold] for s in required),
            'current': sum(s not in history[fold] and s in current['predictions'] for s in required),
            'missing': [s for s in required if s not in history[fold] and s not in current['predictions']]}
    save(out / 'prediction_plan.json', {'required_smiles': required, 'task_structure_pairs': pairs,
        'fold_coverage': coverage, 'historical_reuse': history_audit, 'model_reference': original.MODEL_REFERENCE})
    evaluation = {'experiment': design['experiment'], 'raw_positions': 1024, 'optimizer_runs': 16,
        'policies': list(POLICIES), 'fixed_outputs_no_reselection': True, 'primary_threshold': .5,
        'external_scores_used_in_selection': False, 'model_reference': original.MODEL_REFERENCE,
        'history_groups': [{'reuse_source': label, 'paths': [str(p) for p in paths]} for label, paths in groups],
        'statistics': design['statistics'], 'CPU_only_external_inference': True,
        'named_history_feature_runtime_evidence': runtime_evidence, **reuse_note}
    save(out / 'evaluation_design.json', evaluation)
    if args.mode == 'prepare':
        result = {'raw_source_positions': 1024, 'policy_rows': 5120, 'required_structures': len(required),
                  'missing_fold_vectors': sum(len(c['missing']) for c in coverage.values()), 'model_predictions_this_invocation': 0}
        save(out / 'preparation.json', result)
        return result
    values, inference = original.obtain_predictions(required, history, history_origins, args.assets, out,
        args.mode == 'run', predictor_factory or original.PublishedPredictor)
    prediction_records = {'tasks': list(original.TASKS), 'models': 10,
                          'model_reference': original.MODEL_REFERENCE, **inference, 'predictions': values}
    save(out / 'predictions.json', prediction_records)
    descriptor_reuse = seed_descriptor_cache(out, args.dev, chemistry)
    profiles, similarities, descriptor_work = original.descriptor_records(required, rows, chemistry, out)
    descriptor_enriched = original.enrich(rows, values, profiles, similarities, None)
    enriched, vector_audit = saved_analysis.enrich(descriptor_enriched, prediction_records)
    require(not vector_audit['missing_required_structures'], 'Selected external vectors are incomplete')
    mask_unknown_finals(enriched)
    result = analyze(design, enriched, args.study, bootstrap_replicates=bootstrap_replicates)
    result.update(external_vector_arithmetic=vector_audit, evaluation_design=evaluation,
        execution_and_reuse={'method_new_work': method.get('actual_new_work'), 'historical_predictions': history_audit,
            'inference_this_invocation': inference, 'physical_new_inference_cumulative': adapter.stored_inference_cost(out, original),
            'descriptor_reuse': descriptor_reuse, 'descriptor_work_this_invocation': descriptor_work,
            'training_generation_RF_forward_parent_work_this_external_stage': 0, 'seconds_this_invocation': time.monotonic() - started})
    result['runtime'] = {'python': sys.version, 'executable': sys.executable, 'packages': packages}
    adapter.write_outputs(out, result, enriched, saved_analysis)
    save(out / 'profiles.json', {'rdkit_version': chemistry.version, 'profiles': profiles})
    return {key: result[key] for key in ('raw_source_positions', 'policy_rows', 'primary', 'execution_and_reuse')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('prepare', 'run', 'analyze'), required=True)
    for name in ('study', 'code-dir', 'old-study', 'dev', 'assets'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(execute(args), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
