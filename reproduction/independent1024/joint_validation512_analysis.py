"""Paired run-level analysis of the new fixed-slot validation, without inference.

All raw slots remain in the denominator. The primary interval resamples the
four complete training runs within each fixed task, not individual molecules.
Unknown execution endpoints retain bounds instead of being scored as failures.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import statistics

from joint_validation512_graph import POLICIES, THRESHOLDS, identity

CONTRASTS = [
    ('primary_Mix_minus_B_split', 'Mix_joint6plus6', 'B_joint6plus6'),
    ('attribution_B_minus_C', 'B_joint6plus6', 'C_joint6plus6'),
    ('trigger_joint_minus_route', 'B_joint12', 'B_route12'),
    ('allocation_split_minus_twelve', 'B_joint6plus6', 'B_joint12'),
    ('B_minus_Morgan', 'B_joint6plus6', 'A_joint6plus6'),
    ('Mix_minus_Morgan', 'Mix_joint6plus6', 'A_joint6plus6'),
    ('B_minus_native', 'B_joint6plus6', 'SynTwins_joint6plus6'),
    ('Mix_minus_native', 'Mix_joint6plus6', 'SynTwins_joint6plus6'),
]


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def limits(row, endpoint):
    value = row[endpoint]
    if row['scientific_endpoint_available']:
        if not isinstance(value, bool):
            raise ValueError('Available endpoint must be an explicit boolean')
        return int(value), int(value)
    if value is not None:
        raise ValueError('Unavailable scientific endpoint must remain unknown')
    lower = int(bool(row.get('endpoint_lower_' + endpoint, False)))
    upper = int(bool(row.get('endpoint_upper_' + endpoint, True)))
    if lower > upper:
        raise ValueError('Invalid endpoint bounds')
    return lower, upper


def percentile(values, probability):
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    low = int(position)
    fraction = position - low
    return ordered[low] * (1 - fraction) + ordered[min(low + 1, len(ordered) - 1)] * fraction


def exact_stratified_bootstrap(run_effects):
    """Enumerate the ordinary empirical bootstrap, conditional on two tasks.

    There are 4**4 ordered resamples per task and 65,536 paired combinations.
    This is not an exact randomization test or a population guarantee.
    """
    tasks = sorted(run_effects)
    if len(tasks) != 2 or any(len(run_effects[t]) != 4 for t in tasks):
        raise ValueError('Expected four complete training-run effects in each of two fixed tasks')
    means = []
    for task in tasks:
        values = run_effects[task]
        means.append([statistics.mean(sample) for sample in itertools.product(values, repeat=4)])
    distribution = [(a + b) / 2 for a in means[0] for b in means[1]]
    return {'method': 'Empirical percentile bootstrap of paired complete training runs, '
                      'stratified within the two fixed tasks; all ordered resamples enumerated.',
            'resample_combinations': len(distribution),
            'effect': statistics.mean(statistics.mean(run_effects[t]) for t in tasks),
            'interval_95': [percentile(distribution, .025), percentile(distribution, .975)],
            'scope': 'Conditional on these two task predictors and project assets; only four runs per task. '
                     'Descriptive uncertainty, not a power guarantee or biological validation.'}


def stratum(source):
    if not source['valid']:
        return 'invalid'
    high = source['task_score'] >= THRESHOLDS[source['task']]
    return ('high_task' if high else 'low_task') + ('_route_solved' if source['graph_verified'] else '_route_unresolved')


def policy_summary(selected):
    result = {'raw_n': len(selected), 'execution_or_pending_slots': sum(
              not row['scientific_endpoint_available'] for row in selected),
              'statuses': dict(Counter(row['status'] for row in selected)),
              'replaced': sum(row['replaced'] for row in selected),
              'logical_model_calls': sum(row['logical_model_calls'] for row in selected),
              'logical_search_seconds': sum(row['logical_search_seconds'] for row in selected),
              'unique_joint_task_structures_known': len({(row['task'], row['final_smiles'])
                                                        for row in selected if row['joint'] is True})}
    for endpoint in ('joint', 'solved'):
        bounded = [limits(row, endpoint) for row in selected]
        result[endpoint + '_bounds'] = [sum(pair[0] for pair in bounded), sum(pair[1] for pair in bounded)]
        result[endpoint + '_known_successes'] = sum(row[endpoint] is True for row in selected)
    return result


def contrast(name, first, second, cohorts, rows):
    runs, wins, losses, unobserved = [], 0, 0, 0
    complete = True
    effects = {task: [] for task in THRESHOLDS}
    for spec in cohorts:
        cohort_id = f"{spec['task']}-rep{spec['replicate']}"
        left = sorted((r for r in rows if r['policy'] == first and r['cohort_id'] == cohort_id),
                      key=lambda r: r['position'])
        right = sorted((r for r in rows if r['policy'] == second and r['cohort_id'] == cohort_id),
                       key=lambda r: r['position'])
        if len(left) != 64 or len(right) != 64:
            raise ValueError('Contrast must include every raw slot of every planned run')
        low = high = run_wins = run_losses = run_unknown = 0
        for a, b in zip(left, right):
            if a['sample_id'] != b['sample_id']:
                raise ValueError('Paired source identities differ')
            al, au = limits(a, 'joint')
            bl, bu = limits(b, 'joint')
            low += al - bu
            high += au - bl
            if a['scientific_endpoint_available'] and b['scientific_endpoint_available']:
                run_wins += int(a['joint'] and not b['joint'])
                run_losses += int(b['joint'] and not a['joint'])
            else:
                run_unknown += 1
        if run_unknown:
            complete = False
        else:
            effects[spec['task']].append(low / 64)
        runs.append({'cohort_id': cohort_id, 'task': spec['task'], 'replicate': spec['replicate'],
                     'raw_n': 64, 'first_joint_bounds': policy_summary(left)['joint_bounds'],
                     'second_joint_bounds': policy_summary(right)['joint_bounds'],
                     'net_joint_bounds': [low, high], 'known_wins': run_wins,
                     'known_losses': run_losses, 'unknown_pairs': run_unknown})
        wins += run_wins
        losses += run_losses
        unobserved += run_unknown
    lower, upper = [sum(r['net_joint_bounds'][i] for r in runs) for i in range(2)]
    result = {'name': name, 'first': first, 'second': second, 'runs': runs,
              'raw_n': 512, 'known_wins': wins, 'known_losses': losses, 'unknown_pairs': unobserved,
              'net_joint_bounds': [lower, upper], 'rate_difference_bounds': [lower / 512, upper / 512],
              'all_run_endpoints_available': complete}
    if complete:
        result['run_level_uncertainty'] = exact_stratified_bootstrap(effects)
        result['positive_runs'] = sum(r['net_joint_bounds'][0] > 0 for r in runs)
        result['negative_runs'] = sum(r['net_joint_bounds'][0] < 0 for r in runs)
        result['tied_runs'] = sum(r['net_joint_bounds'][0] == 0 for r in runs)
    return result


def analyze(study, results_paths):
    design, sources = read(study / 'design.json'), read(study / 'sources.json')['rows']
    if design['raw_denominator'] != 512 or len(sources) != 512 or len(design['runs']) != 8:
        raise ValueError('This analysis is for the prespecified eight-run raw512 study')
    source_map = {s['sample_id']: s for s in sources}
    reports = [read(path) for path in results_paths]
    rows = [row for report in reports for row in report['rows']]
    expected = {(s['sample_id'], p) for s in sources for p in POLICIES}
    actual = [(r['sample_id'], r['policy']) for r in rows]
    if set(actual) != expected or len(actual) != len(expected):
        raise ValueError('All512 source slots and all policies must appear exactly once; do not pool selected runs')
    for row in rows:
        source = source_map[row['sample_id']]
        if row['source_identity'] != identity(source):
            raise ValueError('Outcome row belongs to a different original source')
        if row['logical_model_calls'] > 12 or row['logical_model_calls'] < 0:
            raise ValueError('Policy combined call budget invalid')
    policies = {p: policy_summary([r for r in rows if r['policy'] == p]) for p in POLICIES}
    stratified = {}
    for group in sorted({stratum(s) for s in sources}):
        selected_ids = {s['sample_id'] for s in sources if stratum(s) == group}
        stratified[group] = {p: policy_summary([r for r in rows if r['policy'] == p and
                                               r['sample_id'] in selected_ids]) for p in POLICIES}
    task_summaries = {task: {p: policy_summary([r for r in rows if r['task'] == task and r['policy'] == p])
                            for p in POLICIES} for task in THRESHOLDS}
    comparisons = [contrast(name, a, b, design['runs'], rows) for name, a, b in CONTRASTS]
    return {'scope': design['scope'], 'raw_denominator': 512, 'source_training_runs': 8,
            'primary_contrast': comparisons[0], 'secondary_contrasts': comparisons[1:],
            'policies': policies, 'by_task': task_summaries, 'by_source_stratum': stratified,
            'source_counts': {'valid': sum(s['valid'] for s in sources),
                              'original_joint': sum(s['graph_verified'] and s['task_score'] >=
                                                    THRESHOLDS[s['task']] for s in sources),
                              'strata': dict(Counter(stratum(s) for s in sources))},
            'input_reports': [str(p) for p in results_paths],
            'actual_execution_by_worker': [{key: report.get(key) for key in (
                'selected_complete_cohorts', 'stored_actual_new_calls_cumulative',
                'stored_historical12_queries', 'stored_route_audit_seconds_cumulative',
                'stored_saving_seconds_cumulative', 'native_pending')} for report in reports],
            'interpretation': 'Eight training runs are the replication units. No individual-slot '
                              'independence assumption or effect-driven enrollment change. '
                              'The old256 is development evidence, not pooled into this new validation. '
                              'Native pending or execution-limited endpoints remain bounded unknowns.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=Path(
        '/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development/joint_validation512_v1'))
    parser.add_argument('--results', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = analyze(args.study, args.results)
    save(args.output or args.study / 'analysis.json', result)
    print(json.dumps({'raw_denominator': result['raw_denominator'], 'policies': result['policies'],
                      'primary_contrast': {k: v for k, v in result['primary_contrast'].items() if k != 'runs'}}, indent=2))


if __name__ == '__main__':
    main()
