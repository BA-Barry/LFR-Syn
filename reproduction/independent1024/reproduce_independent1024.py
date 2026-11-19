"""Rebuild the completed five-arm 16-run analysis with Python's standard library."""
import argparse
import copy
import json
import math
from pathlib import Path
import platform
import statistics
import sys
import time

# Isolated Python omits the script directory. Import only the explicitly
# delivered project modules next to this entry point, never original paths.
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import evaluate_independent as study_analysis
import graphga_saved_vector_analysis as vectors_analysis
import joint_validation512_external as descriptor_arithmetic
import audit_saved_independent_endpoints as endpoint_checks


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')


def load_inputs(data):
    """Use the actual evaluator's source/row rules, preserving all raw IDs."""
    data = Path(data)
    design = read(data / 'design.json')
    specs = study_analysis.validate_design(design)
    sources = read(data / 'sources_with_parents.json')['rows']
    rows = read(data / 'policy_rows.json')['rows']
    if len(sources) != 1024 or len(rows) != 5120:
        raise ValueError('Complete 1,024-source/five-arm denominator required')
    source_index = {s['sample_id']: s for s in sources}
    if len(source_index) != len(sources):
        raise ValueError('Duplicate raw source ID')
    for source in sources:
        spec = specs.get(source['cohort_id'])
        if spec is None or (source['task'], source['replicate'], source['task_index']) != (
                spec['task'], spec['replicate'], spec['task_index']):
            raise ValueError('Source outside the prescribed sixteen optimizer runs')
        if type(source['position']) is not int or source['position'] not in range(64):
            raise ValueError('Invalid raw draw position')
        requested = source.get('requested_search_seed', source['search_seed'])
        if requested != spec['search_seed_base'] + source['position']:
            raise ValueError('Source requested search seed changed')
    for cohort in specs:
        selected = [s for s in sources if s['cohort_id'] == cohort]
        if len(selected) != 64 or {s['position'] for s in selected} != set(range(64)):
            raise ValueError('Missing or duplicate raw draw position')
    keys = [(r['sample_id'], r['policy']) for r in rows]
    expected = {(sid, policy) for sid in source_index for policy in study_analysis.POLICIES}
    if len(set(keys)) != len(keys) or set(keys) != expected:
        raise ValueError('Missing, duplicate or extra policy membership')
    for row in rows:
        study_analysis.validate_row(row, source_index[row['sample_id']], design)
        row.setdefault('source_valid', source_index[row['sample_id']]['valid'])
    order = list(specs)
    rows.sort(key=lambda r: (order.index(r['cohort_id']), r['position'], study_analysis.POLICIES.index(r['policy'])))
    return design, sources, rows


def rebuild_rows(data, rows):
    external = read(Path(data) / 'predictions.json')
    # Validates every ten-by-two record and reported means; then recomputes
    # means and checkpoint dispersion, rather than trusting derived scalars.
    vectors = vectors_analysis.vector_records(external)
    for record in vectors.values():
        if len(record['fold_provenance']) != 10:
            raise ValueError('Complete checkpoint provenance required')
        record['fold_sd'] = [statistics.pstdev(v[i] for v in record['folds']) for i in range(2)]
    rebuilt_external = {key: external[key] for key in ('tasks', 'models', 'model_reference')}
    rebuilt_external['predictions'] = vectors
    profiles = read(Path(data) / 'profiles.json')
    descriptors = read(Path(data) / 'structure_descriptors.json')
    if profiles['rdkit_version'] != descriptors['rdkit_version'] or profiles['profiles'] != descriptors['profiles']:
        raise ValueError('Saved profile/cache objects differ')
    similarities = {(r['source'], r['final']): r['similarity'] for r in descriptors['source_final_similarities']}
    if len(similarities) != len(descriptors['source_final_similarities']):
        raise ValueError('Duplicate saved source/final similarity pair')
    required = {r[k] for r in rows for k in ('source_smiles', 'final_smiles') if r[k]}
    if not required <= set(vectors):
        raise ValueError('A fixed source/final structure has no complete checkpoint vector')
    if any((r['source_smiles'], r['final_smiles']) not in similarities for r in rows):
        raise ValueError('A fixed source/final pair has no saved similarity')
    enriched = descriptor_arithmetic.enrich(rows, vectors, profiles['profiles'], similarities, None)
    enriched, vector_report = vectors_analysis.enrich(enriched, rebuilt_external)
    if vector_report['missing_required_structures']:
        raise ValueError('Complete selected checkpoint vectors required')
    study_analysis.mask_unknown_finals(enriched)
    # Independent saved-row arithmetic checks add no RDKit, SciPy or model use.
    row_checks = endpoint_checks.row_arithmetic(enriched, rebuilt_external, profiles, descriptors)
    return enriched, vector_report, len(row_checks)


def comparison(actual, expected, path='analysis', differences=None, failures=None, checked=None):
    """Record every floating difference; memberships/counts/text stay exact."""
    differences = [] if differences is None else differences
    failures = [] if failures is None else failures
    checked = [0] if checked is None else checked
    if isinstance(actual, dict):
        if not isinstance(expected, dict):
            failures.append({'path': path, 'reason': 'dictionary type differs'})
        else:
            for key, value in actual.items():
                if key not in expected:
                    failures.append({'path': path + '.' + key, 'reason': 'missing reference field'})
                else:
                    comparison(value, expected[key], path + '.' + key, differences, failures, checked)
    elif isinstance(actual, list):
        if not isinstance(expected, list) or len(actual) != len(expected):
            failures.append({'path': path, 'reason': 'list membership/length differs'})
        else:
            for i, (value, reference) in enumerate(zip(actual, expected)):
                comparison(value, reference, path + '[%d]' % i, differences, failures, checked)
    else:
        checked[0] += 1
        if type(actual) is float and type(expected) in (float, int):
            delta = actual - expected
            if delta != 0.:
                differences.append({'path': path, 'reconstructed': actual, 'reference': expected,
                                    'difference': delta, 'absolute_difference': abs(delta)})
            if not math.isfinite(actual) or not math.isfinite(expected) or abs(delta) > 1e-10:
                failures.append({'path': path, 'reason': 'numeric difference exceeds 1e-10'})
        elif type(actual) is not type(expected) or actual != expected:
            failures.append({'path': path, 'reason': 'exact type/value differs',
                             'reconstructed': actual, 'reference': expected})
    return {'scalar_fields_checked': checked[0], 'floating_differences': differences,
            'maximum_absolute_floating_difference': max((d['absolute_difference'] for d in differences), default=0.),
            'failures': failures}


def run(data, output, reference):
    started = time.monotonic()
    data, output, reference = map(lambda p: Path(p).resolve(), (data, output, reference))
    design, sources, rows = load_inputs(data)
    enriched, vector_report, checked_rows = rebuild_rows(data, rows)
    result = study_analysis.analyze(design, enriched, data, bootstrap_replicates=100000)
    result['external_vector_arithmetic'] = vector_report
    # The reference is first opened after the saved-data reconstruction.
    expected = read(reference / 'analysis.json')
    compared = comparison(result, expected)
    save(output / 'analysis.json', result)
    save(output / 'comparison.json', compared)
    policy_names = study_analysis.POLICIES
    csv = 'policy,raw_n,route_success,RF_route_success,external_route_success,distinct_external_successful_task_structures\n'
    for policy in policy_names:
        summary = result['policies'][policy]
        csv += ','.join(map(str, (policy, summary['raw_n'], summary['solved']['known_successes'],
                         summary['joint']['known_successes'], summary['external_joint']['known_successes'],
                         summary['external_joint']['distinct_task_structure_known_successes']))) + '\n'
    (output / 'policy_yields.csv').write_text(csv, encoding='utf-8')
    scientific_modules = [name for name in sys.modules if name.split('.')[0] in
                          ('rdkit', 'numpy', 'torch', 'sklearn', 'chemprop', 'syntheseus', 'scipy', 'tdc')]
    if scientific_modules:
        raise RuntimeError('Scientific third-party module imported by saved-data analysis')
    receipt = {'scope': 'actual_private_saved_data_analysis_reproduction', 'data_directory': str(data),
        'output_directory': str(output), 'reference_directory': str(reference), 'working_directory': str(Path.cwd()),
        'entry_point': str(Path(__file__).resolve()), 'python': platform.python_version(), 'python_executable': sys.executable,
        'isolated': sys.flags.isolated, 'no_site': sys.flags.no_site, 'source_positions': len(sources),
        'policy_rows': len(rows), 'independent_optimizer_runs': result['independent_optimizer_runs'],
        'complete_checkpoint_vectors': vector_report['supplied_complete_vectors'],
        'checkpoint_predictions_per_structure': 10, 'independent_row_arithmetic_checks': checked_rows,
        'primary_bootstrap_replicates': 100000, 'primary_bootstrap_seed': 2026101801,
        'comparison_scalar_fields': compared['scalar_fields_checked'], 'comparison_failures': len(compared['failures']),
        'floating_difference_count': len(compared['floating_differences']),
        'maximum_absolute_floating_difference': compared['maximum_absolute_floating_difference'],
        'imported_scientific_third_party_modules': scientific_modules,
        'new_scientific_work': {'RF_predictions': 0, 'Graph_model_calls': 0, 'RunReactants': 0,
                                'external_predictions': 0, 'source_generation': 0},
        'seconds': time.monotonic() - started,
        'status': 'reconstructed_and_compared' if not compared['failures'] else 'comparison_failed',
        'interpretation': 'Saved-data arithmetic and 16-run statistics only; original molecular identities, scaffolds, profiles and audited routes are retained measurements, not newly computed chemistry or native-witness audits.'}
    save(output / 'execution_receipt.json', receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.data, args.output, args.reference)
    print(json.dumps(result, indent=2, ensure_ascii=True))
    return int(result['comparison_failures'] != 0)


if __name__ == '__main__':
    raise SystemExit(main())
