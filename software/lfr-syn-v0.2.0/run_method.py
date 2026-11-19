"""Execute fixed-task B/C qualification or explicitly replay historical strict decisions."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def normalize_cli_paths(args):
    """Resolve caller-relative paths before workers change their directory.

    A bare interpreter name uses PATH. Keep interpreter arguments as strings
    until this point so ``./python`` remains distinguishable from ``python``.
    """
    for name in ('sources', 'output', 'cached_records', 'cached_outcomes',
                 'rf_model_dir', 'forward_templates', 'retro_templates',
                 'stock_smiles', 'stock_keys', 'graph_model_dir', 'reuse_study'):
        value = getattr(args, name, None)
        if value is not None:
            setattr(args, name, Path(value).resolve())
    for name in ('rf_python', 'graph_python'):
        value = getattr(args, name, None)
        if value is None:
            continue
        text = str(value)
        if not Path(text).is_absolute() and '/' not in text and '\\' not in text:
            found = shutil.which(text)
            if found is None:
                raise FileNotFoundError(f'Interpreter executable is not on PATH: {text}')
            text = found
        # A venv's Python commonly links to its base interpreter. Keep the
        # executable path so Python can discover that venv's pyvenv.cfg.
        setattr(args, name, Path(text).absolute())
    return args


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def selected_inputs(args):
    if args.smiles is not None:
        if args.task is None:
            raise ValueError('--smiles requires --task gsk3b or jnk3')
        rows = [{'sample_id': 'provided-source', 'raw_smiles': args.smiles,
                 'task': args.task, 'replicate': args.replicate, 'position': args.position}]
        if args.search_seed is not None:
            rows[0]['search_seed'] = args.search_seed
    else:
        path = args.sources or ROOT / 'examples/source_first_raw.json'
        data = read(path)
        rows = data['rows'] if isinstance(data, dict) and 'rows' in data else data
        if isinstance(rows, dict):
            rows = [rows]
        if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
            raise ValueError('Source JSON must be one object, a list, or an object with rows')
    if args.limit < 1:
        raise ValueError('--limit must be positive')
    if not rows:
        raise ValueError('Source input has no rows')
    chosen = rows[:args.limit]
    if len({r.get('sample_id', f'input-{i}') for i, r in enumerate(chosen)}) != len(chosen):
        raise ValueError('Selected source IDs must be unique')
    return {'scope': 'software_method_reproduction_not_new_scientific_validation',
            'rows': chosen, 'input_row_count': len(rows), 'selected_row_count': len(chosen),
            'selection': f'First {args.limit} supplied rows, with invalids and duplicates retained; no outcome selection by CLI.',
            'training': 0, 'sampling': 0}


def environment():
    return dict(os.environ, PYTHONPATH=str(ROOT) + os.pathsep + os.environ.get('PYTHONPATH', ''),
                OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2',
                PYTHONUTF8='1')


def execute(args):
    from lfr_method.qualification import POLICIES_BY_QUALIFICATION
    args.output.mkdir(parents=True, exist_ok=True)
    selected = selected_inputs(args)
    input_path = args.output / 'selected_inputs.json'
    if input_path.exists() and read(input_path) != selected:
        raise ValueError('Existing output contains different selected source inputs; use another output directory')
    save(input_path, selected)
    config = {'scope': selected['scope'], 'package_version': '0.2.0',
              'rf_model_dir': str(args.rf_model_dir), 'forward_templates': str(args.forward_templates),
              'retro_templates': str(args.retro_templates), 'stock_smiles': str(args.stock_smiles),
              'stock_keys': str(args.stock_keys), 'graph_model_dir': str(args.graph_model_dir),
              'rf_python': str(args.rf_python), 'graph_python': str(args.graph_python),
              'reuse_study': str(args.reuse_study) if args.reuse_study else None,
              'device': args.device, 'parent_calls': 12, 'parent_soft_seconds': 45.,
              'candidate_calls': [6, 6], 'candidate_shared_soft_seconds': 45.,
              'candidate_score_slots_per_arm': 8, 'candidate_forward_attempts_per_arm': 16,
              'qualification': args.qualification,
              'policies': list(POLICIES_BY_QUALIFICATION[args.qualification]),
              'training': 0, 'sampling': 0}
    cfg = args.output / 'method_config.json'
    if cfg.exists() and read(cfg) != config:
        raise ValueError('Output configuration differs; retained output is not overwritten')
    save(cfg, config)
    prepare = [str(args.rf_python), '-m', 'lfr_method.rf_worker', 'prepare',
               '--inputs', str(input_path), '--output', str(args.output),
               '--rf-model-dir', str(args.rf_model_dir)]
    run = [str(args.graph_python), '-m', 'lfr_method.graph_worker',
           '--output', str(args.output)]
    for name, command in [('prepare_sources', prepare), ('method_pipeline', run)]:
        log = args.output / f'{name}.log'
        print(json.dumps({'stage': name, 'log': str(log), 'command': command}), flush=True)
        with log.open('a', encoding='utf-8') as stream:
            code = subprocess.call(command, cwd=ROOT, env=environment(), stdout=stream,
                                   stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
        if code:
            raise RuntimeError(f'{name} exited {code}; saved observations retained. Inspect {log}')
    result = read(args.output / 'method_result.json')
    print(json.dumps({'output': str(args.output), 'mode': result['mode'],
                      'policies': result['policies'], 'actual_new_work': result['actual_new_work']}, indent=2))
    return 0


def replay(args):
    # This path imports only the standard-library policy evaluator. It never
    # imports RF, RDKit, Graph2Edits, or third-party checkpoints.
    from lfr_method.graph_core import summary
    from lfr_method.qualification import POLICIES_BY_QUALIFICATION, evaluator, report_qualification
    policies = POLICIES_BY_QUALIFICATION[args.qualification]
    evaluate_policy = evaluator(args.qualification)
    selected = selected_inputs(args)
    records = read(args.cached_records)
    by_id = {r['source']['sample_id']: r for r in records['sources']}
    old_rows = read(args.cached_outcomes)['rows']
    old = {(r['sample_id'], r['policy']): r for r in old_rows}
    rows, comparisons, reused = [], [], []
    scientific_fields = ('final_smiles', 'final_score', 'final_qed', 'solved', 'joint',
                         'replaced', 'status', 'scientific_endpoint_available')
    for supplied in selected['rows']:
        sid = supplied['sample_id']
        if sid not in by_id:
            raise ValueError(f'No cached source/candidate record for {sid}')
        entry = by_id[sid]
        source = entry['source']
        for key in ('task', 'task_index', 'replicate', 'position', 'canonical_smiles',
                    'valid', 'task_score', 'qed', 'search_seed'):
            if key in supplied and supplied[key] != source.get(key):
                raise ValueError(f'Saved input identity differs for {sid}: {key}')
        for policy, arm in zip(policies, ('task_contribution', 'permuted_contribution')):
            expected = old.get((sid, policy))
            if expected is None:
                raise ValueError(f'No cached policy outcome for {sid}/{policy}')
            attempts = list(expected['attempts'])
            used = []
            def cached_query(smiles, seed, calls, seconds):
                index = len(used)
                if index >= len(attempts):
                    raise ValueError('Replay requested an unmeasured saved candidate; no inference is allowed')
                a = attempts[index]
                if (a['candidate']['smiles'] != smiles or a['calls_limit'] != calls
                    or a['soft_seconds_requested'] != seconds):
                    raise ValueError('Saved route query differs from replay branch/budget')
                query = dict(a['query'])
                used.append(index)
                reused.append({'source': sid, 'policy': policy, 'branch': index,
                               'saved_query_path': query.get('result_path'),
                               'requested_seed': seed, 'adopted_seed': query.get('adopted_seed', seed)})
                query.update(actual_new_calls=0, cache_mode='saved_decision_replay')
                return query
            row = report_qualification(evaluate_policy(source, policy, entry['arms'][arm], cached_query), args.qualification)
            mismatch = {field: {'saved': expected.get(field), 'replay': row.get(field)}
                        for field in scientific_fields if expected.get(field) != row.get(field)}
            comparisons.append({'sample_id': sid, 'policy': policy, 'differences': mismatch})
            rows.append(row)
    result = {'scope': selected['scope'], 'mode': 'saved_decision_replay_no_new_queries',
              'selection': selected, 'qualification': args.qualification,
              'rows': rows, 'policies': summary(rows, policies),
              'comparisons': comparisons, 'reused_query_references': reused,
              'actual_new_work': {'training': 0, 'sampling': 0, 'RF_predictions': 0,
                                  'forward_reactions': 0, 'Graph_queries': 0, 'Graph_model_calls': 0},
              'route_evidence': 'Replay accepts original saved query/audit fields; it does not independently deserialize or audit graph witnesses.'}
    save(args.output / 'method_replay.json', result)
    print(json.dumps({'mode': result['mode'], 'policies': result['policies'],
                      'scientific_differences': sum(bool(x['differences']) for x in comparisons)}, indent=2))
    return int(any(x['differences'] for x in comparisons))


def source_options(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--sources', type=Path, help='JSON; default is one already measured raw source')
    group.add_argument('--smiles', help='One provided structure; invalid input remains a failed slot')
    parser.add_argument('--task', choices=('gsk3b', 'jnk3'))
    parser.add_argument('--replicate', type=int, default=0)
    parser.add_argument('--position', type=int, default=0)
    parser.add_argument('--search-seed', type=int)
    parser.add_argument('--limit', type=int, default=1, help='First supplied rows, default 1; no extra sampling')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--qualification', choices=('fixed_task', 'source_preserving'), default='fixed_task',
                        help='Fixed task RF floor (default), or historical max(source RF, task floor); QED and budgets are unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='mode', required=True)
    real = commands.add_parser('execute', help='Real CPU chemistry/RF and bounded Graph execution, with explicit optional study reuse')
    source_options(real)
    for name in ('rf-model-dir', 'forward-templates', 'retro-templates', 'stock-smiles', 'stock-keys', 'graph-model-dir'):
        real.add_argument('--' + name, type=Path, required=True)
    for name in ('rf-python', 'graph-python'):
        real.add_argument('--' + name, required=True)
    real.add_argument('--reuse-study', type=Path, help='Explicit existing study for exactly compatible saved graph queries and task-score reuse')
    real.add_argument('--device', default='cuda:0')
    cached = commands.add_parser('replay', help='Recompute B/C decisions from saved candidate scores/query records; zero inference')
    source_options(cached)
    cached.add_argument('--cached-records', type=Path, default=ROOT / 'examples/cached/candidate_records.json')
    cached.add_argument('--cached-outcomes', type=Path, default=ROOT / 'examples/cached/policy_outcomes.json')
    args = parser.parse_args()
    try:
        normalize_cli_paths(args)
        return execute(args) if args.mode == 'execute' else replay(args)
    except (ValueError, OSError, RuntimeError, KeyError) as error:
        parser.exit(2, f'{type(error).__name__}: {error}\n')


if __name__ == '__main__':
    raise SystemExit(main())
