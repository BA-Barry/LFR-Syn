"""CPU source scoring and retained B/C precursor-generation functions."""
import argparse
import ast
import importlib.metadata
import inspect
import json
import math
from pathlib import Path
import pickle
import re
import time
import textwrap

from rdkit import Chem, RDLogger
from rdkit.Chem import QED
from . import common, retrieval, rf_reuse
from .cpu_costs import summarize_cpu_work


def versions():
    result = {}
    for name in ('rdkit', 'numpy', 'scikit-learn', 'PyTDC'):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    return result


def model_path(directory, task):
    p = directory / f'{task}_current.pkl'
    if not p.is_file():
        raise FileNotFoundError(f'Fixed {task} RF asset is missing: {p}')
    return p


def canonicalize(smiles):
    return common.canon(Chem.MolFromSmiles(smiles)) if isinstance(smiles, str) and smiles else None


def external_study(config, active_versions):
    if not config.get('reuse_study'):
        return None
    # The legacy records predate explicit feature fields. Compare the retained
    # function body directly, rather than assuming any future package still
    # computes the same radius2/2048 achiral vector.
    retained = Path(__file__).resolve().parents[1] / 'sources/v2_common_slot_retrieval_candidates.py'
    old = next(n for n in ast.parse(retained.read_text(encoding='utf-8')).body
               if isinstance(n, ast.FunctionDef) and n.name == 'features')
    current = ast.parse(textwrap.dedent(inspect.getsource(common.features))).body[0]
    if ast.dump(current, include_attributes=False) != ast.dump(old, include_attributes=False):
        raise ValueError('Active RF feature implementation differs from the retained formal512 implementation')
    from rdkit import __version__ as runtime_rdkit_version
    return rf_reuse.Formal512(config['reuse_study'], config['rf_model_dir'], active_versions, canonicalize,
                             runtime_rdkit_version=runtime_rdkit_version)


def prepare(args):
    data = common.read(args.inputs)
    inputs = data['rows']
    active_versions = versions()
    config_path = args.output / 'method_config.json'
    config = common.read(config_path) if config_path.exists() else {'rf_model_dir': str(args.rf_model_dir), 'reuse_study': None}
    if config['rf_model_dir'] != str(args.rf_model_dir):
        raise ValueError('Source scoring RF directory differs from the actual execution configuration')
    complete = args.output / 'sources_scored.json'
    if complete.exists():
        saved = common.read(complete)
        if (saved['supplied_inputs'] != inputs or saved['rf_model_dir'] != str(args.rf_model_dir) or
                saved['versions'] != active_versions or saved.get('reuse_study') != config.get('reuse_study')):
            raise ValueError('Saved normalized sources belong to different inputs/model directory')
        print(json.dumps({'reuse_source_scoring': str(complete), 'new_RF_predictions': 0}))
        return
    models, rows, predictions, reused = {}, [], 0, 0
    begin = time.monotonic()
    oracle = {task: {} for task in common.THRESHOLDS}
    external = external_study(config, active_versions)
    reused_score_paths = external.scores(oracle) if external else []
    for index, supplied in enumerate(inputs):
        task = supplied.get('task')
        if task not in common.THRESHOLDS:
            raise ValueError('Source task must be gsk3b or jnk3')
        task_index = ('gsk3b', 'jnk3').index(task)
        rep, position = int(supplied.get('replicate', 0)), int(supplied.get('position', index))
        if rep < 0 or position < 0 or supplied.get('task_index', task_index) != task_index:
            raise ValueError('Task/replicate/position identity is inconsistent')
        sid = supplied.get('sample_id', f'{task}-rep{rep}-r{position:03d}')
        if not isinstance(sid, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', sid):
            raise ValueError('Source sample_id must be a safe ordinary file-name identifier')
        raw = supplied.get('raw_smiles', supplied.get('SMILES', supplied.get('smiles', supplied.get('canonical_smiles'))))
        mol = Chem.MolFromSmiles(raw) if isinstance(raw, str) and raw else None
        canonical = common.canon(mol)
        valid = canonical is not None
        if 'valid' in supplied and bool(supplied['valid']) != valid:
            raise ValueError(f'Input validity differs under the RF environment: {sid}')
        if supplied.get('canonical_smiles') is not None and supplied['canonical_smiles'] != canonical:
            raise ValueError(f'Saved canonical identity differs: {sid}')
        score, qed = 0., None
        score_provenance = None
        if valid:
            path = model_path(args.rf_model_dir, task)
            cached = oracle[task].get(canonical)
            if cached is not None:
                score, score_provenance = cached['score'], cached['source']
                reused += 1
            else:
                if task not in models:
                    with path.open('rb') as stream:
                        models[task] = pickle.load(stream)
                model = models[task]
                score = float(model.predict_proba(common.features(mol).reshape(1, -1))[0, list(model.classes_).index(1)])
                predictions += 1
                score_provenance = str(complete)
                rf_reuse.add_score(oracle, task, canonical, score, complete, canonicalize)
            qed = float(QED.qed(mol))
            for field, value in [('task_score', score), ('qed', qed)]:
                if field in supplied and supplied[field] is not None and abs(float(supplied[field]) - value) > 1e-8:
                    raise ValueError(f'Saved {field} differs from actual current fixed asset: {sid}')
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError('Fixed forest returned an invalid score')
        row = dict(supplied, sample_id=sid, task=task, task_index=task_index, replicate=rep,
                   position=position, cohort_id=supplied.get('cohort_id', f'{task}-rep{rep}'),
                   raw_smiles=raw, canonical_smiles=canonical, valid=valid, task_score=score, qed=qed,
                   score_provenance=score_provenance,
                   search_seed=int(supplied.get('search_seed', 10302000 + task_index*1000000 + rep*10000 + position)),
                   graph_verified=False, graph_result=None, graph_result_path=None)
        rows.append(row)
    report = {'scope': data['scope'], 'rows': rows, 'supplied_inputs': inputs,
              'rf_model_dir': str(args.rf_model_dir), 'versions': active_versions,
              'actual_new_source_predictions': predictions, 'seconds': time.monotonic() - begin,
              'reused_source_predictions': reused, 'external_cached_score_paths': reused_score_paths,
              'reuse_study': config.get('reuse_study'),
              'training': 0, 'sampling': 0}
    common.save(complete, report)
    print(json.dumps({'source_rows': len(rows), 'new_RF_predictions': predictions}))


def candidates(args):
    config = common.read(args.output / 'method_config.json')
    rows = common.read(args.output / 'sources_with_parents.json')['rows']
    final = args.output / 'candidate_stage.json'
    if final.exists():
        saved = common.read(final)
        if [r['source'] for r in saved['sources']] != rows or saved['versions'] != versions():
            raise ValueError('Completed candidate source identities differ')
        print(json.dumps({'reuse_candidate_stage': str(final), 'new_RF_predictions': 0}))
        return
    from lfr_syn.shared_stock import InchiKeyInventory
    paths = {key: Path(config[key]) for key in ('forward_templates', 'retro_templates', 'stock_smiles', 'stock_keys', 'rf_model_dir')}
    for key in ('forward_templates', 'retro_templates', 'stock_smiles', 'stock_keys'):
        if not paths[key].is_file():
            raise FileNotFoundError(f'{key} asset is missing: {paths[key]}')
    begin = time.monotonic()
    active_versions = versions()
    external = external_study(config, active_versions)
    inventory, forward, retros = common.initialize_chemistry(paths['forward_templates'], paths['retro_templates'], paths['stock_smiles'])
    stock = InchiKeyInventory.load_from_file(paths['stock_keys'])
    observed_keys = {Chem.MolToInchiKey(entry[0]) for entry in inventory.values()}
    if observed_keys != set(stock.stock_inchikeys):
        raise ValueError('SMILES retrieval stock and InChIKey planner stock have different membership')
    retrieval.unified.NATIVE['stock_keys'] = stock.stock_inchikeys
    retrieval.OUT = args.output
    cache, roots, models = common.ForwardCache(), {}, {}
    missing_sides = lambda precursors: [i for i, s in enumerate(precursors)
                                       if Chem.MolToInchiKey(Chem.MolFromSmiles(s)) not in stock.stock_inchikeys]
    reused_chemistry = external.chemistry(config, cache, roots, forward, missing_sides) if external else None
    cache_path = args.output / 'forward_cache.json'
    if cache_path.exists():
        for r in common.read(cache_path)['records']:
            rf_reuse.add_forward(cache, r, forward, canonicalize)
    for path in sorted((args.output / 'source_records').glob('*.plan.json')):
        rf_reuse.add_roots(roots, common.read(path), path, forward, missing_sides, canonicalize)
    oracle = {task: {} for task in common.THRESHOLDS}
    ledger_path = args.output / 'new_scores.json'
    ledger = common.read(ledger_path) if ledger_path.exists() else {'models': {}, 'records': {task: [] for task in common.THRESHOLDS}, 'score_batches': [], 'versions': active_versions}
    if ledger.get('versions') != active_versions:
        raise ValueError('Resumed candidate RF environment differs')
    for task, records in ledger['records'].items():
        if records and ledger['models'].get(task) != str(paths['rf_model_dir'] / f'{task}_current.pkl'):
            raise ValueError('Resumed candidate RF model path differs')
        for record in records:
            rf_reuse.add_score(oracle, task, record['smiles'], record['score'], ledger_path, canonicalize)
    reused_score_paths = external.scores(oracle) if external else []
    scored_sources = common.read(args.output / 'sources_scored.json')
    if scored_sources['rf_model_dir'] != config['rf_model_dir'] or scored_sources['versions'] != active_versions:
        raise ValueError('Prepared source RF model/environment differs')
    for source in scored_sources['rows']:
        if source['valid']:
            rf_reuse.add_score(oracle, source['task'], source['canonical_smiles'], source['task_score'],
                               source.get('score_provenance') or args.output / 'sources_scored.json', canonicalize)
    results = []
    for source in rows:
        prefix = args.output / 'source_records' / source['sample_id']
        plan_path, complete_path = prefix.with_suffix('.plan.json'), prefix.with_suffix('.complete.json')
        failed = source['valid'] and (not source['graph_verified'] or source['task_score'] < common.THRESHOLDS[source['task']])
        if complete_path.exists():
            result = common.read(complete_path)
            if result['source'] != source:
                raise ValueError('Resumed complete candidate source differs')
        elif not failed or (source.get('graph_result') or {}).get('status') not in ('solved', 'unsolved_within_budget'):
            result = common.retained_source(source)
            result['status'] = ('retained_invalid' if not source['valid'] else 'retained_joint_success' if not failed else 'parent_execution_limitation')
            common.save(complete_path, result)
        else:
            task = source['task']
            if task not in models:
                with model_path(paths['rf_model_dir'], task).open('rb') as stream:
                    models[task] = pickle.load(stream)
            if plan_path.exists():
                plan = common.read(plan_path)
                if plan['source'] != source:
                    raise ValueError('Saved candidate plan source differs')
            else:
                plan = retrieval.generate_plan(source, models[task], inventory, forward, retros, cache, roots)
                common.save(cache_path, {'records': cache.dump()})
                common.save(plan_path, plan)
            # Logical score slots include cached predictions exactly as in the
            # original workflow. Pool membership was already written above.
            result = common.score_plan(plan, models[task], oracle, ledger, args.output, paths['rf_model_dir'])
            for arm in common.ARMS:
                for candidate in result['arms'][arm]['scored_candidates']:
                    candidate['task_pass'] = candidate['task_score'] >= max(source['task_score'], common.THRESHOLDS[task])
                result['arms'][arm]['eligible_candidates'] = sorted(
                    [c for c in result['arms'][arm]['scored_candidates'] if c['task_pass']],
                    key=lambda c: (-c['task_score'], -c['qed'], c['smiles']))
            common.save(complete_path, result)
        results.append(result)
        print(json.dumps({'source': source['sample_id'], 'status': result['status']}), flush=True)
    report = {'scope': config['scope'], 'sources': results, 'versions': active_versions,
              'elapsed_seconds': time.monotonic() - begin, 'inventory_structures': len(inventory),
              'inventory_keys': len(stock), 'external_cached_score_paths': reused_score_paths,
              'external_cached_chemistry': reused_chemistry,
              'actual_new_candidate_predictions': sum(len(rs) for rs in ledger['records'].values()),
              'actual_new_forward_attempts': sum(r['arms'][a]['costs'].get('new_forward_attempts', 0) for r in results for a in common.ARMS),
              'actual_new_forward_attempts_scope': 'Candidate formation only; excludes source replay forward calls.',
              'actual_new_cpu_work': summarize_cpu_work(results),
              'logical_score_slots': {a: sum(r['arms'][a]['costs'].get('logical_score_slots', 0) for r in results) for a in common.ARMS},
              'training': 0, 'sampling': 0, 'Graph_queries': 0}
    common.save(final, report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--rf-model-dir', type=Path, required=True)
    c = sub.add_parser('candidates')
    c.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    if args.mode == 'prepare':
        prepare(args)
    else:
        candidates(args)


if __name__ == '__main__':
    main()
