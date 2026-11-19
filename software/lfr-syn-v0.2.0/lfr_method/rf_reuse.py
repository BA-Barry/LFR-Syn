"""Read only the explicitly selected formal512 study's CPU observations.

The legacy ledgers lack general feature/asset metadata. Consequently this reader
supports that retained study and its original fixed paths, not arbitrary studies
or relocated models. It never follows its lists of other historical directories.
"""
import json
import math
from pathlib import Path
import shlex


FORMAL_SCOPE = 'new_training_sources_fixed_slot_joint_failure_and_candidate_complementarity_validation'
FORMAL_RF_DIR = '/home/chu/lfr-assets/tdc-oracle-interface-20260908/oracle'
FORMAL_RF_PYTHON = '/home/chu/lfr-envs/pytdc-1.1.15-sklearn-1.2.2/bin/python'
FORMAL_DEV = '/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development'
FORMAL_STOCK = '/home/chu/trace_lfr_g0_mapping_20260721/raw_external/paroutes_v2_7341155/stock_n1.txt'
FORMAL_KEYS = '/home/chu/trace_lfr_g0_stock_identity_20260721/runs/n1_aiz_20260721/stock_inchikeys.txt'
TASKS = ('gsk3b', 'jnk3')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def rdkit_version_key(value):
    """Compare numeric release components, including metadata's zero removal."""
    if not isinstance(value, str):
        return None
    parts = value.split('.')
    return tuple(int(part) for part in parts) if all(part.isdigit() for part in parts) else None


def add_score(oracle, task, smiles, score, provenance, canonicalize):
    if task not in TASKS or not isinstance(smiles, str) or canonicalize(smiles) != smiles:
        raise ValueError(f'RF observation has incompatible task/canonical identity: {provenance}')
    score = float(score)
    if not math.isfinite(score) or not 0 <= score <= 1:
        raise ValueError(f'RF observation has an invalid score: {provenance}')
    previous = oracle[task].get(smiles)
    if previous is not None and abs(previous['score'] - score) > 1e-10:
        raise ValueError(f'Conflicting RF observations for {task}/{smiles}; no best-score selection is allowed')
    oracle[task].setdefault(smiles, {'score': score, 'source': str(provenance)})


class Formal512:
    def __init__(self, study, model_dir, versions, canonicalize, *, runtime_rdkit_version=None):
        self.study, self.canonicalize = Path(study), canonicalize
        design = read(self.study / 'design.json')
        if design.get('experiment') != 'joint_validation512_v1' or design.get('scope') != FORMAL_SCOPE:
            raise ValueError('--reuse-study CPU reuse currently supports the retained joint_validation512_v1 layout only')
        if str(model_dir) != FORMAL_RF_DIR:
            raise ValueError('Legacy CPU observations require their original unchanged fixed RF model paths')
        audit = read(self.study / 'joint_validation512_training_score_overlap_20261002.json')
        expected_rdkit = rdkit_version_key('2023.09.6')
        if (audit.get('python_executable') != FORMAL_RF_PYTHON or
                versions.get('scikit-learn') != '1.2.2' or
                rdkit_version_key(versions.get('rdkit')) != expected_rdkit or
                rdkit_version_key(audit.get('rdkit_version')) != expected_rdkit or
                (runtime_rdkit_version is not None and rdkit_version_key(runtime_rdkit_version) != expected_rdkit)):
            raise ValueError('Saved formal512 RF/features environment differs from the active RF environment')
        self.models = {task: f'{FORMAL_RF_DIR}/{task}_current.pkl' for task in TASKS}
        self.raw = read(self.study / 'raw_scores.json')
        if self.raw.get('scope') != FORMAL_SCOPE:
            raise ValueError('Raw-source score scope differs from the explicitly selected study')
        self.runs = read(self.study / 'generation_complete.json')['runs']
        expected = {(r['task'], r['replicate']) for r in design['runs']}
        observed = {(r['task'], r['replicate']) for r in self.runs}
        if observed != expected or len(observed) != len(self.runs):
            raise ValueError('Selected study has incomplete or duplicate source-training runs')

    def scores(self, oracle):
        paths = []
        path = self.study / 'new_scores.json'
        ledger = read(path)
        for task, records in ledger['records'].items():
            if task not in TASKS or (records and ledger['models'].get(task) != self.models[task]):
                raise ValueError(f'Cached candidate RF model differs: {task}')
            for record in records:
                add_score(oracle, task, record['smiles'], record['score'], path, self.canonicalize)
        paths.append(str(path))
        # Use the per-cohort model-qualified score files, and cross-check raw row
        # identities. The aggregate raw_scores file alone does not name a model.
        for run in self.runs:
            task, cohort = run['task'], run['cohort_id']
            if task not in TASKS or cohort != f"{task}-rep{run['replicate']}":
                raise ValueError('Training task/cohort identity differs')
            directory = self.study / cohort / 'oracle-only'
            path = directory / 'raw_task_scores.json'
            scores = read(path)
            selected = sorted([r for r in self.raw['rows'] if r['cohort_id'] == cohort], key=lambda r: r['position'])
            if (not scores.get('complete') or scores.get('task') != task or scores.get('cohort_id') != cohort or
                    scores.get('model_path') != self.models[task] or
                    scores.get('raw_smiles') != [r['raw_smiles'] for r in selected] or
                    scores.get('canonical_smiles') != [r['canonical_smiles'] for r in selected]):
                raise ValueError(f'Raw source RF ledger identity differs: {cohort}')
            for row in selected:
                if not row['valid']:
                    continue
                smiles = row['canonical_smiles']
                if smiles not in scores['values'] or abs(float(scores['values'][smiles])-float(row['task_score'])) > 1e-10:
                    raise ValueError(f'Raw source RF ledger value differs: {row["sample_id"]}')
                add_score(oracle, task, smiles, row['task_score'], path, self.canonicalize)
            paths.append(str(path))
            resolved = read(directory / 'resolved.json')
            components = resolved['stage'][0]['scoring']['component']
            endpoints = [e for c in components for e in c.get('ExternalProcess', {}).get('endpoint', [])]
            expected_audit = run['directory'].rstrip('/') + '/feedback.jsonl'
            for name in ('task', 'feedback'):
                found = [e for e in endpoints if e['name'] == name]
                if len(found) != 1:
                    raise ValueError('Training feedback does not use the retained task-only endpoints')
                params = found[0]['params']
                command = shlex.split(params['args'])
                if (params.get('executable') != FORMAL_RF_PYTHON or params.get('property') != name or
                        command != [f'{FORMAL_DEV.rsplit("/", 1)[0]}/code/v2_score_task_only.py', '--task', task,
                                    '--audit', expected_audit]):
                    raise ValueError(f'Training feedback scorer/task differs: {cohort}')
            path = directory / 'feedback.jsonl'
            for line in path.read_text(encoding='utf-8').splitlines():
                batch = json.loads(line)
                raw, values, outcomes = batch['smiles'], batch['raw_oracle'], batch['results']
                if len(raw) != len(values) or len(raw) != len(outcomes):
                    raise ValueError('Training feedback batch lengths differ')
                for smiles, value, result in zip(raw, values, outcomes):
                    if (result.get('reward_rule') != 'task_only' or result.get('applied') is not False or
                            abs(float(result['score'])-float(value)) > 1e-10):
                        raise ValueError('Training feedback is not an unchanged task RF observation')
                    canonical = self.canonicalize(smiles)
                    if canonical is not None:
                        add_score(oracle, task, canonical, value, path, self.canonicalize)
            paths.append(str(path))
        # Some already measured candidate values were cache hits from earlier
        # studies, so this study's new_scores ledger does not contain them. Read
        # its complete observations against its own pre-score plans; never follow
        # their external provenance paths or import only the selected winners.
        for path in sorted((self.study / 'source_records').glob('*.complete.json')):
            result = read(path)
            if result.get('status') != 'candidate_scoring_complete':
                continue
            plan_path = path.with_name(path.name.replace('.complete.json', '.plan.json'))
            plan = read(plan_path)
            task = result['source']['task']
            if task not in TASKS or result['source'] != plan['source']:
                raise ValueError('Completed candidate score source differs from its pre-score plan')
            for arm, observations in result['arms'].items():
                planned = plan['arms'][arm]['planned_candidates']
                scored = observations['scored_candidates']
                if (observations['planned_candidates'] != planned or
                        [r['smiles'] for r in scored] != [r['smiles'] for r in planned]):
                    raise ValueError('Completed candidate scores differ from the enrolled pre-score pool')
                for record in scored:
                    add_score(oracle, task, record['smiles'], record['task_score'], path, self.canonicalize)
            paths.append(str(path))
        return paths

    def chemistry(self, config, cache, roots, forward_refs, missing_sides):
        expected = {'forward_templates': f'{FORMAL_DEV}/external_baselines/SynTwins/data/templates.csv',
                    'retro_templates': f'{FORMAL_DEV}/external_baselines/SynTwins/data/retro_templates.csv',
                    'stock_smiles': FORMAL_STOCK, 'stock_keys': FORMAL_KEYS}
        if any(config.get(key) != value for key, value in expected.items()):
            raise ValueError('Legacy root/forward observations require the original template and stock paths')
        count, paths = 0, []
        path = self.study / 'forward_cache.json'
        if path.exists():
            for record in read(path)['records']:
                add_forward(cache, record, forward_refs, self.canonicalize)
                count += 1
            paths.append(str(path))
        for path in sorted((self.study / 'source_records').glob('*.plan.json')):
            add_roots(roots, read(path), path, forward_refs, missing_sides, self.canonicalize)
            paths.append(str(path))
        return {'forward_records': count, 'root_sources': len(roots), 'paths': paths,
                'scope': 'Exact canonical source roots and ordered-reactant forward observations; no saved ranking or candidate enrollment reused.'}


def add_forward(cache, record, forward_refs, canonicalize):
    if record['rxn_ref'] not in forward_refs or len(record['reactants']) != 2:
        raise ValueError('Saved forward reaction template/arity differs')
    for smiles in record['reactants'] + record['products']:
        if canonicalize(smiles) != smiles:
            raise ValueError('Saved forward reactant/product canonical identity differs')
    excluded = record.get('excluded_source')
    if excluded is not None and canonicalize(excluded) != excluded:
        raise ValueError('Saved filtered forward observation has invalid source identity')
    key = (record['rxn_ref'], tuple(record['reactants']), excluded)
    previous = cache.records.get(key)
    if previous is not None and (previous['products'] != sorted(set(record['products'])) or
                                 previous.get('forward_error') != record.get('forward_error')):
        raise ValueError('Conflicting forward observations; no outcome-based history selection is allowed')
    if previous is None:
        cache.add(record['rxn_ref'], record['reactants'], record['products'], record['provenance'],
                  excluded, record.get('forward_error'))


def add_roots(roots, plan, path, forward_refs, missing_sides, canonicalize):
    source = plan['source']['canonical_smiles']
    if canonicalize(source) != source or 'all_exact_replay_roots' not in plan.get('costs', {}):
        raise ValueError('Saved root discovery has no compatible source/replay ledger')
    selected = [{key: root[key] for key in ('rxn_ref', 'precursors', 'missing_sides')} for root in plan['roots']]
    replay = plan['costs']['all_exact_replay_roots']
    for root in replay:
        if root['rxn_ref'] not in forward_refs or len(root['precursors']) != 2:
            raise ValueError('Saved root reaction template/arity differs')
        if any(canonicalize(s) != s for s in root['precursors']) or root['missing_sides'] != missing_sides(root['precursors']):
            raise ValueError('Saved root canonical/stock membership differs')
    groups = {n: [] for n in (1, 2)}
    for n in groups:
        seen = set()
        for root in sorted(replay, key=lambda r: (r['rxn_ref'], tuple(r['precursors']))):
            if len(root['missing_sides']) == n and root['rxn_ref'] not in seen:
                groups[n].append(root)
                seen.add(root['rxn_ref'])
    expected = groups[1][:2] + groups[2][:2]
    expected += (groups[1][2:] + groups[2][2:])[:4-len(expected)]
    expected.sort(key=lambda r: (len(r['missing_sides']), r['rxn_ref'], tuple(r['precursors'])))
    if expected != selected:
        raise ValueError('Saved selected roots differ from the retained deterministic root rule')
    previous = roots.get(source)
    if previous is not None and (previous['selected'] != selected or previous['costs']['all_exact_replay_roots'] != replay):
        raise ValueError('Conflicting source root observations; no successful-root selection is allowed')
    roots.setdefault(source, {'selected': selected, 'costs': plan['costs'], 'provenance': str(path)})
