"""Small ledger/cache tests; no real RF weights, route models or SSH.

The worker tests use chemistry/predictor doubles to verify which observations
are executed. Historical reader tests operate entirely on ordinary JSON files.
"""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

from lfr_method import rf_reuse as reuse


ROOT = Path(__file__).resolve().parent
VERSIONS = {'rdkit': '2023.09.6', 'scikit-learn': '1.2.2', 'numpy': '1.23.5', 'PyTDC': '1.1.15'}


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding='utf-8')


def canonical(smiles):
    if not isinstance(smiles, str) or not smiles or smiles == 'invalid':
        return None
    return 'CCO' if smiles == 'OCC' else smiles


def fixture(directory):
    study = Path(directory) / 'study'
    task, cohort = 'gsk3b', 'gsk3b-rep0'
    original_directory = f'{reuse.FORMAL_DEV}/joint_validation512_v1/{cohort}/oracle-only'
    run = {'task': task, 'replicate': 0, 'cohort_id': cohort, 'directory': original_directory}
    save(study / 'design.json', {'experiment': 'joint_validation512_v1', 'scope': reuse.FORMAL_SCOPE, 'runs': [run]})
    save(study / 'generation_complete.json', {'runs': [run]})
    save(study / 'joint_validation512_training_score_overlap_20261002.json',
         {'python_executable': reuse.FORMAL_RF_PYTHON, 'rdkit_version': VERSIONS['rdkit']})
    source = {'sample_id': 'gsk3b-rep0-r000', 'task': task, 'cohort_id': cohort, 'position': 0,
              'raw_smiles': 'CCO', 'canonical_smiles': 'CCO', 'valid': True, 'task_score': .2}
    save(study / 'raw_scores.json', {'scope': reuse.FORMAL_SCOPE, 'rows': [source]})
    models = {t: f'{reuse.FORMAL_RF_DIR}/{t}_current.pkl' for t in reuse.TASKS}
    save(study / 'new_scores.json', {'models': models,
         'records': {'gsk3b': [{'smiles': 'CCC', 'score': .3}], 'jnk3': [{'smiles': 'CCC', 'score': .6}]}})
    directory = study / cohort / 'oracle-only'
    save(directory / 'raw_task_scores.json', {'complete': True, 'task': task, 'cohort_id': cohort,
         'model_path': models[task], 'raw_smiles': ['CCO'], 'canonical_smiles': ['CCO'], 'values': {'CCO': .2}})
    endpoints = [{'name': name, 'params': {'property': name, 'executable': reuse.FORMAL_RF_PYTHON,
                 'args': f'{reuse.FORMAL_DEV.rsplit("/",1)[0]}/code/v2_score_task_only.py --task {task} --audit {original_directory}/feedback.jsonl'}}
                 for name in ('task', 'feedback')]
    save(directory / 'resolved.json', {'stage': [{'scoring': {'component': [{'ExternalProcess': {'endpoint': endpoints}}]}}]})
    save(directory / 'feedback.jsonl', {'smiles': ['CCN', 'OCC'], 'raw_oracle': [.4, .2],
         'results': [{'score': value, 'applied': False, 'reward_rule': 'task_only'} for value in (.4, .2)]})
    return study


class FakeCache:
    def __init__(self):
        self.records = {}

    def add(self, ref, reactants, products, provenance, excluded_source=None, forward_error=None):
        self.records[(ref, tuple(reactants), excluded_source)] = {
            'rxn_ref': ref, 'reactants': list(reactants), 'products': sorted(set(products)),
            'provenance': provenance, 'excluded_source': excluded_source, 'forward_error': forward_error}

    def get(self, ref, reactants, source):
        return self.records.get((ref, tuple(reactants), None)) or self.records.get((ref, tuple(reactants), source))


def chemistry_config():
    return {'forward_templates': f'{reuse.FORMAL_DEV}/external_baselines/SynTwins/data/templates.csv',
            'retro_templates': f'{reuse.FORMAL_DEV}/external_baselines/SynTwins/data/retro_templates.csv',
            'stock_smiles': reuse.FORMAL_STOCK, 'stock_keys': reuse.FORMAL_KEYS}


def root_plan():
    root = {'rxn_ref': 'R1', 'precursors': ['CCO', 'CCN'], 'missing_sides': [0]}
    return {'source': {'canonical_smiles': 'CCCC'}, 'roots': [{**root, 'sides': ['old ranking deliberately ignored']}],
            'costs': {'all_exact_replay_roots': [root], 'source_replay_forward_calls': 2,
                      'retro_template_calls': 9, 'selected_root_counts': {'1': 1}}}


class HistoryTests(unittest.TestCase):
    def test_all_score_ledgers_reused_without_cross_task_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            oracle = {t: {} for t in reuse.TASKS}
            paths = reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical).scores(oracle)
            self.assertEqual({s: r['score'] for s, r in oracle['gsk3b'].items()}, {'CCC': .3, 'CCO': .2, 'CCN': .4})
            self.assertEqual(oracle['jnk3']['CCC']['score'], .6)
            self.assertEqual(len(paths), 3)
            self.assertTrue(oracle['gsk3b']['CCO']['source'].endswith('raw_task_scores.json'))

    def test_changed_models_or_feature_environment_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            with self.assertRaisesRegex(ValueError, 'fixed RF model paths'):
                reuse.Formal512(study, '/a/different/model', VERSIONS, canonical)
            for key, value in [('rdkit', '2025.9.2'), ('scikit-learn', '1.7.0')]:
                with self.assertRaisesRegex(ValueError, 'environment differs'):
                    reuse.Formal512(study, reuse.FORMAL_RF_DIR, dict(VERSIONS, **{key: value}), canonical)

    def test_equivalent_distribution_runtime_rdkit_version_spelling_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            reader = reuse.Formal512(study, reuse.FORMAL_RF_DIR, dict(VERSIONS, rdkit='2023.9.6'),
                                     canonical, runtime_rdkit_version='2023.09.6')
            oracle = {task: {} for task in reuse.TASKS}
            reader.scores(oracle)
            self.assertEqual(oracle['gsk3b']['CCO']['score'], .2)

    def test_different_loaded_rdkit_runtime_is_rejected_despite_matching_distribution(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            with self.assertRaisesRegex(ValueError, 'environment differs'):
                reuse.Formal512(study, reuse.FORMAL_RF_DIR, dict(VERSIONS, rdkit='2023.9.6'),
                                canonical, runtime_rdkit_version='2025.09.2')

    def test_wrong_raw_model_is_not_accepted_via_aggregate_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            path = study / 'gsk3b-rep0/oracle-only/raw_task_scores.json'
            data = reuse.read(path)
            data['model_path'] = f'{reuse.FORMAL_RF_DIR}/jnk3_current.pkl'
            save(path, data)
            with self.assertRaisesRegex(ValueError, 'Raw source RF ledger identity'):
                reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical).scores({t: {} for t in reuse.TASKS})

    def test_conflicts_rejected_instead_of_choosing_higher_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            path = study / 'gsk3b-rep0/oracle-only/feedback.jsonl'
            data = reuse.read(path)
            data['raw_oracle'][1] = data['results'][1]['score'] = .9
            save(path, data)
            with self.assertRaisesRegex(ValueError, 'Conflicting RF observations'):
                reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical).scores({t: {} for t in reuse.TASKS})

    def test_training_feedback_requires_actual_task_only_endpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            path = study / 'gsk3b-rep0/oracle-only/resolved.json'
            data = reuse.read(path)
            data['stage'][0]['scoring']['component'][0]['ExternalProcess']['endpoint'][0]['params']['args'] = 'different.py --task jnk3'
            save(path, data)
            with self.assertRaisesRegex(ValueError, 'scorer/task differs'):
                reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical).scores({t: {} for t in reuse.TASKS})

    def test_no_other_history_directory_is_followed(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            data = reuse.read(study / 'raw_scores.json')
            data['score_history_files'] = ['/must/not/be/opened/new_scores.json']
            save(study / 'raw_scores.json', data)
            oracle = {t: {} for t in reuse.TASKS}
            reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical).scores(oracle)
            self.assertEqual(len(oracle['gsk3b']), 3)

    def test_complete_candidate_ledger_reuses_all_observations_not_only_winner(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            source = {'task': 'gsk3b', 'canonical_smiles': 'CCO', 'sample_id': 's'}
            planned = [{'smiles': 'CCCC'}, {'smiles': 'CCCO'}]
            plan = {'source': source, 'arms': {'task_contribution': {'planned_candidates': planned}}}
            complete = {'source': source, 'status': 'candidate_scoring_complete', 'arms': {'task_contribution': {
                'planned_candidates': planned, 'scored_candidates': [
                    {'smiles': 'CCCC', 'task_score': .01, 'score_provenance': '/unavailable/earlier.json'},
                    {'smiles': 'CCCO', 'task_score': .95, 'score_provenance': '/unavailable/earlier.json'}],
                'best_candidate': {'smiles': 'CCCO', 'task_score': .95}}}}
            save(study / 'source_records/s.plan.json', plan)
            save(study / 'source_records/s.complete.json', complete)
            oracle = {t: {} for t in reuse.TASKS}
            paths = reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical).scores(oracle)
            self.assertEqual(oracle['gsk3b']['CCCC']['score'], .01)
            self.assertEqual(oracle['gsk3b']['CCCO']['score'], .95)
            self.assertEqual(len(paths), 4)
            self.assertTrue(oracle['gsk3b']['CCCC']['source'].endswith('s.complete.json'))
            complete['arms']['task_contribution']['scored_candidates'].reverse()
            save(study / 'source_records/s.complete.json', complete)
            with self.assertRaisesRegex(ValueError, 'pre-score pool'):
                reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical).scores({t: {} for t in reuse.TASKS})

    def test_filtered_forward_result_is_not_promoted_to_full_product_set(self):
        cache = FakeCache()
        record = {'rxn_ref': 'R1', 'reactants': ['CCO', 'CCN'], 'products': ['CCC'],
                  'provenance': 'old', 'excluded_source': 'CCCC'}
        reuse.add_forward(cache, record, {'R1'}, canonical)
        self.assertIsNotNone(cache.get('R1', ['CCO', 'CCN'], 'CCCC'))
        self.assertIsNone(cache.get('R1', ['CCO', 'CCN'], 'CCCCC'))
        self.assertIsNone(cache.get('R1', ['CCN', 'CCO'], 'CCCC'))

    def test_forward_template_canonical_identity_and_conflicts(self):
        record = {'rxn_ref': 'R1', 'reactants': ['CCO', 'CCN'], 'products': ['CCC'], 'provenance': 'old'}
        cache = FakeCache()
        with self.assertRaisesRegex(ValueError, 'template/arity'):
            reuse.add_forward(cache, record, {'R2'}, canonical)
        with self.assertRaisesRegex(ValueError, 'canonical identity'):
            reuse.add_forward(cache, dict(record, reactants=['OCC', 'CCN']), {'R1'}, canonical)
        reuse.add_forward(cache, record, {'R1'}, canonical)
        with self.assertRaisesRegex(ValueError, 'Conflicting forward'):
            reuse.add_forward(cache, dict(record, products=['CCCO']), {'R1'}, canonical)

    def test_cached_roots_do_not_carry_saved_ranking_and_match_current_stock(self):
        roots = {}
        reuse.add_roots(roots, root_plan(), 'old.plan.json', {'R1'}, lambda _: [0], canonical)
        self.assertNotIn('sides', roots['CCCC']['selected'][0])
        with self.assertRaisesRegex(ValueError, 'stock membership differs'):
            reuse.add_roots({}, root_plan(), 'old', {'R1'}, lambda _: [0, 1], canonical)
        changed = root_plan()
        changed['roots'] = []
        with self.assertRaisesRegex(ValueError, 'deterministic root rule'):
            reuse.add_roots({}, changed, 'old', {'R1'}, lambda _: [0], canonical)

    def test_actual_selected_study_chemistry_import_and_asset_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = fixture(tmp)
            save(study / 'source_records/s.plan.json', root_plan())
            save(study / 'forward_cache.json', {'records': [{'rxn_ref': 'R1', 'reactants': ['CCO', 'CCN'],
                 'products': ['CCC'], 'provenance': 'old', 'excluded_source': 'CCCC'}]})
            history = reuse.Formal512(study, reuse.FORMAL_RF_DIR, VERSIONS, canonical)
            cache, roots = FakeCache(), {}
            report = history.chemistry(chemistry_config(), cache, roots, {'R1'}, lambda _: [0])
            self.assertEqual((report['forward_records'], report['root_sources']), (1, 1))
            self.assertEqual(roots['CCCC']['selected'][0]['precursors'], ['CCO', 'CCN'])
            with self.assertRaisesRegex(ValueError, 'template and stock paths'):
                history.chemistry(dict(chemistry_config(), stock_smiles='other_stock'), cache, roots, {'R1'}, lambda _: [0])

    def test_empty_replayed_root_set_is_reused_without_rediscovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            plan = {'source': {'canonical_smiles': 'CCCC'}, 'roots': [],
                    'costs': {'all_exact_replay_roots': [], 'selected_root_counts': {}}}
            roots = {}
            reuse.add_roots(roots, plan, 'empty.plan.json', {'R1'}, lambda _: [0], canonical)
            module, unified = load_retrieval_double()
            selected, costs = module.find_roots({'canonical_smiles': 'CCCC'}, {}, {}, [], roots)
            self.assertEqual(selected, [])
            self.assertEqual(costs['retro_template_calls'], 0)
            self.assertEqual(costs['source_replay_forward_calls'], 0)
            self.assertEqual(costs['root_discovery_provenance'], 'empty.plan.json')
            unified.find_roots.assert_not_called()


def load_worker_double():
    # Import the real worker with minimal chemistry doubles, rather than making
    # RDKit or real model assets a prerequisite for these ledger semantics tests.
    package = types.ModuleType('_rf_reuse_test_pkg')
    package.__path__ = [str(ROOT / 'lfr_method')]
    common = types.ModuleType('_rf_reuse_test_pkg.common')
    common.THRESHOLDS = {'gsk3b': .1, 'jnk3': .05}
    common.read, common.save = reuse.read, save
    common.canon = canonical
    common.features = lambda mol: types.SimpleNamespace(reshape=lambda *_: [[1]])
    retrieval = types.ModuleType('_rf_reuse_test_pkg.retrieval')
    rdkit = types.ModuleType('rdkit')
    chem = types.ModuleType('rdkit.Chem')
    chem.MolFromSmiles = canonical
    chem.QED = types.SimpleNamespace(qed=lambda _: .5)
    rdkit.Chem, rdkit.RDLogger = chem, types.SimpleNamespace(DisableLog=lambda _: None)
    modules = {'_rf_reuse_test_pkg': package, '_rf_reuse_test_pkg.common': common,
               '_rf_reuse_test_pkg.retrieval': retrieval, '_rf_reuse_test_pkg.rf_reuse': reuse,
               'rdkit': rdkit, 'rdkit.Chem': chem}
    spec = importlib.util.spec_from_file_location('_rf_reuse_test_pkg.rf_worker', ROOT / 'lfr_method/rf_worker.py')
    module = importlib.util.module_from_spec(spec)
    with patch.dict(sys.modules, modules):
        spec.loader.exec_module(module)
    return module


def load_retrieval_double():
    from unittest.mock import Mock
    package = types.ModuleType('_retrieval_reuse_test_pkg')
    package.__path__ = [str(ROOT / 'lfr_method')]
    common = types.ModuleType('_retrieval_reuse_test_pkg.common')
    common.ARMS = ('task_contribution', 'permuted_contribution')
    roots = types.ModuleType('_retrieval_reuse_test_pkg.roots')
    roots.find_roots = Mock(side_effect=AssertionError('Old root discovery must not be executed'))
    rdkit = types.ModuleType('rdkit')
    chem = types.ModuleType('rdkit.Chem')
    chem.QED = types.SimpleNamespace()
    rdkit.Chem, rdkit.DataStructs = chem, types.SimpleNamespace()
    modules = {'_retrieval_reuse_test_pkg': package, '_retrieval_reuse_test_pkg.common': common,
               '_retrieval_reuse_test_pkg.roots': roots, 'rdkit': rdkit, 'rdkit.Chem': chem}
    spec = importlib.util.spec_from_file_location('_retrieval_reuse_test_pkg.retrieval', ROOT / 'lfr_method/retrieval.py')
    module = importlib.util.module_from_spec(spec)
    with patch.dict(sys.modules, modules):
        spec.loader.exec_module(module)
    return module, roots


def load_scoring_double():
    # Execute the real score_plan implementation with one-bit test features.
    # Its predictor records exactly the missing structures supplied to it.
    import numpy as np
    rdkit = types.ModuleType('rdkit')
    chem = types.ModuleType('rdkit.Chem')
    chem.MolFromSmiles = lambda s: s
    for name in ('AllChem', 'QED', 'rdChemReactions'):
        setattr(chem, name, types.SimpleNamespace())
    rdkit.Chem, rdkit.DataStructs = chem, types.SimpleNamespace()
    spec = importlib.util.spec_from_file_location('_common_reuse_test', ROOT / 'lfr_method/common.py')
    module = importlib.util.module_from_spec(spec)
    with patch.dict(sys.modules, {'rdkit': rdkit, 'rdkit.Chem': chem}):
        spec.loader.exec_module(module)
    module.features = lambda s: np.array([len(s)], dtype=np.float32)
    return module


class ExecutionTests(unittest.TestCase):
    def test_completed_candidate_resume_preserves_original_cost_records_without_model_work(self):
        worker = load_worker_double()
        source = {'sample_id': 'source', 'task': 'gsk3b', 'valid': True}
        record = {'source': source, 'status': 'candidate_scoring_complete',
                  'costs': {'source_replay_forward_calls': 121, 'retro_template_calls': 1071,
                            'source_RF_path_attributions': 1},
                  'arms': {'task_contribution': {'costs': {'new_forward_attempts': 6}},
                           'permuted_contribution': {'costs': {'new_forward_attempts': 7}}}}
        for has_summary in (False, True):
            with self.subTest(has_summary=has_summary), tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp)
                save(out / 'method_config.json', {'scope': 'saved_cpu_cost_test'})
                save(out / 'sources_with_parents.json', {'rows': [source]})
                stage = {'sources': [record], 'versions': VERSIONS,
                         'actual_new_forward_attempts': 13}
                if has_summary:
                    stage['actual_new_cpu_work'] = worker.summarize_cpu_work(stage['sources'])
                path = out / 'candidate_stage.json'
                save(path, stage)
                original = path.read_bytes()
                with patch.object(worker, 'versions', return_value=VERSIONS), \
                        patch.object(worker.pickle, 'load') as load:
                    worker.candidates(types.SimpleNamespace(output=out))
                    load.assert_not_called()
                self.assertEqual(path.read_bytes(), original)
                summary = worker.summarize_cpu_work(reuse.read(path)['sources'])
                self.assertEqual(summary['forward_reactions_total']['count'], 134)

    def test_candidate_parent_unknown_is_retained_without_rf_or_chemistry_queries(self):
        worker = load_worker_double()
        worker.common.ARMS = ('task_contribution', 'permuted_contribution')
        worker.common.ForwardCache = FakeCache
        worker.common.initialize_chemistry = lambda *_: ({}, {}, [])
        worker.common.retained_source = lambda source: {'source': source, 'costs': {},
            'arms': {a: {'scored_candidates': [], 'costs': {}} for a in worker.common.ARMS}}
        worker.retrieval.unified = types.SimpleNamespace(NATIVE={})
        shared_stock = types.ModuleType('lfr_syn.shared_stock')
        shared_stock.InchiKeyInventory = types.SimpleNamespace(load_from_file=lambda _: type(
            'EmptyStock', (), {'stock_inchikeys': set(), '__len__': lambda _: 0})())
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            config = {'scope': 'test', 'reuse_study': None}
            for name in ('forward_templates', 'retro_templates', 'stock_smiles', 'stock_keys', 'rf_model_dir'):
                asset = out / name
                asset.touch()
                config[name] = str(asset)
            source = {'sample_id': 'unknown', 'task': 'gsk3b', 'valid': True, 'canonical_smiles': 'CCO',
                      'task_score': .05, 'graph_verified': False, 'graph_result': None}
            save(out / 'method_config.json', config)
            save(out / 'sources_with_parents.json', {'rows': [source]})
            save(out / 'sources_scored.json', {'rows': [source], 'rf_model_dir': config['rf_model_dir'], 'versions': VERSIONS})
            with patch.dict(sys.modules, {'lfr_syn.shared_stock': shared_stock}), \
                    patch.object(worker, 'versions', return_value=VERSIONS), \
                    patch.object(worker.pickle, 'load') as loads:
                worker.candidates(types.SimpleNamespace(output=out))
                loads.assert_not_called()
            report = reuse.read(out / 'candidate_stage.json')
            self.assertEqual(report['sources'][0]['status'], 'parent_execution_limitation')
            self.assertEqual(report['actual_new_candidate_predictions'], 0)
            self.assertEqual(report['actual_new_forward_attempts'], 0)

    def test_candidate_predictions_are_missing_only_but_all_logical_slots_remain(self):
        import numpy as np
        common = load_scoring_double()
        observed = []
        class Predictor:
            classes_ = [0, 1]
            def predict_proba(self, rows):
                observed.extend(rows[:, 0].tolist())
                return np.array([[.2, .8]] * len(rows))
        oracle = {t: {} for t in reuse.TASKS}
        reuse.add_score(oracle, 'gsk3b', 'CCC', .3, 'old/new_scores.json', canonical)
        reuse.add_score(oracle, 'gsk3b', 'CCN', .4, 'old/feedback.jsonl', canonical)
        ledger = {'models': {}, 'records': {t: [] for t in reuse.TASKS}, 'score_batches': []}
        plan = {'source': {'task': 'gsk3b', 'sample_id': 'source', 'task_score': .2},
                'arms': {a: {'planned_candidates': [{'smiles': s, 'qed': .5} for s in ('CCC', 'CCN', 'CCCO')],
                             'costs': {}} for a in common.ARMS}, 'costs': {}}
        with tempfile.TemporaryDirectory() as tmp:
            result = common.score_plan(plan, Predictor(), oracle, ledger, Path(tmp), Path(reuse.FORMAL_RF_DIR))
            self.assertEqual(observed, [4.])
            self.assertEqual([r['smiles'] for r in ledger['records']['gsk3b']], ['CCCO'])
            self.assertEqual(result['costs']['new_task_predictions_this_completion'], 1)
            self.assertEqual(result['costs']['reused_unique_task_predictions_this_completion'], 2)
            self.assertEqual([result['arms'][a]['costs']['logical_score_slots'] for a in common.ARMS], [3, 3])
            common.score_plan(plan, Predictor(), oracle, ledger, Path(tmp), Path(reuse.FORMAL_RF_DIR))
            self.assertEqual(observed, [4.])

    def test_prepare_executes_only_missing_predictions_and_keeps_invalid_duplicates(self):
        worker = load_worker_double()
        predicted = []
        class Predictor:
            classes_ = [0, 1]
            def predict_proba(self, rows):
                predicted.append(rows)
                return {(0, 1): .7}
        class History:
            def scores(self, oracle):
                reuse.add_score(oracle, 'gsk3b', 'CCO', .2, 'old/raw_task_scores.json', canonical)
                return ['old/raw_task_scores.json']
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            model_dir, output = directory / 'models', directory / 'output'
            model_dir.mkdir()
            (model_dir / 'gsk3b_current.pkl').touch()
            inputs = directory / 'inputs.json'
            rows = [{'task': 'gsk3b', 'raw_smiles': s} for s in ('OCC', 'CCC', 'CCC', 'invalid')]
            save(inputs, {'scope': 'test', 'rows': rows})
            save(output / 'method_config.json', {'rf_model_dir': str(model_dir), 'reuse_study': 'explicit_study'})
            args = types.SimpleNamespace(inputs=inputs, output=output, rf_model_dir=model_dir)
            with patch.object(worker, 'versions', return_value=VERSIONS), \
                    patch.object(worker, 'external_study', return_value=History()) as selected, \
                    patch.object(worker.pickle, 'load', return_value=Predictor()) as loads:
                worker.prepare(args)
            result = reuse.read(output / 'sources_scored.json')
            self.assertEqual(len(predicted), 1)
            self.assertEqual(loads.call_count, 1)
            self.assertEqual(selected.call_args.args[0]['reuse_study'], 'explicit_study')
            self.assertEqual(result['actual_new_source_predictions'], 1)
            self.assertEqual(result['reused_source_predictions'], 2)
            self.assertEqual([r['task_score'] for r in result['rows']], [.2, .7, .7, 0.])
            self.assertEqual([r['valid'] for r in result['rows']], [True, True, True, False])
            self.assertEqual(len(result['rows']), 4)
            with patch.object(worker, 'versions', return_value=VERSIONS), patch.object(worker.pickle, 'load') as loads:
                worker.prepare(args)
                loads.assert_not_called()


if __name__ == '__main__':
    unittest.main()
