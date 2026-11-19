"""Exercise relocated CLI subprocess paths without chemistry or GPU models."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


PACKAGE = Path(__file__).resolve().parent

FAKE_RF = '''import argparse, json, sys
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('mode')
p.add_argument('--inputs', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--rf-model-dir', type=Path, required=True)
a = p.parse_args()
assert a.mode == 'prepare'
assert a.inputs.is_file() and a.output.is_dir() and a.rf_model_dir.is_dir()
data = json.loads(a.inputs.read_text(encoding='utf-8'))
record = {'mode': 'fake_worker_no_scientific_calls', 'cwd': str(Path.cwd()),
          'inputs': str(a.inputs), 'output': str(a.output),
          'rf_model_dir': str(a.rf_model_dir), 'rows': len(data['rows']),
          'executable': sys.executable, 'prefix': sys.prefix, 'base_prefix': sys.base_prefix}
(a.output / 'fake_rf_observation.json').write_text(json.dumps(record), encoding='utf-8')
'''

FAKE_GRAPH = '''import argparse, json, sys
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
config = json.loads((a.output / 'method_config.json').read_text(encoding='utf-8'))
assert (a.output / 'fake_rf_observation.json').is_file()
for name in ('forward_templates', 'retro_templates', 'stock_smiles', 'stock_keys'):
    assert Path(config[name]).is_absolute() and Path(config[name]).is_file(), name
for name in ('rf_model_dir', 'graph_model_dir'):
    assert Path(config[name]).is_absolute() and Path(config[name]).is_dir(), name
result = {'mode': 'fake_worker_no_scientific_calls', 'policies': {},
          'actual_new_work': {'RF_predictions': 0, 'forward_reactions': 0,
                              'Graph_queries': 0, 'Graph_model_calls': 0},
          'cwd': str(Path.cwd()), 'output': str(a.output),
          'executable': sys.executable, 'prefix': sys.prefix, 'base_prefix': sys.base_prefix}
(a.output / 'method_result.json').write_text(json.dumps(result), encoding='utf-8')
'''


class RelocatedCliPaths(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='lfr_cli_paths_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.package = self.root / 'relocated package 空格'
        self.caller = self.root / 'caller working directory 空格'
        self.package.mkdir()
        self.caller.mkdir()
        shutil.copy2(PACKAGE / 'run_method.py', self.package / 'run_method.py')
        module_dir = self.package / 'lfr_method'
        module_dir.mkdir()
        (module_dir / '__init__.py').write_text('', encoding='utf-8')
        (module_dir / 'rf_worker.py').write_text(FAKE_RF, encoding='utf-8')
        (module_dir / 'graph_worker.py').write_text(FAKE_GRAPH, encoding='utf-8')
        shutil.copy2(PACKAGE / 'lfr_method/graph_core.py', module_dir / 'graph_core.py')
        for name in ('qualification.py', 'fixedtask_kernel.py'):
            shutil.copy2(PACKAGE / 'lfr_method' / name, module_dir / name)
        shutil.copytree(PACKAGE / 'examples', self.package / 'examples')
        assets = self.caller / 'assets with spaces'
        assets.mkdir()
        for name in ('rf', 'graph'):
            (assets / name).mkdir()
        for name in ('forward.csv', 'retro.csv', 'stock.txt', 'keys.txt'):
            (assets / name).write_text('fake asset: never consumed by a model\n', encoding='utf-8')
        self.source_path = self.caller / 'sources.json'
        self.source_path.write_text(json.dumps({'rows': [
            {'sample_id': 'software-input', 'task': 'gsk3b', 'raw_smiles': 'C'}
        ]}), encoding='utf-8')

    def invoke(self, arguments, *, extra_env=None):
        env = dict(os.environ, PYTHONUTF8='1', PYTHONNOUSERSITE='1')
        for name in ('PYTHONHOME', 'PYTHONPATH'):
            env.pop(name, None)
        if extra_env:
            env.update(extra_env)
        return subprocess.run([sys.executable, str(self.package / 'run_method.py'), *arguments],
                              cwd=self.caller, env=env, capture_output=True,
                              text=True, encoding='utf-8', timeout=30)

    def execute_arguments(self, rf_python, graph_python):
        return ['execute', '--sources', 'sources.json', '--output', 'result with spaces',
                '--rf-python', str(rf_python), '--graph-python', str(graph_python),
                '--rf-model-dir', 'assets with spaces/rf',
                '--graph-model-dir', 'assets with spaces/graph',
                '--forward-templates', 'assets with spaces/forward.csv',
                '--retro-templates', 'assets with spaces/retro.csv',
                '--stock-smiles', 'assets with spaces/stock.txt',
                '--stock-keys', 'assets with spaces/keys.txt']

    def assert_fake_execution(self, result, *, expected_python=None):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = self.caller / 'result with spaces'
        config = json.loads((output / 'method_config.json').read_text(encoding='utf-8'))
        observed = json.loads((output / 'fake_rf_observation.json').read_text(encoding='utf-8'))
        final = json.loads((output / 'method_result.json').read_text(encoding='utf-8'))
        self.assertEqual(Path(observed['cwd']).resolve(), self.package.resolve())
        self.assertEqual(Path(final['cwd']).resolve(), self.package.resolve())
        self.assertEqual(Path(observed['inputs']).resolve(), output / 'selected_inputs.json')
        self.assertEqual(Path(observed['output']).resolve(), output)
        self.assertEqual(Path(final['output']).resolve(), output)
        self.assertEqual(observed['rows'], 1)
        expected_python = Path(expected_python or sys.executable).absolute()
        for name in ('rf_python', 'graph_python'):
            self.assertEqual(Path(config[name]).resolve(), expected_python.resolve())
        self.assertEqual(Path(config['rf_model_dir']).resolve(), self.caller / 'assets with spaces/rf')
        self.assertEqual(final['mode'], 'fake_worker_no_scientific_calls')
        self.assertEqual(config['qualification'], 'fixed_task')
        self.assertEqual(config['policies'], ['B_fixedtask6plus6', 'C_fixedtask6plus6'])
        self.assertFalse((self.package / 'result with spaces').exists())
        return observed, final

    @unittest.skipUnless(os.name == 'posix', 'Requires Unix virtual-environment symlink semantics')
    def test_symlinked_venv_interpreters_keep_worker_environment(self):
        venv = self.caller / 'stdlib virtual environment with spaces'
        env = dict(os.environ, PYTHONUTF8='1', PYTHONNOUSERSITE='1')
        for name in ('PYTHONHOME', 'PYTHONPATH'):
            env.pop(name, None)
        created = subprocess.run([sys.executable, '-m', 'venv', '--without-pip', str(venv)],
                                 env=env, capture_output=True, text=True,
                                 encoding='utf-8', timeout=30)
        self.assertEqual(created.returncode, 0, created.stdout + created.stderr)
        python = venv / 'bin/python'
        self.assertTrue(python.is_symlink(), 'Regression must use an actual symlinked venv interpreter')
        self.assertNotEqual(python, python.resolve())
        self.assertIn('include-system-site-packages = false',
                      (venv / 'pyvenv.cfg').read_text(encoding='utf-8'))
        variants = {
            'absolute': str(python),
            'caller_relative': os.path.relpath(python, self.caller),
            'PATH_lookup': python.name,
        }
        search_path = str(python.parent) + os.pathsep + os.environ.get('PATH', '')
        for name, interpreter in variants.items():
            with self.subTest(path_form=name):
                observed, final = self.assert_fake_execution(
                    self.invoke(self.execute_arguments(interpreter, interpreter),
                                extra_env={'PATH': search_path}), expected_python=python)
                for worker in (observed, final):
                    self.assertEqual(Path(worker['prefix']), venv)
                    self.assertNotEqual(worker['prefix'], worker['base_prefix'])
                    self.assertEqual(Path(worker['executable']), python)

    def test_caller_relative_sources_output_and_assets_cross_real_subprocesses(self):
        self.assert_fake_execution(self.invoke(self.execute_arguments(sys.executable, sys.executable)))

    def test_caller_relative_interpreters_cross_real_subprocesses(self):
        try:
            relative = os.path.relpath(sys.executable, self.caller)
        except ValueError:
            self.skipTest('Interpreter and temporary directory use different filesystem drives')
        self.assert_fake_execution(self.invoke(self.execute_arguments(relative, relative)))

    def test_bare_interpreter_names_keep_path_lookup(self):
        name = Path(sys.executable).name
        search_path = str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', '')
        self.assert_fake_execution(self.invoke(self.execute_arguments(name, name),
                                              extra_env={'PATH': search_path}))

    def test_missing_bare_interpreter_reports_the_name(self):
        result = self.invoke(self.execute_arguments('lfr_missing_interpreter_20261004', sys.executable))
        self.assertEqual(result.returncode, 2)
        self.assertIn('Interpreter executable is not on PATH: lfr_missing_interpreter_20261004', result.stderr)
        self.assertFalse((self.caller / 'result with spaces').exists())

    def test_relocated_standard_library_replay_accepts_caller_relative_files(self):
        inputs = self.caller / 'replay input'
        shutil.copytree(self.package / 'examples', inputs)
        for source in ('source_first_raw.json', 'source_candidate_branch.json'):
            result = self.invoke(['replay', '--qualification', 'source_preserving', '--sources', f'replay input/{source}',
                                  '--cached-records', 'replay input/cached/candidate_records.json',
                                  '--cached-outcomes', 'replay input/cached/policy_outcomes.json',
                                  '--output', f'replay {source}'])
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            data = json.loads((self.caller / f'replay {source}/method_replay.json').read_text(encoding='utf-8'))
            self.assertTrue(all(not row['differences'] for row in data['comparisons']))
            self.assertTrue(all(value == 0 for value in data['actual_new_work'].values()))

    def test_explicit_dot_interpreter_path_does_not_use_path_lookup(self):
        spec = importlib.util.spec_from_file_location('lfr_cli_path_probe', self.package / 'run_method.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        previous = Path.cwd()
        try:
            os.chdir(self.caller)
            args = argparse.Namespace(output=Path('output'), rf_python='./python')
            with patch.object(module.shutil, 'which') as lookup:
                module.normalize_cli_paths(args)
                lookup.assert_not_called()
            self.assertEqual(args.rf_python, self.caller / 'python')
            self.assertEqual(args.output, self.caller / 'output')
        finally:
            os.chdir(previous)

    def test_default_fixed_mode_replays_actual_saved_new_study_branches(self):
        result = self.invoke(['replay', '--sources', str(self.package / 'examples/fixed_task/source_saved_independent.json'),
            '--cached-records', str(self.package / 'examples/fixed_task/candidate_records.json'),
            '--cached-outcomes', str(self.package / 'examples/fixed_task/policy_outcomes.json'),
            '--output', 'fixed saved replay 空格'])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        saved = json.loads((self.caller / 'fixed saved replay 空格/method_replay.json').read_text(encoding='utf-8'))
        self.assertEqual(saved['qualification'], 'fixed_task')
        self.assertTrue(all(not row['differences'] for row in saved['comparisons']))
        self.assertTrue(all(value == 0 for value in saved['actual_new_work'].values()))
        self.assertTrue(all(not row['guiding_RF_nondegradation_guarantee'] for row in saved['rows']))


if __name__ == '__main__':
    unittest.main()
