"""Extract retained project functions; no models, third-party source, or searches.

The original full project files in sources/ are local provenance copies acquired
read-only from A6000 on 2026-10-03. This is a packaging operation, not an experiment.
"""
import ast
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'sources'
PACKAGE = ROOT / 'lfr_method'


def extract(filename, names):
    text = (SOURCE / filename).read_text(encoding='utf-8')
    lines = text.splitlines(keepends=True)
    tree = ast.parse(text)
    nodes = {n.name: n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    records, snippets = [], []
    for name in names:
        n = nodes[name]
        begin = min([n.lineno] + [d.lineno for d in n.decorator_list])
        snippets.append(''.join(lines[begin - 1:n.end_lineno]).rstrip())
        records.append({'source_file': filename, 'name': name,
                        'first_line': begin, 'last_line': n.end_lineno})
    return '\n\n\n'.join(snippets) + '\n', records


def main():
    PACKAGE.mkdir(exist_ok=True)
    (PACKAGE / '__init__.py').write_text('"""Local method execution package, version 0.2.0."""\n', encoding='utf-8')
    specifications = [
        ('common.py', 'v2_common_slot_retrieval_candidates.py',
         ['json_default', 'read', 'save', 'canon', 'fingerprints', 'features',
          'local_contributions', 'ForwardCache', 'initialize_chemistry', 'score_plan', 'retained_source'],
         'import csv, json, random, time\nfrom collections import Counter\nfrom pathlib import Path\n'
         'import numpy as np\nfrom rdkit import Chem, DataStructs\nfrom rdkit.Chem import AllChem, QED, rdChemReactions\n'
         'ARMS = ("task_contribution", "permuted_contribution")\nTHRESHOLDS = {"gsk3b": 0.10, "jnk3": 0.05}\n'),
        ('roots.py', 'v2_unified_raw8_candidates.py', ['find_roots'],
         'from collections import Counter\nimport time\nfrom rdkit import Chem\nfrom . import common\nNATIVE = {}\n'),
        ('retrieval.py', 'independent_training256_candidates.py', ['find_roots', 'generate_plan'],
         'from collections import Counter\nfrom pathlib import Path\nimport time\nfrom rdkit import Chem, DataStructs\n'
         'from rdkit.Chem import QED\nfrom . import common, roots as unified\n'
         'from .permutation import permutation\nARMS = common.ARMS\nOUT = Path(".")\n'),
        ('permutation.py', 'joint_validation512_candidates.py', ['permutation'], 'import random\n'),
        ('graph_core.py', 'joint_validation512_graph.py',
         ['read', 'save', 'identity', 'candidate_pool', 'evaluate_policy', 'ResidentGraph', 'QueryStore', 'summary'],
         'import json, math, os, pickle, time, traceback\nfrom pathlib import Path\n'
         'POLICIES = ("B_joint6plus6", "C_joint6plus6")\n'
         'THRESHOLDS = {"gsk3b": .10, "jnk3": .05}\nNORMAL = {"solved", "unsolved_within_budget"}\n'),
    ]
    origins = []
    for output, filename, names, imports in specifications:
        body, records = extract(filename, names)
        modifications = []
        if output == 'common.py':
            body = body.replace('def initialize_chemistry(repo, stock_path):',
                                'def initialize_chemistry(forward_templates, retro_templates, stock_path):')
            body = body.replace("(repo / 'data/templates.csv')", 'forward_templates')
            body = body.replace("(repo / 'data/retro_templates.csv')", 'retro_templates')
            modifications.append('Initialize chemistry takes explicit forward/retro CSV paths; function body is otherwise retained.')
        if output == 'graph_core.py':
            # Candidate target-stock exclusion remains unchanged. The shared
            # search implementation now also permits the original parent-search
            # operation, which removes the target key even if initially stocked.
            body = body.replace('if key in self.stock.stock_inchikeys:\n            raise ValueError',
                                'if key in self.stock.stock_inchikeys and request.get("target_role", "candidate") != "parent":\n            raise ValueError')
            body = body.replace("'graph_nodes': len(graph), 'target_in_stock': False", 
                                "'graph_nodes': len(graph), 'target_in_stock': key in self.stock.stock_inchikeys")
            modifications.append('Parent request permits initially stocked targets while removing their key; candidate stock guard and all search/audit settings retained.')
        doc = f'"""Core extracted from retained {filename}; see code_origins.json."""\n'
        (PACKAGE / output).write_text(doc + imports + '\n\n' + body, encoding='utf-8')
        origins.append({'module': 'lfr_method/' + output, 'functions': records, 'packaging_changes': modifications})
    stock_dir = ROOT / 'lfr_syn'
    stock_dir.mkdir(exist_ok=True)
    (stock_dir / '__init__.py').write_text('"""Project stock adapter only; no model assets."""\n', encoding='utf-8')
    stock = (SOURCE / 'shared_stock.py').read_text(encoding='utf-8')
    (stock_dir / 'shared_stock.py').write_text(stock, encoding='utf-8')
    origins.append({'module': 'lfr_syn/shared_stock.py', 'source_file': 'shared_stock.py',
                    'packaging_changes': [], 'role': 'Retained project stock adapter, complete source copy.'})
    fixed = (SOURCE / 'fixedtask_kernel.py').read_bytes()
    (PACKAGE / 'fixedtask_kernel.py').write_bytes(fixed)
    origins.append({'module': 'lfr_method/fixedtask_kernel.py', 'source_file': 'fixedtask_kernel.py',
                    'packaging_changes': [],
                    'role': 'Exact validated fixed-task wrapper; chemistry and strict execution kernels are retained.'})
    origins.append({'module': 'lfr_method/qualification.py',
                    'role': 'CLI mode dispatch and output annotations; no scientific kernel implementation.'})
    data = {'package_version': '0.2.0', 'acquired_date': '2026-10-03',
            'source_locations': json.loads((ROOT / 'source_locations.json').read_text(encoding='utf-8')),
            'modules': origins,
            'interpretation': 'Retained current project source used by completed workflows; new parameterized packaging is not asserted to be the original full execution controller or a new scientific result.'}
    (ROOT / 'code_origins.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'extracted_modules': len(origins), 'package_version': '0.2.0'}))


if __name__ == '__main__':
    main()
