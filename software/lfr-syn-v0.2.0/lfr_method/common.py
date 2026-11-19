"""Core extracted from retained v2_common_slot_retrieval_candidates.py; see code_origins.json."""
import csv, json, random, time
from collections import Counter
from pathlib import Path
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, QED, rdChemReactions
ARMS = ("task_contribution", "permuted_contribution")
THRESHOLDS = {"gsk3b": 0.10, "jnk3": 0.05}


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f'Not JSON serializable: {type(value).__name__}')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, default=json_default, allow_nan=False), encoding='utf-8')
    temporary.replace(path)


def canon(mol):
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        if len(Chem.GetMolFrags(mol)) != 1:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def fingerprints(mol):
    return (AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
            set(AllChem.GetMorganFingerprint(mol, 2).GetNonzeroElements()))


def features(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    row = np.zeros(2048, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, row)
    return row


def local_contributions(model, row):
    """Exact class-1 RF path decomposition, using the previous implementation."""
    positive = list(model.classes_).index(1)
    contrib = np.zeros(2048, dtype=float)
    bias = prediction = 0.0
    for estimator in model.estimators_:
        tree = estimator.tree_
        values = tree.value[:, 0, :]
        probability = values[:, positive] / values.sum(axis=1)
        node = 0
        bias += probability[node]
        while tree.children_left[node] != tree.children_right[node]:
            feature = tree.feature[node]
            child = (tree.children_left[node] if row[feature] <= tree.threshold[node]
                     else tree.children_right[node])
            contrib[feature] += probability[child] - probability[node]
            node = child
        prediction += probability[node]
    count = len(model.estimators_)
    bias, prediction, contrib = bias / count, prediction / count, contrib / count
    if abs(bias + contrib.sum() - prediction) > 1e-8:
        raise ValueError('RF path decomposition does not sum to source prediction')
    return float(prediction), float(bias), contrib


class ForwardCache:
    """Retain source-specific filtered results separately from full product sets."""
    def __init__(self):
        self.records = {}

    def add(self, ref, reactants, products, provenance, excluded_source=None, forward_error=None):
        key = (ref, tuple(reactants), excluded_source)
        self.records[key] = {'rxn_ref': ref, 'reactants': list(reactants),
                             'products': sorted(set(products)), 'provenance': provenance,
                             'excluded_source': excluded_source, 'forward_error': forward_error}

    def get(self, ref, reactants, source_smiles):
        return (self.records.get((ref, tuple(reactants), None)) or
                self.records.get((ref, tuple(reactants), source_smiles)))

    def dump(self):
        return sorted(self.records.values(), key=lambda row: (
            row['rxn_ref'], row['reactants'], row['excluded_source'] or ''))


def initialize_chemistry(forward_templates, retro_templates, stock_path):
    inventory = {}
    with stock_path.open(encoding='utf-8') as stream:
        for line in stream:
            molecule = Chem.MolFromSmiles(line.strip())
            smiles = canon(molecule)
            if smiles and smiles not in inventory:
                fp, environments = fingerprints(molecule)
                inventory[smiles] = (molecule, fp, environments)
    with forward_templates.open(newline='', encoding='utf-8') as stream:
        smarts = {row['rxn_ref']: row['rxn_smirks'] for row in csv.DictReader(stream)}
    forward = {ref: rdChemReactions.ReactionFromSmarts(smarts) for ref, smarts in smarts.items()}
    for reaction in forward.values():
        reaction.Initialize()
    retros = []
    with retro_templates.open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            if row['rxn_ref'] in forward:
                reaction = rdChemReactions.ReactionFromSmarts(row['retro_template'])
                reaction.Initialize()
                retros.append((row['rxn_ref'], reaction))
    return inventory, forward, retros


def score_plan(plan, model, oracle, ledger, out, model_path):
    source = plan['source']
    task = source['task']
    pool = sorted({candidate['smiles'] for arm in ARMS for candidate in plan['arms'][arm]['planned_candidates']})
    missing = [smiles for smiles in pool if smiles not in oracle[task]]
    preexisting = len(pool) - len(missing)
    # The calling code has already persisted the plan and all eight-slot limits.
    started = time.monotonic()
    predicted = (model.predict_proba(np.stack([features(Chem.MolFromSmiles(smiles)) for smiles in missing]))[
        :, list(model.classes_).index(1)] if missing else [])
    score_seconds = time.monotonic() - started
    records = [{'smiles': smiles, 'score': float(value), 'first_source_id': source['sample_id']}
               for smiles, value in zip(missing, predicted)]
    ledger['models'][task] = str(model_path / f'{task}_current.pkl')
    ledger['records'][task].extend(records)
    batch = {'source_id': source['sample_id'], 'task': task, 'new_predictions': len(records),
             'new_prediction_smiles': missing, 'reused_unique_task_predictions': preexisting,
             'score_seconds': score_seconds, 'source_prediction_identity_checked': True}
    ledger['score_batches'].append(batch)
    save(out / 'new_scores.json', ledger)
    for record in records:
        oracle[task][record['smiles']] = {'score': record['score'], 'source': str(out / 'new_scores.json')}
    for arm in ARMS:
        scored = [{**candidate, 'task_score': oracle[task][candidate['smiles']]['score'],
                   'score_provenance': oracle[task][candidate['smiles']]['source']}
                  for candidate in plan['arms'][arm]['planned_candidates']]
        eligible = [candidate for candidate in scored
                    if candidate['task_score'] >= max(source['task_score'], THRESHOLDS[task])]
        eligible.sort(key=lambda row: (-row['task_score'], -row['qed'], row['smiles']))
        plan['arms'][arm]['scored_candidates'] = scored
        plan['arms'][arm]['best_candidate'] = eligible[0] if eligible else None
        plan['arms'][arm]['costs']['logical_score_slots'] = len(scored)
    plan['status'] = 'candidate_scoring_complete'
    plan['costs']['score_seconds_this_completion'] = score_seconds
    plan['costs']['new_task_predictions_this_completion'] = len(records)
    plan['costs']['reused_unique_task_predictions_this_completion'] = preexisting
    return plan


def retained_source(source):
    return {'source': source, 'status': 'retained_solved' if source['valid'] else 'retained_invalid',
            'arms': {arm: {'scored_candidates': [], 'best_candidate': None,
                           'costs': {'logical_forward_attempts': 0, 'new_forward_attempts': 0,
                                     'logical_score_slots': 0, 'ranking_seconds': 0.0,
                                     'new_forward_seconds': 0.0, 'logical_source_attributions': 0}}
                     for arm in ARMS}, 'costs': {}}
