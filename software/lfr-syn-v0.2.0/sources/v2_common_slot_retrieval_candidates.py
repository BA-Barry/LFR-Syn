"""Bounded precursor retrieval on common raw slots; no sampling or Graph search.

The source list and design are prepared externally before this script runs.
Each unsolved source receives the same three retrieval arms.  Candidate order
is stock rank, root index, canonical product SMILES; only the first eight
unique quality-passing products may be scored, including cache hits.  A source
plan is saved before candidate scores are consulted for selection or predicted.
"""
import argparse
import csv
import json
import pickle
import random
import time
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, QED, rdChemReactions


DEFAULT_DEV = Path('/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development')
DEFAULT_STOCK = Path('/home/chu/trace_lfr_g0_mapping_20260721/raw_external/paroutes_v2_7341155/stock_n1.txt')
DEFAULT_MODEL = Path('/home/chu/lfr-assets/tdc-oracle-interface-20260908/oracle')
ARMS = ('morgan', 'task_contribution', 'permuted_contribution')
THRESHOLDS = {'gsk3b': 0.10, 'jnk3': 0.05}
ROOT_LIMIT = STOCK_LIMIT = 4
SCORE_LIMIT = 8


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


def source_identity(source):
    # Ordinary field comparison prevents a resumed source file being silently reused
    # for another input. No digest or separate contract is introduced.
    return {key: source.get(key) for key in ('sample_id', 'task', 'task_index',
            'position', 'canonical_smiles', 'valid', 'task_score', 'search_seed',
            'graph_result_path', 'graph_verified', 'graph_result')}


def validate_sources(rows):
    if len(rows) != 8:
        raise ValueError(f'Expected all eight raw slots, received {len(rows)}')
    if Counter(row['task'] for row in rows) != Counter({'gsk3b': 4, 'jnk3': 4}):
        raise ValueError('Expected four raw slots per task')
    if len({row['sample_id'] for row in rows}) != len(rows):
        raise ValueError('Duplicate raw source sample_id')
    for task_index, task in enumerate(THRESHOLDS):
        selected = [row for row in rows if row['task'] == task]
        if sorted(row['position'] for row in selected) != list(range(4)):
            raise ValueError(f'Expected original positions 0..3 for {task}')
        for source in selected:
            if source['task_index'] != task_index:
                raise ValueError(f'Unexpected task_index: {source["sample_id"]}')
            if source['valid']:
                canonical = canon(Chem.MolFromSmiles(source['canonical_smiles']))
                if canonical != source['canonical_smiles'] or canonical is None:
                    raise ValueError(f'Source SMILES is not valid canonical: {source["sample_id"]}')
                result = source.get('graph_result')
                if not isinstance(result, dict) or not isinstance(result.get('solved'), bool):
                    raise ValueError(f'Missing completed source Graph result: {source["sample_id"]}')
                if result.get('status') in ('error', 'execution_error', 'invalid'):
                    raise ValueError(f'Source Graph error is not an unsolved result: {source["sample_id"]}')


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


def load_forward_history(dev, out):
    cache = ForwardCache()
    saved_roots = {}
    pilot_path = dev / 'syntwins_style_root_rescue_pilot_v1/result.json'
    if pilot_path.exists():
        for row in read(pilot_path)['rows']:
            for alternative in row['alternatives']:
                cache.add(row['rxn_ref'], alternative['reactants'], alternative['novel_nonstock_products'],
                          str(pilot_path), row['source_smiles'])
    replicate_path = dev / 'syntwins_style_root_rescue_replication_v1/result.json'
    if replicate_path.exists():
        for target in read(replicate_path)['targets']:
            for row in target['rows']:
                for alternative in row['alternatives']:
                    cache.add(row['rxn_ref'], alternative['reactants'], alternative['novel_nonstock_products'],
                              str(replicate_path), target['source']['smiles'])
    previous_path = dev / 'task_attribution_retrieval_v1/generated.json'
    if previous_path.exists():
        for target in read(previous_path)['targets']:
            source_smiles = target['source']['raw_smiles']
            saved_roots[source_smiles] = {'roots': [
                {key: root[key] for key in ('rxn_ref', 'precursors', 'missing_index')}
                for root in target['roots']], 'available_root_pairs': target['available_root_pairs'],
                'provenance': str(previous_path)}
            for root in target['roots']:
                for choices in root['choices'].values():
                    for alternative in choices:
                        # generated.json stores novel, non-stock products, even
                        # where its upstream in-memory cache held complete products.
                        cache.add(root['rxn_ref'], alternative['reactants'], alternative['products'],
                                  str(previous_path), source_smiles)
    current_path = out / 'forward_cache.json'
    if current_path.exists():
        for row in read(current_path)['records']:
            cache.add(row['rxn_ref'], row['reactants'], row['products'],
                      row['provenance'], row.get('excluded_source'), row.get('forward_error'))
    return cache, saved_roots


def load_oracle_history(dev, out):
    lookup = {task: {} for task in THRESHOLDS}

    def add(task, smiles, score, provenance):
        score = float(score)
        previous = lookup[task].get(smiles)
        if previous is not None and abs(previous['score'] - score) > 1e-10:
            raise ValueError(f'Conflicting saved {task} oracle score: {smiles}')
        if previous is None:
            lookup[task][smiles] = {'score': score, 'source': provenance}

    prior_path = dev / 'task_attribution_retrieval_v1/prior_oracle_cache.json'
    for task, scores in read(prior_path)['lookup'].items():
        if task in lookup:
            for smiles, score in scores.items():
                add(task, smiles, score, str(prior_path))
    old_path = dev / 'task_attribution_retrieval_v1/new_scores.json'
    for task, records in read(old_path)['records'].items():
        for record in records:
            add(task, record['smiles'], record['score'], str(old_path))
    raw_score_path = out / 'raw_scores.json'
    if raw_score_path.exists():
        for record in read(raw_score_path)['rows']:
            if record.get('valid'):
                add(record['task'], record['canonical_smiles'], record['task_score'], str(raw_score_path))
    path = out / 'new_scores.json'
    ledger = read(path) if path.exists() else {
        'models': {}, 'records': {task: [] for task in THRESHOLDS}, 'score_batches': []}
    for task, records in ledger['records'].items():
        for record in records:
            add(task, record['smiles'], record['score'], str(path))
    return lookup, ledger


def initialize_chemistry(repo, stock_path):
    inventory = {}
    with stock_path.open(encoding='utf-8') as stream:
        for line in stream:
            molecule = Chem.MolFromSmiles(line.strip())
            smiles = canon(molecule)
            if smiles and smiles not in inventory:
                fp, environments = fingerprints(molecule)
                inventory[smiles] = (molecule, fp, environments)
    with (repo / 'data/templates.csv').open(newline='', encoding='utf-8') as stream:
        smarts = {row['rxn_ref']: row['rxn_smirks'] for row in csv.DictReader(stream)}
    forward = {ref: rdChemReactions.ReactionFromSmarts(smarts) for ref, smarts in smarts.items()}
    for reaction in forward.values():
        reaction.Initialize()
    retros = []
    with (repo / 'data/retro_templates.csv').open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            if row['rxn_ref'] in forward:
                reaction = rdChemReactions.ReactionFromSmarts(row['retro_template'])
                reaction.Initialize()
                retros.append((row['rxn_ref'], reaction))
    return inventory, forward, retros


def find_roots(source_smiles, molecule, inventory, forward, retros, saved_roots):
    if source_smiles in saved_roots:
        old = saved_roots[source_smiles]
        return old['roots'], {'available_root_pairs': old['available_root_pairs'],
                             'source_replay_forward_calls': 0, 'retro_template_calls': 0,
                             'root_seconds': 0.0, 'root_discovery_provenance': old['provenance']}
    started = time.monotonic()
    roots, replay_calls, retro_calls, retro_errors, replay_errors = {}, 0, 0, 0, 0
    for ref, retro in retros:
        reaction = forward[ref]
        if reaction.GetNumReactantTemplates() != 2:
            continue
        retro_calls += 1
        try:
            precursor_sets = retro.RunReactants((molecule,))
        except Exception:
            retro_errors += 1
            continue
        for precursor_set in precursor_sets:
            precursors = [canon(Chem.Mol(mol)) for mol in precursor_set]
            if (len(precursors) != 2 or not all(precursors)
                    or sum(smiles in inventory for smiles in precursors) != 1):
                continue
            key = (ref, tuple(precursors))
            if key in roots:
                continue
            replay_calls += 1
            try:
                recovered = any(source_smiles in {canon(Chem.Mol(product)) for product in products}
                                for products in reaction.RunReactants(tuple(
                                    Chem.MolFromSmiles(smiles) for smiles in precursors)))
            except Exception:
                recovered = False
                replay_errors += 1
            if recovered:
                roots[key] = {'rxn_ref': ref, 'precursors': precursors,
                              'missing_index': [smiles not in inventory for smiles in precursors].index(True)}
    return [roots[key] for key in sorted(roots)[:ROOT_LIMIT]], {
        'available_root_pairs': len(roots), 'source_replay_forward_calls': replay_calls,
        'retro_template_calls': retro_calls, 'retro_errors': retro_errors,
        'source_replay_errors': replay_errors, 'root_seconds': time.monotonic() - started,
        'root_discovery_provenance': 'this_experiment'}


def generate_plan(source, model, inventory, forward, retros, cache, saved_roots):
    started = time.monotonic()
    source_smiles = source['canonical_smiles']
    molecule = Chem.MolFromSmiles(source_smiles)
    source_qed = float(QED.qed(molecule))
    _, source_environments = fingerprints(molecule)
    attribution_start = time.monotonic()
    prediction, bias, contributions = local_contributions(model, features(molecule))
    if abs(prediction - source['task_score']) > 1e-8:
        raise ValueError(f'Source RF differs from saved oracle: {source["sample_id"]} {prediction}')
    bit_environments = {}
    for environment in source_environments:
        bit_environments.setdefault(environment % 2048, []).append(environment)
    weights = {environment: float(max(0.0, contributions[environment % 2048]) /
                                 len(bit_environments[environment % 2048]))
               for environment in source_environments}
    attribution_seconds = time.monotonic() - attribution_start
    roots, root_costs = find_roots(source_smiles, molecule, inventory, forward, retros, saved_roots)
    arm_costs = {arm: {'logical_forward_attempts': 0, 'new_forward_attempts': 0,
                       'cached_forward_attempts': 0, 'logical_score_slots': 0,
                       'ranking_seconds': 0.0, 'new_forward_seconds': 0.0,
                       'quality_seconds': 0.0, 'new_forward_errors': 0,
                       'logical_source_attributions': int(arm != 'morgan')}
                 for arm in ARMS}
    root_records, matching_seconds = [], 0.0
    for index, root in enumerate(roots):
        reaction = forward[root['rxn_ref']]
        missing_index = root['missing_index']
        missing = Chem.MolFromSmiles(root['precursors'][missing_index])
        missing_fp, missing_environments = fingerprints(missing)
        eligible = sorted(source_environments & missing_environments)
        real_weights = [weights[environment] for environment in eligible]
        shuffled_weights = list(real_weights)
        seed = 9273100 + 100 * source['task_index'] + 10 * source['position'] + index
        random.Random(seed).shuffle(shuffled_weights)
        permuted = dict(zip(eligible, shuffled_weights))
        matching_start = time.monotonic()
        query = reaction.GetReactantTemplate(missing_index)
        matches = [(smiles, float(DataStructs.TanimotoSimilarity(missing_fp, fp)), environments)
                   for smiles, (mol, fp, environments) in inventory.items() if mol.HasSubstructMatch(query)]
        match_seconds = time.monotonic() - matching_start
        matching_seconds += match_seconds
        choices = {}
        for arm in ARMS:
            rank_start = time.monotonic()
            if arm == 'morgan':
                ranked = sorted(((smiles, similarity, 0.0) for smiles, similarity, envs in matches),
                                key=lambda row: (-row[1], row[0]))
            else:
                use_weights = weights if arm == 'task_contribution' else permuted
                ranked = sorted(((smiles, similarity, float(sum(use_weights[env] for env in eligible if env in envs)))
                                 for smiles, similarity, envs in matches),
                                key=lambda row: (-row[2], -row[1], row[0]))
            arm_costs[arm]['ranking_seconds'] += time.monotonic() - rank_start
            choices[arm] = [{'stock': smiles, 'similarity': similarity, 'retained_weight': retained}
                            for smiles, similarity, retained in ranked[:STOCK_LIMIT]]
        root_records.append({**root, 'root_index': index, 'matching_stock_count': len(matches),
                             'matching_seconds': match_seconds,
                             'product_template_count': reaction.GetNumProductTemplates(),
                             'shared_environment_count': len(eligible),
                             'positive_weight_count': sum(weight > 0 for weight in real_weights),
                             'weight_sum': float(sum(real_weights)), 'permutation_seed': seed,
                             'weights': [{'environment': env, 'weight': weights[env],
                                          'permuted_weight': permuted[env]} for env in eligible],
                             'choices': choices})
    arms, quality_cache = {}, {}
    for arm in ARMS:
        seen, planned, attempts = set(), [], []
        for stock_rank in range(STOCK_LIMIT):
            if len(planned) >= SCORE_LIMIT:
                break
            for root in root_records:
                if len(planned) >= SCORE_LIMIT:
                    break
                choices = root['choices'][arm]
                if stock_rank >= len(choices):
                    continue
                choice = choices[stock_rank]
                reactants = root['precursors'].copy()
                reactants[root['missing_index']] = choice['stock']
                arm_costs[arm]['logical_forward_attempts'] += 1
                cached = cache.get(root['rxn_ref'], reactants, source_smiles)
                error = None
                if cached is None:
                    forward_start = time.monotonic()
                    products = set()
                    try:
                        for product_set in forward[root['rxn_ref']].RunReactants(tuple(
                                Chem.MolFromSmiles(smiles) for smiles in reactants)):
                            for product in product_set:
                                smiles = canon(Chem.Mol(product))
                                if smiles:
                                    products.add(smiles)
                    except Exception as exception:
                        error = str(exception)
                        arm_costs[arm]['new_forward_errors'] += 1
                    elapsed = time.monotonic() - forward_start
                    arm_costs[arm]['new_forward_attempts'] += 1
                    arm_costs[arm]['new_forward_seconds'] += elapsed
                    provenance = ('this_experiment_partial_products_after_forward_error' if error
                                  else 'this_experiment_complete_products')
                    cache.add(root['rxn_ref'], reactants, products, provenance, forward_error=error)
                    cached = cache.get(root['rxn_ref'], reactants, source_smiles)
                    cache_hit = False
                else:
                    arm_costs[arm]['cached_forward_attempts'] += 1
                    elapsed = 0.0
                    cache_hit = True
                    error = cached.get('forward_error')
                products = sorted(smiles for smiles in cached['products']
                                  if smiles != source_smiles and smiles not in inventory)
                evaluated = []
                for smiles in products:
                    if len(planned) >= SCORE_LIMIT:
                        break
                    if smiles in seen:
                        evaluated.append({'smiles': smiles, 'disposition': 'duplicate_within_arm'})
                        continue
                    seen.add(smiles)
                    quality_start = time.monotonic()
                    if smiles not in quality_cache:
                        quality_cache[smiles] = float(QED.qed(Chem.MolFromSmiles(smiles)))
                    qed = quality_cache[smiles]
                    arm_costs[arm]['quality_seconds'] += time.monotonic() - quality_start
                    quality_pass = qed >= source_qed - 0.05
                    candidate = {'smiles': smiles, 'qed': qed, 'quality_pass': quality_pass,
                                 'root_index': root['root_index'], 'stock_rank': stock_rank,
                                 'rxn_ref': root['rxn_ref'], 'reactants': reactants,
                                 'retrieval_similarity': choice['similarity'],
                                 'retained_weight': choice['retained_weight']}
                    if quality_pass:
                        candidate['score_slot'] = len(planned)
                        planned.append(candidate)
                    evaluated.append({**candidate, 'disposition': 'planned_score' if quality_pass else 'quality_rejected'})
                attempts.append({'root_index': root['root_index'], 'stock_rank': stock_rank,
                                 'rxn_ref': root['rxn_ref'], 'reactants': reactants,
                                 'cache_hit': cache_hit, 'forward_provenance': cached['provenance'],
                                 'cached_excluded_source': cached['excluded_source'],
                                 'new_forward_seconds': elapsed, 'forward_error': error,
                                 'products': products, 'evaluated_products': evaluated})
        arm_costs[arm]['logical_score_slots'] = len(planned)
        arms[arm] = {'planned_candidates': planned, 'forward_attempts': attempts,
                     'quota_reached': len(planned) == SCORE_LIMIT, 'costs': arm_costs[arm]}
    return {'source': source, 'status': 'candidate_plan_saved_before_oracle_lookup',
            'source_qed': source_qed, 'prediction_from_tree_paths': prediction,
            'source_bias': bias, 'source_folded_bit_collisions': sum(len(envs) > 1 for envs in bit_environments.values()),
            'roots': root_records, 'arms': arms, 'costs': {
                **root_costs, 'source_RF_path_attributions': 1, 'attribution_seconds': attribution_seconds,
                'matching_seconds': matching_seconds, 'plan_generation_seconds': time.monotonic() - started}}


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dev', type=Path, default=DEFAULT_DEV)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--sources', type=Path)
    parser.add_argument('--stock', type=Path, default=DEFAULT_STOCK)
    parser.add_argument('--models', type=Path, default=DEFAULT_MODEL)
    args = parser.parse_args()
    out = args.output or args.dev / 'common_slot_retrieval_v1'
    source_path = args.sources or out / 'sources.json'
    RDLogger.DisableLog('rdApp.*')
    rows = read(source_path)['rows']
    validate_sources(rows)
    rows = sorted(rows, key=lambda source: (source['task_index'], source['position']))
    design = read(out / 'design.json')
    final_path = out / 'candidate_stage.json'
    if final_path.exists():
        previous = read(final_path)
        if [source_identity(item['source']) for item in previous['sources']] != [source_identity(source) for source in rows]:
            raise ValueError('Completed output belongs to different raw sources; use another output directory')
        print(json.dumps({'existing_candidate_stage': str(final_path), 'new_work': False}), flush=True)
        return
    out.mkdir(parents=True, exist_ok=True)
    initialization_start = time.monotonic()
    inventory, forward, retros = initialize_chemistry(args.dev / 'external_baselines/SynTwins', args.stock)
    cache, saved_roots = load_forward_history(args.dev, out)
    oracle, score_ledger = load_oracle_history(args.dev, out)
    initialization_seconds = time.monotonic() - initialization_start
    models = {}
    started = time.monotonic()
    for task in THRESHOLDS:
        if any(source['task'] == task and source['valid'] and not source['graph_result']['solved'] for source in rows):
            with (args.models / f'{task}_current.pkl').open('rb') as stream:
                models[task] = pickle.load(stream)
    model_load_seconds = time.monotonic() - started
    invocation_path = out / 'candidate_invocations.json'
    invocations = read(invocation_path) if invocation_path.exists() else {'invocations': []}
    invocations['invocations'].append({'invocation_index': len(invocations['invocations']),
        'initialization_seconds': initialization_seconds, 'model_load_seconds': model_load_seconds,
        'inventory_size': len(inventory), 'loaded_forward_cache_records': len(cache.records),
        'loaded_oracle_scores': {task: len(scores) for task, scores in oracle.items()}})
    save(invocation_path, invocations)
    completed = []
    for source in rows:
        prefix = out / 'source_records' / f'{source["task"]}-{source["position"]:03d}'
        plan_path, completed_path = prefix.with_suffix('.plan.json'), prefix.with_suffix('.complete.json')
        if completed_path.exists():
            result = read(completed_path)
            if source_identity(result['source']) != source_identity(source):
                raise ValueError(f'Resumed source identity differs: {source["sample_id"]}')
            completed.append(result)
            print(json.dumps({'source': source['sample_id'], 'resumed_completed': True}), flush=True)
            continue
        if not source['valid'] or source['graph_result']['solved']:
            result = retained_source(source)
        else:
            if plan_path.exists():
                plan = read(plan_path)
                if source_identity(plan['source']) != source_identity(source):
                    raise ValueError(f'Resumed candidate plan identity differs: {source["sample_id"]}')
            else:
                plan = generate_plan(source, models[source['task']], inventory, forward, retros, cache, saved_roots)
                # Save reusable products first, so a persisted plan always has
                # its corresponding shared forward cache available on resume.
                save(out / 'forward_cache.json', {'records': cache.dump()})
                save(plan_path, plan)
            result = score_plan(plan, models[source['task']], oracle, score_ledger, out, args.models)
        save(completed_path, result)
        completed.append(result)
        save(out / 'candidate_progress.json', {'sources': completed, 'raw_denominator': len(rows)})
        print(json.dumps({'source': source['sample_id'], 'status': result['status'],
                          'logical_score_slots': {arm: result['arms'][arm]['costs']['logical_score_slots'] for arm in ARMS},
                          'selected': {arm: result['arms'][arm]['best_candidate'] for arm in ARMS}},
                         default=json_default), flush=True)
    requests = []
    for result in completed:
        source = result['source']
        for arm in ARMS:
            candidate = result['arms'][arm]['best_candidate']
            if candidate is not None:
                requests.append({'sample_id': source['sample_id'] + '-' + arm,
                                 'source_id': source['sample_id'], 'arm': arm, 'task': source['task'],
                                 'position': source['position'], 'SMILES': candidate['smiles'],
                                 'search_seed': source['search_seed'], 'candidate': candidate})
    summary = {arm: {
        'candidate_eligible_slots': sum(result['arms'][arm]['best_candidate'] is not None for result in completed),
        'logical_forward_attempts': sum(result['arms'][arm]['costs']['logical_forward_attempts'] for result in completed),
        'logical_score_slots': sum(result['arms'][arm]['costs']['logical_score_slots'] for result in completed),
        'ranking_seconds': sum(result['arms'][arm]['costs']['ranking_seconds'] for result in completed)} for arm in ARMS}
    new_work = {
        'source_RF_path_attributions': sum(result['costs'].get('source_RF_path_attributions', 0) for result in completed),
        'retro_template_calls': sum(result['costs'].get('retro_template_calls', 0) for result in completed),
        'source_replay_forward_calls': sum(result['costs'].get('source_replay_forward_calls', 0) for result in completed),
        'candidate_forward_attempts': sum(result['arms'][arm]['costs']['new_forward_attempts'] for result in completed for arm in ARMS),
        'candidate_forward_errors': sum(result['arms'][arm]['costs'].get('new_forward_errors', 0) for result in completed for arm in ARMS),
        'candidate_task_predictions': {task: len(records) for task, records in score_ledger['records'].items()},
        'candidate_forward_seconds': sum(result['arms'][arm]['costs']['new_forward_seconds'] for result in completed for arm in ARMS),
        'ranking_seconds': sum(result['arms'][arm]['costs']['ranking_seconds'] for result in completed for arm in ARMS),
        'candidate_score_seconds': sum(batch['score_seconds'] for batch in score_ledger['score_batches']),
        'attribution_seconds': sum(result['costs'].get('attribution_seconds', 0.0) for result in completed),
        'root_seconds': sum(result['costs'].get('root_seconds', 0.0) for result in completed),
        'matching_seconds': sum(result['costs'].get('matching_seconds', 0.0) for result in completed),
        'initialization_seconds_all_invocations': sum(row['initialization_seconds'] for row in invocations['invocations']),
        'model_load_seconds_all_invocations': sum(row['model_load_seconds'] for row in invocations['invocations']),
        'training': 0, 'Graph': 0}
    stage = {'design': design, 'source_file': str(source_path), 'raw_denominator': len(rows),
             'sources': completed, 'arms': summary, 'requests': requests, 'new_work': new_work,
             'cost_interpretation': 'Matched per-source candidate and Graph ceilings, not independent equal end-to-end wallclock. Shared initialization, source attribution, stock matching, forward and oracle caches are actual joint execution costs; logical score slots count cache hits.',
             'candidate_order': 'stock_rank outer loop, selected root inner loop, sorted canonical product SMILES; stop each arm once eight unique quality-passing products are planned',
             'graph_replacement_rule': 'Requests are unverified candidates. Replace a raw slot only after the separate same-seed Graph search solves and saved-route audit verifies it; otherwise keep the original slot.'}
    save(final_path, stage)
    print(json.dumps({'candidate_stage': str(final_path), 'arms': summary,
                      'requests': requests, 'new_work': new_work}, default=json_default), flush=True)


if __name__ == '__main__':
    main()
