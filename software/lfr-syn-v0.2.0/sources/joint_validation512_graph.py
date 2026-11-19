"""Resident Graph2Edits evaluation of the prespecified joint-validation policies.

Run with the existing Graph2Edits environment and PYTHONPATH=/home/chu/lfr-syn/src.
Imports are standard-library-only until the production backend is constructed;
control-flow tests use an injected executor and never load a model.

Use one writer per output directory. Independent --cohort workers must have
separate output directories. Search witnesses, errors and interrupted queries
are retained. A saved incomplete query is never automatically searched again.
"""
import argparse
import importlib.metadata
import json
import math
import os
import pickle
import re
import signal
import time
import traceback
from pathlib import Path

BASE = Path('/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914')
DEV = BASE / 'development'
STUDY = DEV / 'joint_validation512_v1'
STOCK = '/home/chu/trace_lfr_g0_stock_identity_20260721/runs/n1_aiz_20260721/stock_inchikeys.txt'
MODEL = '/home/chu/trace_lfr_syntheseus_smoke_20260721/model-cache/Graph2Edits_backward'
POLICIES = (
    'B_route12', 'B_joint12', 'B_joint6plus6', 'Mix_joint6plus6',
    'C_joint6plus6', 'A_joint6plus6', 'SynTwins_joint6plus6',
)
ARMS = {'B_route12': 'task_contribution', 'B_joint12': 'task_contribution',
        'B_joint6plus6': 'task_contribution', 'Mix_joint6plus6': 'mix',
        'C_joint6plus6': 'permuted_contribution', 'A_joint6plus6': 'morgan',
        'SynTwins_joint6plus6': 'syntwins_depth3_k10'}
THRESHOLDS = {'gsk3b': .10, 'jnk3': .05}
COMPATIBLE = {'stock': STOCK, 'model_dir': MODEL, 'num_top_results': 20,
              'max_expansion_depth': 12,
              'target_exclusion': 'remove_target_standard_inchikey_from_inventory_per_search',
              'cache': 'reset_between_every_target_including_duplicates'}
NORMAL = {'solved', 'unsolved_within_budget'}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def identity(source):
    return {k: source.get(k) for k in ('sample_id', 'task', 'task_index', 'replicate',
        'cohort_id', 'position', 'canonical_smiles', 'valid', 'task_score', 'qed',
        'search_seed', 'graph_verified', 'graph_result_path')}


def ordered_policies(cohort_index):
    shift = cohort_index % len(POLICIES)
    return list(POLICIES[shift:] + POLICIES[:shift])


def candidate_pool(source, arm):
    scored = arm.get('scored_candidates', [])
    if len(scored) > 8:
        raise ValueError('More than eight scored candidates in a policy source')
    seen, pool = set(), []
    floor = max(float(source['task_score']), THRESHOLDS[source['task']])
    for c in scored:
        if c['smiles'] in seen:
            raise ValueError('Duplicate scored candidate identity')
        seen.add(c['smiles'])
        score, qed = float(c['task_score']), float(c['qed'])
        if not math.isfinite(score) or not math.isfinite(qed):
            raise ValueError('Non-finite candidate score/QED')
        eligible = (bool(c.get('task_pass')) and score >= floor
                    and qed >= float(source['qed']) - .05
                    and c['smiles'] != source['canonical_smiles'])
        if c.get('task_pass') and not eligible:
            raise ValueError('Candidate task_pass conflicts with prespecified rules')
        if eligible:
            pool.append(c)
    return sorted(pool, key=lambda c: (-c['task_score'], -c['qed'], c['smiles']))[:2]


def evaluate_policy(source, policy, arm, execute, saved=None, persist=None):
    """Evaluate one slot with saved query costs, not walltime of cache lookup.

    execute(smiles, seed, calls, soft_seconds) returns a query record containing
    status, verified_solved, model_calls, search_seconds and actual_new_calls.
    saved attempts permit continuation of a *policy decision* after interruption;
    they never permit repeating a chemical search.
    """
    if policy not in POLICIES:
        raise ValueError('Unknown policy')
    source_solved = bool(source.get('graph_verified'))
    source_joint = source_solved and source['task_score'] >= THRESHOLDS[source['task']]
    row = {'sample_id': source['sample_id'], 'cohort_id': source['cohort_id'],
        'task': source['task'], 'replicate': source['replicate'], 'position': source['position'],
        'source_identity': identity(source), 'policy': policy, 'source_solved': source_solved,
        'source_joint': source_joint, 'source_smiles': source.get('canonical_smiles'),
        'source_score': source['task_score'], 'source_qed': source.get('qed'),
        'final_smiles': source.get('canonical_smiles'), 'final_score': source['task_score'],
        'final_qed': source.get('qed'), 'solved': source_solved, 'joint': source_joint,
        'fallback_solved': source_solved, 'fallback_joint': source_joint,
        'replaced': False, 'attempts': [], 'status': 'processing', 'planning_soft_seconds': 45.,
        'logical_model_calls': 0, 'actual_new_calls': 0, 'logical_search_seconds': 0.,
        'logical_cost_observation_complete': True,
        'soft_budget_overrun_seconds': 0., 'scientific_endpoint_available': True}
    def unknown(status):
        row.update(status=status, scientific_endpoint_available=False, solved=None, joint=None,
                   endpoint_lower_solved=row['fallback_solved'], endpoint_lower_joint=row['fallback_joint'],
                   endpoint_upper_solved=True, endpoint_upper_joint=True)
    if saved:
        if saved['source_identity'] != identity(source) or saved['policy'] != policy:
            raise ValueError('Saved policy source identity differs')
        if saved.get('status') != 'processing':
            return saved
    if not source['valid'] or source_joint or (policy == 'B_route12' and source_solved):
        row['status'] = 'source_retained'
        if persist:
            persist(row)
        return row
    if source.get('graph_result', {}).get('status') not in NORMAL:
        unknown('parent_execution_limitation')
        if persist:
            persist(row)
        return row
    if arm is None:
        unknown('pending_candidate_stage')
        if persist:
            persist(row)
        return row
    native_status = arm.get('native_status')
    if native_status in ('execution_limitation', 'execution_error', 'interrupted'):
        unknown('candidate_execution_limitation')
        if persist:
            persist(row)
        return row
    pool = candidate_pool(source, arm)
    row['qualified_candidate_count_up_to_two'] = len(pool)
    if not pool:
        row['status'] = 'no_qualified_candidate'
        if persist:
            persist(row)
        return row
    is_split = policy.endswith('6plus6')
    row['prespecified_candidates'] = pool[:2 if is_split else 1]
    prior = saved.get('attempts', []) if saved else []
    for branch, candidate in enumerate(row['prespecified_candidates']):
        remaining = max(0., 45. - row['logical_search_seconds'])
        if remaining <= 0.:
            row['second_skipped'] = 'no_remaining_planning_time'
            break
        cap = 6 if is_split else 12
        if branch < len(prior):
            previous = prior[branch]
            if (previous['candidate']['smiles'] != candidate['smiles']
                or previous['calls_limit'] != cap
                or previous['soft_seconds_requested'] != remaining):
                raise ValueError('Saved candidate/branch/time budget differs')
            query = previous['query']
        else:
            query = execute(candidate['smiles'], int(source['search_seed']), cap, remaining)
        used = float(query.get('search_seconds', 0.))
        calls = int(query.get('model_calls', 0))
        if used < 0 or not math.isfinite(used) or calls < 0:
            raise ValueError('Invalid search time or call expenditure')
        attempt = {'branch': branch, 'candidate': candidate, 'calls_limit': cap,
                   'soft_seconds_requested': remaining, 'query': query}
        row['attempts'].append(attempt)
        row['logical_model_calls'] += calls
        row['actual_new_calls'] += int(query.get('actual_new_calls', 0))
        row['logical_search_seconds'] += used
        row['logical_cost_observation_complete'] &= query.get('cost_observation_complete', True)
        row['soft_budget_overrun_seconds'] = max(0., row['logical_search_seconds'] - 45.)
        if persist:
            persist(row)
        if calls > cap:
            row['model_calls_overrun'] = calls - cap
            unknown('reaction_call_budget_limitation')
            break
        if query.get('status') not in NORMAL:
            unknown('search_execution_limitation')
            break
        if query.get('verified_solved'):
            row.update(final_smiles=candidate['smiles'], final_score=candidate['task_score'],
                       final_qed=candidate['qed'], solved=True, joint=True,
                       replaced=True, chosen_branch=branch, status='replacement_solved')
            break
        if query.get('status') == 'solved':
            unknown('route_audit_limitation')
            break
    else:
        row['status'] = 'candidate_budget_unresolved'
    if row['status'] == 'processing':
        row['status'] = 'candidate_budget_unresolved'
    if row['logical_model_calls'] > 12:
        row['combined_model_calls_overrun'] = row['logical_model_calls'] - 12
        unknown('reaction_call_budget_limitation')
    if persist:
        persist(row)
    return row


class ResidentGraph:
    def __init__(self, stock, model_dir, device):
        from rdkit import Chem
        from syntheseus import Molecule
        from syntheseus.reaction_prediction.inference.graph2edits import Graph2EditsModel
        from syntheseus.reaction_prediction.utils.misc import set_random_seed
        from syntheseus.search.algorithms.best_first.retro_star import RetroStarSearch
        from syntheseus.search.node_evaluation.common import ConstantNodeEvaluator, ReactionModelLogProbCost
        from syntheseus.search.analysis.solution_time import get_first_solution_time
        from lfr_syn.shared_stock import InchiKeyInventory
        self.Chem, self.Molecule, self.seed = Chem, Molecule, set_random_seed
        self.Search, self.Constant, self.Cost = RetroStarSearch, ConstantNodeEvaluator, ReactionModelLogProbCost
        self.solution_time, self.Inventory = get_first_solution_time, InchiKeyInventory
        begin = time.monotonic()
        self.stock = InchiKeyInventory.load_from_file(stock)
        self.model = Graph2EditsModel(model_dir=str(model_dir), device=device,
                                     use_cache=True, default_num_results=20)
        self.load_seconds = time.monotonic() - begin
        self.device = device

    def search(self, request):
        setup_started = time.monotonic()
        Chem = self.Chem
        mol = Chem.MolFromSmiles(request['smiles'])
        if mol is None or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError('Invalid candidate passed to Graph search')
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        if canonical != request['smiles']:
            raise ValueError('Candidate is not canonical')
        key = Chem.MolToInchiKey(mol)
        if key in self.stock.stock_inchikeys:
            raise ValueError('Candidate is in stock and should have been excluded before scoring')
        inventory = self.Inventory(self.stock.stock_inchikeys - {key})
        self.seed(request['seed'])
        alg = self.Search(reaction_model=self.model, mol_inventory=inventory,
            time_limit_s=request['soft_seconds'], limit_reaction_model_calls=request['calls_limit'],
            limit_iterations=10000, max_expansion_depth=12, prevent_repeat_mol_in_trees=True,
            stop_on_first_solution=True, value_function=self.Constant(constant=0.),
            and_node_cost_fn=self.Cost())
        begin = time.monotonic()
        setup_seconds = begin - setup_started
        graph = None
        reset_done = False
        try:
            alg.reset()
            reset_done = True
            alg.time_limit_s = max(0., request['soft_seconds'] - (time.monotonic() - begin))
            # Holding the graph locally preserves a partial witness on interruption.
            # This is the installed run_from_mol implementation expanded verbatim.
            graph = alg.create_graph(self.Molecule(canonical))
            alg.run_from_graph(graph)
            elapsed = time.monotonic() - begin
            solved = bool(graph.root_node.has_solution)
            for node in graph.nodes():
                node.data['analysis_time'] = node.data['num_calls_rxn_model']
            first = self.solution_time(graph) if solved else None
            if solved and (not math.isfinite(float(first)) or not 0 <= first <= self.model.num_calls()):
                raise ValueError('Invalid first-solution call counter')
            return {'status': 'solved' if solved else 'unsolved_within_budget',
                    'solved': solved, 'model_calls': int(self.model.num_calls()),
                    'search_seconds': elapsed, 'query_setup_seconds': setup_seconds,
                    'first_solution_calls': first,
                    'graph_nodes': len(graph), 'target_in_stock': False}, graph
        except BaseException as exc:
            # The partial graph may be available even when run_from_mol raised.
            elapsed = time.monotonic() - begin
            return {'status': 'execution_interrupted' if not isinstance(exc, Exception) else 'execution_error',
                    'solved': False, 'model_calls': int(self.model.num_calls()) if reset_done else 0,
                    'search_seconds': elapsed, 'query_setup_seconds': setup_seconds,
                    'error': traceback.format_exc(),
                    'interrupted': not isinstance(exc, Exception)}, graph

    def audit(self, graph, smiles):
        import networkx as nx
        from syntheseus.search.analysis.route_extraction import iter_routes_time_order
        from syntheseus.search.graph.and_or import AndNode
        Chem = self.Chem
        begin = time.monotonic()
        graph.assert_validity()
        nodes = next(iter_routes_time_order(graph, max_routes=1, max_time_s=30), None)
        if nodes is None:
            raise ValueError('no_native_route_extracted')
        synthesis = graph.to_synthesis_graph(nodes)
        synthesis.assert_validity()
        topology = nx.DiGraph()
        topology.add_nodes_from(nodes)
        topology.add_edges_from((n, c) for n in nodes for c in graph.successors(n) if c in nodes)
        if not nx.is_directed_acyclic_graph(topology) or not nx.is_weakly_connected(topology):
            raise ValueError('cyclic_or_disconnected_route')
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
        if Chem.MolToSmiles(Chem.MolFromSmiles(synthesis.root_mol.smiles), isomericSmiles=True) != canonical:
            raise ValueError('target_identity_mismatch')
        key = Chem.MolToInchiKey(Chem.MolFromSmiles(canonical))
        leaves = sorted(m.smiles for m in synthesis.get_starting_molecules())
        if not leaves or any(not self.stock.is_smiles_purchasable(s)
                             or Chem.MolToInchiKey(Chem.MolFromSmiles(s)) == key for s in leaves):
            raise ValueError('unstocked_or_target_leaf')
        depth = nx.dag_longest_path_length(topology)
        if depth > 12:
            raise ValueError('route_depth_exceeded')
        return {'verified_saved_graph': True, 'target': canonical, 'leaves': leaves,
                'alternating_depth': depth, 'seconds': time.monotonic() - begin,
                'steps': [{'product': n.reaction.product.smiles,
                           'reactants': [m.smiles for m in n.reaction.reactants]}
                          for n in sorted((n for n in nodes if isinstance(n, AndNode)),
                                          key=lambda n: n.reaction.reaction_smiles)]}


class QueryStore:
    """Exact chemical/settings deduplication; first saved seed, no old12 truncation."""
    def __init__(self, output, config, history_root, backend_factory):
        self.output, self.config = Path(output), config
        self.backend_factory, self.backend = backend_factory, None
        self.output.mkdir(parents=True, exist_ok=True)
        cfg = self.output / 'config.json'
        if cfg.exists() and read(cfg) != config:
            raise ValueError('Existing output settings differ')
        save(cfg, config)
        self.queries = self.output / 'queries'
        self.queries.mkdir(exist_ok=True)
        self.cache, self.new_calls_this_invocation = {}, 0
        self.actual_search_seconds_this_invocation = 0.
        for path in sorted(self.queries.glob('query-*.result.json')):
            record = read(path)
            key = self.key(record['request'])
            if key in self.cache:
                raise ValueError('Duplicate stored query key')
            self.cache[key] = (path, record)
        self.next_row = max([read(p)['row'] for p in self.queries.glob('query-*.result.json')] + [-1]) + 1
        hp = self.output / 'historical12_index.json'
        if hp.exists():
            data = read(hp)
            if data['settings'] != config['compatible12_settings']:
                raise ValueError('History settings differ')
            self.history = data['index']
        else:
            self.history = self.scan_history(Path(history_root))
            save(hp, {'settings': config['compatible12_settings'], 'index': self.history,
                      'selection': 'lexicographically first exact-compatible terminal record; never ranked by solve',
                      'created_unix': time.time()})

    @staticmethod
    def key(request):
        return (request['smiles'], int(request['seed']), int(request['calls_limit']), float(request['soft_seconds']))

    @staticmethod
    def chemical_key(request):
        return (request['smiles'], int(request['calls_limit']), float(request['soft_seconds']))

    def first_current_query(self, request):
        wanted = self.chemical_key(request)
        matches = [(path, record) for path, record in self.cache.values()
                   if self.chemical_key(record['request']) == wanted]
        return min(matches, key=lambda pair: str(pair[0])) if matches else None

    def scan_history(self, root):
        index, configs = {}, {}
        wanted = self.config['compatible12_settings']
        for path in sorted(root.rglob('*.result.json')):
            if self.output in path.parents or 'policy_graph' in path.parts:
                continue
            cfg = path.parent / 'config.json'
            try:
                if cfg not in configs:
                    configs[cfg] = read(cfg) if cfg.exists() else {}
                current = configs[cfg]
                if any(current.get(k) != v for k, v in wanted.items()):
                    continue
                record = read(path)
                if (record.get('status') not in NORMAL or 'canonical_smiles' not in record
                    or 'seed' not in record or 'old_model_calls' in record or 'new_model_calls' in record
                    or record.get('solved') != (record['status'] == 'solved')
                    or not path.with_name(path.name.replace('.result.json', '.graph.pkl')).exists()):
                    continue
                elapsed = record.get('search_seconds', record.get('total_seconds'))
                if (not 0 <= int(record['model_calls']) <= 12 or elapsed is None
                    or not math.isfinite(float(elapsed)) or elapsed < 0):
                    continue
                index.setdefault(record['canonical_smiles'], []).append(str(path))
            except (OSError, ValueError, TypeError, KeyError):
                continue
        return {s: sorted(paths) for s, paths in index.items()}

    @staticmethod
    def historical_audit(path, old):
        """Check the lex-selected witness; never switch to a more favourable path."""
        path = Path(path)
        audit_path = path.parent / 'graph_route_audit.json'
        if not audit_path.exists():
            return None
        report = read(audit_path)
        matches = [r for r in report['records']
                   if r.get('row') == old['row'] and r.get('directory') == str(path.parent)]
        if not matches:
            return None
        if len(matches) != 1:
            raise ValueError('Preselected history lacks a unique corresponding audit row')
        row = matches[0]
        if (row.get('search_status') != old['status']
            or row.get('reported_solved') is not old['solved']
            or row.get('verified_saved_graph') is not old['solved']
            or row.get('audit_error')
            or (old['solved'] and row.get('target') != old['canonical_smiles'])):
            raise ValueError('Preselected historical audit failed or disagrees with saved query')
        return {'path': str(audit_path), 'record': row}

    def ensure_backend(self):
        if self.backend is None:
            self.backend = self.backend_factory()
        return self.backend

    def finish_audit(self, path, record):
        graph_path = path.with_name(path.name.replace('.result.json', '.graph.pkl'))
        if record.get('status') not in NORMAL:
            return record
        if record.get('solved'):
            try:
                reading_started = time.monotonic()
                with graph_path.open('rb') as stream:
                    graph = pickle.load(stream)
                record['audit_graph_read_seconds'] = time.monotonic() - reading_started
                backend = self.ensure_backend()
                audit_started = time.monotonic()
                try:
                    record['audit'] = backend.audit(graph, record['request']['smiles'])
                finally:
                    record['audit_seconds'] = time.monotonic() - audit_started
                record['verified_solved'] = bool(record['audit']['verified_saved_graph'])
            except BaseException as exc:
                record.update(status='route_audit_error', verified_solved=False,
                              audit_error=traceback.format_exc(), interrupted=not isinstance(exc, Exception))
        else:
            record['audit'] = {'reported_solved': False, 'verified_saved_graph': False, 'seconds': 0.}
            record['audit_seconds'] = 0.
            record['verified_solved'] = False
        record['audit_complete'] = True
        if record.get('history_path') and record.get('historical_audit_missing') and record['status'] in NORMAL:
            record['historical_audit_completed_from_saved_graph'] = True
        save(path, record)
        return record

    def __call__(self, smiles, seed, calls, seconds):
        begin = time.monotonic()
        requested = {'smiles': smiles, 'seed': int(seed), 'calls_limit': int(calls), 'soft_seconds': float(seconds)}
        historical = self.history.get(smiles, []) if calls == 12 and seconds == 45. else []
        reuse = historical[0] if historical else None
        history_error = None
        try:
            old = read(reuse) if reuse else None
            if old and (old.get('canonical_smiles') != smiles or old.get('status') not in NORMAL):
                raise ValueError('Preselected historical identity/status changed')
        except (OSError, ValueError, TypeError) as exc:
            old = None
            history_error = repr(exc)
        request = dict(requested, seed=int(old['seed'])) if old else requested
        current = self.first_current_query(request)
        if current:
            # Stable path order includes unsolved, error and incomplete records.
            # A different source's requested seed never triggers duplicate work.
            path, record = current
            request = record['request']
            key = self.key(request)
            if record.get('status') == 'search_complete':
                record['status'] = record['search_status']
                record = self.finish_audit(path, record)
                self.cache[key] = (path, record)
            elif record.get('status') == 'running':
                # Never infer scientific failure or silently repeat a partially observed search.
                record = dict(record, status='saved_incomplete_query', verified_solved=False,
                              cost_observation_complete=False,
                              cost_note='An unobserved process interruption has unknown spent cost; saved initial counters are lower bounds, not measured zero work')
            return dict(record, cache_mode='current_study', actual_new_calls=0,
                        requested_seed=seed, adopted_seed=request['seed'],
                        current_study_seed_adopted=int(seed) != int(request['seed']),
                        current_study_selected_query_path=str(path),
                        retrieval_seconds=time.monotonic() - begin)
        key = self.key(request)
        row = self.next_row
        self.next_row += 1
        path = self.queries / f'query-{row:05d}.result.json'
        record = {'row': row, 'request': request, 'requested_seed': seed, 'adopted_seed': request['seed'],
                  'current_study_seed_adopted': False,
                  'status': 'running', 'started_unix': time.time(), 'pid': os.getpid(),
                  'model_calls': 0, 'search_seconds': 0., 'verified_solved': False,
                  'cost_observation_complete': False,
                  'history_path': reuse, 'result_path': str(path), 'cache_mode': 'historical12' if reuse else 'new'}
        save(path, record)
        self.cache[key] = (path, record)
        graph_path = path.with_name(path.name.replace('.result.json', '.graph.pkl'))
        if history_error:
            record.update(status='history_read_error', error=history_error, solved=False, history_unavailable=True)
            save(path, record)
            self.cache[key] = (path, record)
            return dict(record, actual_new_calls=0, retrieval_seconds=time.monotonic() - begin)
        if old:
            # Keep original metadata; audit its saved witness here, without modifying history.
            record.update(status=old['status'], solved=bool(old['solved']),
                          model_calls=int(old['model_calls']),
                          search_seconds=float(old.get('search_seconds', old.get('total_seconds', 0.))),
                          legacy_time_includes_nonsearch_overhead='search_seconds' not in old,
                          first_solution_calls=old.get('first_solution_calls'), original_result=old)
            record['cost_observation_complete'] = True
            try:
                record['historical_audit'] = self.historical_audit(reuse, old)
                record['historical_audit_missing'] = record['historical_audit'] is None
            except (OSError, ValueError, TypeError, KeyError):
                record.update(status='historical_audit_incompatible',
                              error=traceback.format_exc(), history_unavailable=True)
                save(path, record)
                self.cache[key] = (path, record)
                return dict(record, actual_new_calls=0, retrieval_seconds=time.monotonic() - begin)
            original_graph = Path(reuse).with_name(Path(reuse).name.replace('.result.json', '.graph.pkl'))
            reading_started = time.monotonic()
            try:
                with original_graph.open('rb') as stream:
                    graph = pickle.load(stream)
            except BaseException as exc:
                record.update(status='history_graph_read_error', error=traceback.format_exc(), solved=False,
                              interrupted=not isinstance(exc, Exception))
                graph = None
            finally:
                record['history_graph_read_seconds'] = time.monotonic() - reading_started
        else:
            try:
                result, graph = self.ensure_backend().search(request)
                record.update(result)
            except BaseException as exc:
                record.update(status='execution_error', error=traceback.format_exc(), solved=False,
                              setup_failed=True, interrupted=not isinstance(exc, Exception))
                graph = None
            self.new_calls_this_invocation += record['model_calls']
            self.actual_search_seconds_this_invocation += record['search_seconds']
            record['cost_observation_complete'] = True
        saving_started = time.monotonic()
        if graph is not None:
            try:
                with graph_path.open('wb') as stream:
                    pickle.dump(graph, stream)
                record['graph_path'] = str(graph_path)
            except BaseException as exc:
                record.update(status='graph_save_error', verified_solved=False, save_error=traceback.format_exc(),
                              interrupted=not isinstance(exc, Exception))
        record['saving_seconds'] = time.monotonic() - saving_started
        if record['status'] in NORMAL:
            record['search_status'] = record['status']
            record['status'] = 'search_complete'
        save(path, record)
        if record['status'] == 'search_complete':
            record['status'] = record['search_status']
            record = self.finish_audit(path, record)
        self.cache[key] = (path, record)
        returned = dict(record, result_path=str(path), actual_new_calls=0 if old else record['model_calls'],
                        retrieval_seconds=time.monotonic() - begin)
        if record.get('interrupted'):
            raise KeyboardInterrupt('Interrupted search preserved at ' + str(path))
        return returned


def validate_sources(sources, design):
    if len(sources) != design['raw_denominator']:
        raise ValueError('Need all512 original slots, not an outcome-selected subset')
    expected = {(r['task'], r['replicate'], i) for r in design['runs'] for i in range(64)}
    actual = {(s['task'], s['replicate'], s['position']) for s in sources}
    if expected != actual or len({s['sample_id'] for s in sources}) != len(sources):
        raise ValueError('Incomplete or duplicate source coordinate coverage')
    for s in sources:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', s['sample_id']):
            raise ValueError('Unsafe sample_id file name')
        if bool(s['valid']) != bool(s.get('canonical_smiles')):
            raise ValueError('Source validity/canonical identity disagreement')
        if s['cohort_id'] != f"{s['task']}-rep{s['replicate']}":
            raise ValueError('Source cohort identity disagrees with task/replicate coordinate')
        if s['valid'] and s.get('graph_result') is None:
            raise ValueError('Parent Graph record missing; do not infer route failure')


def stage_index(path, sources):
    rows = read(path)['sources']
    indexed = {r['source']['sample_id']: r for r in rows}
    if len(indexed) != len(rows) or set(indexed) != {s['sample_id'] for s in sources}:
        raise ValueError('Candidate stage must retain every source slot')
    for s in sources:
        if identity(indexed[s['sample_id']]['source']) != identity(s):
            raise ValueError('Candidate stage source identity differs')
    return indexed


def summary(rows, policies):
    result = {}
    for policy in policies:
        chosen = [r for r in rows if r['policy'] == policy]
        result[policy] = {'raw_n': len(chosen), 'solved': sum(r['solved'] is True for r in chosen),
            'joint': sum(r['joint'] is True for r in chosen), 'replaced': sum(r['replaced'] for r in chosen),
            'unique_joint_task_structures': len({(r['task'], r['final_smiles']) for r in chosen if r['joint']}),
            'execution_or_pending_limits': sum(not r['scientific_endpoint_available'] for r in chosen),
            'cost_observation_incomplete_rows': sum(not r.get('logical_cost_observation_complete', True) for r in chosen),
            'current_study_seed_adoptions': sum(a['query'].get('current_study_seed_adopted', False)
                                               for r in chosen for a in r['attempts']),
            'joint_lower': sum(r['joint'] is True or r.get('endpoint_lower_joint', False) for r in chosen),
            'joint_upper': sum(r['joint'] is True or not r['scientific_endpoint_available'] for r in chosen),
            'route_lower': sum(r['solved'] is True or r.get('endpoint_lower_solved', False) for r in chosen),
            'route_upper': sum(r['solved'] is True or not r['scientific_endpoint_available'] for r in chosen),
            'logical_calls': sum(r['logical_model_calls'] for r in chosen),
            'logical_search_seconds': sum(r['logical_search_seconds'] for r in chosen)}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, default=STUDY)
    parser.add_argument('--sources', type=Path)
    parser.add_argument('--candidate-stage', type=Path)
    parser.add_argument('--native', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cohort', action='append', help='Whole cohort_id; repeat for multiple complete cohorts')
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--stock', type=Path, default=Path(STOCK))
    parser.add_argument('--model-dir', type=Path, default=Path(MODEL))
    parser.add_argument('--history-root', type=Path, default=BASE)
    args = parser.parse_args()
    def interrupted_signal(signum, frame):
        raise KeyboardInterrupt('Termination signal ' + str(signum))
    signal.signal(signal.SIGTERM, interrupted_signal)
    os.environ['OMP_NUM_THREADS'] = '2'
    os.environ['OPENBLAS_NUM_THREADS'] = '2'
    os.environ['MKL_NUM_THREADS'] = '2'
    design = read(args.study / 'design.json')
    sources = read(args.sources or args.study / 'sources.json')['rows']
    validate_sources(sources, design)
    candidates = stage_index(args.candidate_stage or args.study / 'candidate_stage.json', sources)
    native = stage_index(args.native, sources) if args.native else None
    cohorts = [f"{r['task']}-rep{r['replicate']}" for r in design['runs']]
    if {s['cohort_id'] for s in sources} != set(cohorts):
        raise ValueError('Cohort IDs differ from complete design')
    selected = set(args.cohort) if args.cohort else set(cohorts)
    if not selected <= set(cohorts):
        raise ValueError('Unknown cohort')
    output = args.output or args.study / 'policy_graph'
    versions = {p: importlib.metadata.version(p) for p in ('syntheseus', 'torch', 'rdkit')}
    compatible = dict(COMPATIBLE, stock=str(args.stock), model_dir=str(args.model_dir),
                      reaction_model_calls_limit=12, time_limit_s=45, versions=versions)
    config = {'scope': design['scope'], 'study': str(args.study), 'versions': versions,
        'stock': str(args.stock), 'model_dir': str(args.model_dir), 'device': args.device,
        'compatible12_settings': compatible, 'policies': list(POLICIES),
        'shared_soft_seconds': 45, 'native_enumeration_not_in_planning_budget': True,
        'time_note': 'search_seconds after resident model; audit/save/cache-retrieval time separately reported; cached query search duration counts logically',
        'candidate_seed_rule': 'First request uses source.search_seed unless historical12 adopts the lex-selected saved seed. Within this study, identical smiles/calls_limit/exact soft_seconds adopts the recorded seed of the lex-first query path, including unsolved/error/incomplete queries, without replacing original query metadata. Requested and adopted seeds are reported; new training replication comes from eight checkpoints, not search seeds.',
        'current_study_reuse_selection': 'lexicographically first query path for exact smiles/calls_limit/soft_seconds; never ranked by status or solved; distinct second-branch remaining times are distinct settings',
        'cohort_policy_orders': {c: ordered_policies(i) for i, c in enumerate(cohorts)}}
    store = QueryStore(output, config, args.history_root,
                       lambda: ResidentGraph(str(args.stock), str(args.model_dir), args.device))
    state_path = output / 'execution.json'
    started = time.time()
    state = {'pid': os.getpid(), 'status': 'running', 'started_unix': started,
             'cohorts_requested': sorted(selected), 'native_stage_supplied': bool(native)}
    save(state_path, state)
    rows = []
    try:
        for ci, cohort in enumerate(cohorts):
            if cohort not in selected:
                continue
            for source in sorted((s for s in sources if s['cohort_id'] == cohort), key=lambda s: s['position']):
                for policy in ordered_policies(ci):
                    path = output / 'outcomes' / policy / (source['sample_id'] + '.json')
                    previous = read(path) if path.exists() else None
                    if previous and previous.get('status') == 'pending_candidate_stage' and native:
                        previous = None
                    stage = native if policy.startswith('SynTwins') else candidates
                    arm = None
                    if stage:
                        entry = stage[source['sample_id']]
                        if ARMS[policy] in entry.get('arms', {}):
                            arm = dict(entry['arms'][ARMS[policy]])
                            if policy.startswith('SynTwins'):
                                arm['native_status'] = entry.get('native_status', arm.get('native_status'))
                        elif policy.startswith('SynTwins') and 'scored_candidates' in entry:
                            arm = dict(entry)
                    row = evaluate_policy(source, policy, arm, store, saved=previous,
                                          persist=lambda r, p=path: save(p, r))
                    rows.append(row)
                save(output / 'progress.json', {'last_source': source['sample_id'],
                    'evaluated_policy_rows': len(rows), 'actual_new_calls_this_invocation': store.new_calls_this_invocation})
                print(json.dumps({'source': source['sample_id'], 'policy_rows': len(rows),
                                  'new_calls': store.new_calls_this_invocation}), flush=True)
        report = {'scope': design['scope'], 'full_design_raw_denominator': len(sources),
            'selected_complete_cohorts': sorted(selected), 'rows': rows,
            'policies': summary(rows, POLICIES), 'new_calls_this_invocation': store.new_calls_this_invocation,
            'new_search_seconds_this_invocation': store.actual_search_seconds_this_invocation,
            'stored_actual_new_calls_cumulative': sum(r.get('model_calls', 0) for _, r in store.cache.values() if r.get('cache_mode') == 'new'),
            'stored_cost_observation_incomplete_queries': sum(not r.get('cost_observation_complete', True) for _, r in store.cache.values()),
            'stored_historical12_queries': sum(r.get('history_path') is not None for _, r in store.cache.values()),
            'stored_saving_seconds_cumulative': sum(r.get('saving_seconds', 0.) for _, r in store.cache.values()),
            'stored_query_setup_seconds_cumulative': sum(r.get('query_setup_seconds', 0.) for _, r in store.cache.values()),
            'stored_route_audit_seconds_cumulative': sum(r.get('audit_seconds', 0.) for _, r in store.cache.values()),
            'stored_graph_read_seconds_cumulative': sum(r.get('audit_graph_read_seconds', 0.) + r.get('history_graph_read_seconds', 0.) for _, r in store.cache.values()),
            'resident_model_load_seconds_this_invocation': store.backend.load_seconds if store.backend else 0.,
            'native_pending': args.native is None, 'completed_unix': time.time(),
            'interpretation': 'Shared chemical measurements and historical12 reuse; logical policy budgets are not isolated deployment wallclock. Execution-limited endpoints remain unknown.'}
        save(output / 'policy_results.json', report)
        has_limits = any(not r['scientific_endpoint_available'] and r['status'] != 'pending_candidate_stage' for r in rows)
        status = 'complete_internal_native_pending' if args.native is None else 'complete'
        if has_limits:
            status += '_with_execution_limits'
        state.update(status=status, ended_unix=time.time(), actual_new_calls=store.new_calls_this_invocation)
        save(state_path, state)
        return int(has_limits)
    except BaseException:
        state.update(status='interrupted_or_failed', error=traceback.format_exc(), ended_unix=time.time(),
                     actual_new_calls=store.new_calls_this_invocation)
        save(state_path, state)
        raise


if __name__ == '__main__':
    raise SystemExit(main())
