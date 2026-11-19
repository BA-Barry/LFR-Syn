"""Orchestration/persistence tests using injected searches; no model or GPU work."""
import copy
import pickle
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from lfr_method import graph_worker as worker
from lfr_method.graph_core import ResidentGraph, read, save

VERSIONS = {"syntheseus": "test", "torch": "test", "rdkit": "test"}


def source(sid="gsk3b-rep0-r005", valid=True):
    return dict(sample_id=sid, task="gsk3b", task_index=0, replicate=0,
                position=5, cohort_id="gsk3b-rep0", canonical_smiles="parent" if valid else None,
                valid=valid, task_score=.02 if valid else 0., qed=.6 if valid else None,
                search_seed=10302005, graph_verified=False, graph_result=None,
                graph_result_path=None)


def arm():
    return {"scored_candidates": [dict(smiles="first", task_score=.3, qed=.6, task_pass=True),
                                  dict(smiles="second", task_score=.29, qed=.6, task_pass=True)]}


class FakeBackend:
    load_seconds = .01

    def __init__(self, responses, audit_error=False):
        self.responses, self.requests, self.audited = list(responses), [], []
        self.audit_error = audit_error

    def search(self, request):
        self.requests.append(copy.deepcopy(request))
        response = self.responses.pop(0)
        status, calls, seconds = response[:3]
        result = {"status": status, "solved": status == "solved", "model_calls": calls,
                  "search_seconds": seconds, "cost_observation_complete": True}
        if len(response) > 3:
            result["interrupted"] = response[3]
        return result, {"saved_marker": request["smiles"], "status": status}

    def audit(self, graph, smiles):
        self.audited.append((graph, smiles))
        if self.audit_error:
            raise ValueError("synthetic witness failure")
        if graph != {"saved_marker": smiles, "status": "solved"}:
            raise ValueError("Graph is not the recorded saved witness")
        return {"verified_saved_graph": True, "target": smiles, "leaves": ["stock-leaf"]}


class WorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.output = Path(self.temp.name) / "run"
        self.output.mkdir()
        self.config = dict(scope="software_test_only", package_version="0.1.0",
                           rf_model_dir="models", graph_model_dir="graph-model", stock_keys="stock.txt",
                           rf_python="configured-rf-python", graph_python="configured-graph-python",
                           forward_templates="forward.csv", retro_templates="retro.csv", stock_smiles="stock.csv",
                           device="cuda:0", reuse_study=None, parent_calls=12, parent_soft_seconds=45.,
                           candidate_calls=[6, 6], candidate_shared_soft_seconds=45.,
                           candidate_score_slots_per_arm=8, candidate_forward_attempts_per_arm=16,
                           policies=list(worker.POLICIES))
        self.setup_rows([source()])

    def setup_rows(self, rows):
        self.rows = rows
        save(self.output / "method_config.json", self.config)
        save(self.output / "selected_inputs.json", {"rows": rows})
        save(self.output / "sources_scored.json", {"rows": rows, "supplied_inputs": rows,
             "rf_model_dir": "models", "actual_new_source_predictions": sum(s["valid"] for s in rows)})

    def candidate_runner(self, output, config):
        parents = read(output / "sources_with_parents.json")["rows"]
        save(output / "candidate_stage.json", {"sources": [
            {"source": s, "arms": {name: arm() for name in worker.ARMS.values()}} for s in parents],
            "actual_new_candidate_predictions": 2, "actual_new_forward_attempts": 4})

    def run_backend(self, backend, candidate_runner=None):
        return worker.run(self.output, backend_factory=lambda: backend,
                          candidate_runner=candidate_runner or self.candidate_runner,
                          version_provider=lambda: VERSIONS)

    def test_real_orchestration_roles_shared_time_and_saved_witness(self):
        backend = FakeBackend([("unsolved_within_budget", 12, 21.),
                               ("unsolved_within_budget", 6, 17.), ("solved", 2, 4.)])
        self.assertEqual(self.run_backend(backend), 0)
        self.assertEqual([r["target_role"] for r in backend.requests], ["parent", "candidate", "candidate"])
        self.assertEqual([r["calls_limit"] for r in backend.requests], [12, 6, 6])
        self.assertEqual([r["soft_seconds"] for r in backend.requests], [45., 45., 28.])
        result = read(self.output / "method_result.json")
        self.assertEqual(result["mode"], "real_method_execution")
        self.assertEqual(result["actual_new_work"]["Graph_queries"], 3)
        self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 20)
        self.assertEqual(result["actual_new_work"]["RF_predictions"], 3)
        self.assertEqual(len(backend.audited), 1)
        for policy in worker.POLICIES:
            self.assertEqual(result["policies"][policy]["joint"], 1)
            self.assertEqual(result["policies"][policy]["logical_calls"], 8)
        graph_path = next(p for p in (self.output / "graph_measurements/queries").glob("*.graph.pkl")
                          if pickle.loads(p.read_bytes())["status"] == "solved")
        self.assertEqual(pickle.loads(graph_path.read_bytes())["saved_marker"], "second")

    def test_completed_resume_does_not_repeat_search_or_remove_actual_cost(self):
        backend = FakeBackend([("unsolved_within_budget", 12, 2.), ("solved", 1, 1.)])
        self.run_backend(backend)
        never = FakeBackend([])
        self.assertEqual(self.run_backend(never), 0)
        self.assertFalse(never.requests)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["new_graph_calls_this_invocation"], 0)
        self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 13)

    def test_full_cpu_costs_survive_completed_resume_from_legacy_stage_records(self):
        def saved_cpu(output, config):
            path = output / "candidate_stage.json"
            if path.exists():
                return
            self.candidate_runner(output, config)
            stage = read(path)
            record = stage["sources"][0]
            record.update(status="candidate_scoring_complete", costs={
                "source_replay_forward_calls": 121, "retro_template_calls": 1071,
                "source_RF_path_attributions": 1})
            record["arms"]["task_contribution"]["costs"] = {"new_forward_attempts": 6}
            record["arms"]["permuted_contribution"]["costs"] = {"new_forward_attempts": 7}
            stage["actual_new_forward_attempts"] = 13
            # No new CPU summary is present: old complete source records must
            # still supply cumulative component costs without rerunning RF.
            save(path, stage)
        backend = FakeBackend([("unsolved_within_budget", 12, 2.), ("solved", 1, 1.)])
        self.run_backend(backend, saved_cpu)
        first = read(self.output / "method_result.json")
        never = FakeBackend([])
        self.run_backend(never, saved_cpu)
        resumed = read(self.output / "method_result.json")
        self.assertFalse(never.requests)
        self.assertEqual(first["actual_new_cpu_work"], resumed["actual_new_cpu_work"])
        self.assertEqual(resumed["actual_new_work"]["forward_reactions"], 13)
        self.assertEqual(resumed["actual_new_work"]["source_replay_forward_reactions"], 121)
        self.assertEqual(resumed["actual_new_work"]["candidate_forward_reactions"], 13)
        self.assertEqual(resumed["actual_new_work"]["forward_reactions_total"], 134)
        self.assertEqual(resumed["actual_new_work"]["retro_template_calls"], 1071)
        self.assertEqual(resumed["actual_new_work"]["source_RF_path_attributions"], 1)

    def test_missing_legacy_cpu_counters_do_not_become_a_zero_complete_total(self):
        backend = FakeBackend([("unsolved_within_budget", 12, 2.), ("solved", 1, 1.)])
        self.run_backend(backend)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["actual_new_work"]["forward_reactions"], 4)
        self.assertIsNone(result["actual_new_work"]["source_replay_forward_reactions"])
        self.assertIsNone(result["actual_new_work"]["forward_reactions_total"])
        self.assertFalse(result["actual_new_cpu_work"]["cost_observation_complete"])

    def test_invalid_raw_slot_retained_without_search(self):
        self.setup_rows([source(valid=False)])
        backend = FakeBackend([])
        self.assertEqual(self.run_backend(backend), 0)
        self.assertFalse(backend.requests)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["selected_raw_n"], 1)
        self.assertTrue(all(p["raw_n"] == 1 and p["joint"] == 0 for p in result["policies"].values()))

    def test_parent_joint_success_skips_candidate_search(self):
        s = source()
        s["task_score"] = .2
        self.setup_rows([s])
        backend = FakeBackend([("solved", 3, 2.)])
        self.assertEqual(self.run_backend(backend), 0)
        result = read(self.output / "method_result.json")
        self.assertEqual(len(backend.requests), 1)
        self.assertEqual(result["parent_verified_routes"], 1)
        self.assertTrue(all(r["status"] == "source_retained" for r in result["rows"]))

    def test_parent_error_is_unknown_and_never_retried(self):
        backend = FakeBackend([("execution_error", 2, 4.)])
        self.assertEqual(self.run_backend(backend), 1)
        self.assertEqual(self.run_backend(FakeBackend([])), 1)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["parent_unknown_raw_slots"], 1)
        self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 2)
        self.assertTrue(all(r["joint"] is None and r["status"] == "parent_execution_limitation"
                            for r in result["rows"]))

    def test_parent_call_overrun_is_unknown_preserving_observed_spend(self):
        backend = FakeBackend([("solved", 13, 4.)])
        self.assertEqual(self.run_backend(backend), 1)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["parent_verified_routes"], 0)
        self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 13)
        parent = read(self.output / "sources_with_parents.json")["rows"][0]
        self.assertEqual(parent["graph_result"]["status"], "reaction_call_budget_limitation")
        self.assertEqual(parent["reported_parent_result"]["status"], "solved")

    def test_candidate_audit_error_cannot_be_a_positive_route(self):
        backend = FakeBackend([("unsolved_within_budget", 12, 2.), ("solved", 2, 3.)], audit_error=True)
        self.assertEqual(self.run_backend(backend), 1)
        result = read(self.output / "method_result.json")
        self.assertTrue(all(r["joint"] is None and r["status"] == "search_execution_limitation"
                            for r in result["rows"]))
        self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 14)

    def test_saved_interruption_is_not_researched_on_resume(self):
        backend = FakeBackend([("execution_interrupted", 2, 3., True)])
        with self.assertRaises(KeyboardInterrupt):
            self.run_backend(backend)
        self.assertEqual(read(self.output / "method_execution.json")["status"], "interrupted_or_failed")
        self.assertEqual(self.run_backend(FakeBackend([])), 1)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["parent_unknown_raw_slots"], 1)
        self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 2)

    def test_no_history_reuse_is_the_default(self):
        unrelated = Path(self.temp.name) / "old" / "queries"
        unrelated.mkdir(parents=True)
        save(unrelated / "config.json", worker.graph_config(self.config, VERSIONS)["compatible12_settings"])
        save(unrelated / "row-000.result.json", {"row": 0, "canonical_smiles": "parent", "seed": 1,
             "status": "solved", "solved": True, "model_calls": 1, "search_seconds": 1.})
        (unrelated / "row-000.graph.pkl").write_bytes(pickle.dumps({"saved_marker": "parent", "status": "solved"}))
        backend = FakeBackend([("unsolved_within_budget", 12, 2.), ("solved", 1, 1.)])
        self.run_backend(backend)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["historical12_queries"], 0)
        self.assertEqual(backend.requests[0]["smiles"], "parent")

    def test_explicit_compatible_history_uses_first_negative_and_saved_seed(self):
        history = Path(self.temp.name) / "old" / "parent"
        history.mkdir(parents=True)
        save(history / "config.json", worker.graph_config(self.config, VERSIONS)["compatible12_settings"])
        for row, status in [(0, "unsolved_within_budget"), (1, "solved")]:
            save(history / f"row-{row:03d}.result.json", {"row": row, "canonical_smiles": "parent", "seed": 99 + row,
                 "status": status, "solved": status == "solved", "model_calls": 3, "search_seconds": 5.})
            (history / f"row-{row:03d}.graph.pkl").write_bytes(pickle.dumps({"saved_marker": "parent", "status": status}))
        self.config["reuse_study"] = str(history.parent)
        self.setup_rows(self.rows)
        backend = FakeBackend([("solved", 1, 1.)])
        self.run_backend(backend)
        parent = read(self.output / "sources_with_parents.json")["rows"][0]
        self.assertEqual(parent["search_seed"], 99)
        self.assertFalse(parent["graph_verified"])
        self.assertTrue(parent["prior_graph_reused"])
        self.assertEqual(backend.requests[0]["seed"], 99)
        self.assertEqual(backend.requests[0]["calls_limit"], 6)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["historical12_queries"], 1)
        self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 1)

    def test_candidate_identity_mismatch_stops_before_candidate_queries(self):
        def mismatch(output, config):
            self.candidate_runner(output, config)
            stage = read(output / "candidate_stage.json")
            stage["sources"][0]["source"]["search_seed"] += 1
            save(output / "candidate_stage.json", stage)
        backend = FakeBackend([("unsolved_within_budget", 12, 2.)])
        with self.assertRaisesRegex(ValueError, "Candidate stage source identity"):
            self.run_backend(backend, mismatch)
        self.assertEqual(len(backend.requests), 1)

    def test_rf_worker_uses_explicit_separate_interpreter(self):
        from unittest.mock import patch
        with patch.object(worker.subprocess, "call", return_value=0) as call:
            worker.run_rf_candidates(self.output, self.config)
        command = call.call_args.args[0]
        self.assertEqual(command[:4], ["configured-rf-python", "-m", "lfr_method.rf_worker", "candidates"])


class ResidentSearchSemanticsTests(unittest.TestCase):
    """Exercise the real retained search method against fake installed APIs."""
    def backend(self):
        model = SimpleNamespace(calls=99)
        model.num_calls = lambda: model.calls
        observations = {"inventories": [], "resets": 0}
        class Chem:
            MolFromSmiles = staticmethod(lambda smiles: {"smiles": smiles})
            GetMolFrags = staticmethod(lambda mol: ((1,),))
            MolToSmiles = staticmethod(lambda mol, **kwargs: mol["smiles"])
            MolToInchiKey = staticmethod(lambda mol: "target" if mol["smiles"] == "parent" else "other")
        class Graph:
            root_node = SimpleNamespace(has_solution=False)
            def nodes(self):
                return [SimpleNamespace(data={"num_calls_rxn_model": 3})]
            def __len__(self):
                return 1
        class Search:
            def __init__(self, **kwargs):
                self.time_limit_s = kwargs["time_limit_s"]
                self.settings = kwargs
            def reset(self):
                observations["resets"] += 1
                model.calls = 0
            def create_graph(self, molecule):
                return Graph()
            def run_from_graph(self, graph):
                model.calls += 3
        def inventory(keys):
            observations["inventories"].append(set(keys))
            return keys
        result = object.__new__(ResidentGraph)
        result.Chem, result.Molecule, result.seed = Chem, lambda s: s, lambda seed: None
        result.Search, result.Constant, result.Cost = Search, lambda **kwargs: None, lambda: None
        result.Inventory, result.stock, result.model = inventory, SimpleNamespace(stock_inchikeys={"target", "leaf"}), model
        return result, observations

    def test_parent_in_stock_allowed_but_target_excluded_and_reset_every_request(self):
        backend, seen = self.backend()
        role = worker.RoleAwareGraph(backend, {"target_role": "parent"})
        for unused in range(2):
            result, graph = role.search({"smiles": "parent", "seed": 7, "calls_limit": 12, "soft_seconds": 45.})
            self.assertEqual(result["model_calls"], 3)
            self.assertTrue(result["target_in_stock"])
        self.assertEqual(seen["inventories"], [{"leaf"}, {"leaf"}])
        self.assertEqual(seen["resets"], 2)

    def test_candidate_in_stock_rejected_before_search(self):
        backend, seen = self.backend()
        role = worker.RoleAwareGraph(backend, {"target_role": "candidate"})
        with self.assertRaisesRegex(ValueError, "Candidate is in stock"):
            role.search({"smiles": "parent", "seed": 7, "calls_limit": 6, "soft_seconds": 45.})
        self.assertFalse(seen["inventories"])
        self.assertEqual(seen["resets"], 0)


if __name__ == "__main__":
    unittest.main()
