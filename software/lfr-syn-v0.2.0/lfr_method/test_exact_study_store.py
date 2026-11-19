"""Exact reuse tests with synthetic witnesses; no RDKit, model or GPU calls.

The two included software examples supply real source/candidate identities and
saved logical costs. Test graphs/auditors are injected, so this suite does not
claim to audit their actual scientific graph witnesses.
"""
import copy
from pathlib import Path
import pickle
import tempfile
import unittest

from lfr_method import graph_worker as worker
from lfr_method.exact_study_store import ExactStudyStore, index_key
from lfr_method.graph_core import read, save
from lfr_method.test_graph_worker import FakeBackend, VERSIONS

ROOT = Path(__file__).resolve().parents[1]


class FakeAuditor:
    setup_seconds = 0.

    def __init__(self, fail=False):
        self.calls, self.fail = [], fail

    def audit(self, graph, smiles):
        self.calls.append((copy.deepcopy(graph), smiles))
        if self.fail:
            raise AssertionError("synthetic saved-witness audit failure")
        if graph != {"witness": smiles}:
            raise ValueError("Wrong original saved graph was read")
        return {"verified_saved_graph": True, "target": smiles, "leaves": ["test-stock-leaf"]}


class ExactStudyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.study = self.root / "explicit_study"
        self.study.mkdir()
        self.config = worker.graph_config(dict(scope="software_test_only", package_version="0.1.0",
            stock_keys="stock.txt", graph_model_dir="graph-model", device="cuda:0", reuse_study=str(self.study),
            parent_calls=12, parent_soft_seconds=45., candidate_calls=[6, 6], candidate_shared_soft_seconds=45.,
            candidate_score_slots_per_arm=8, candidate_forward_attempts_per_arm=16,
            policies=list(worker.POLICIES)), VERSIONS)
        self.backend, self.auditor = FakeBackend([]), FakeAuditor()

    def modern(self, name="query-00000", *, smiles="candidate", cap=6, seconds=45., seed=7,
               status="solved", calls=1, elapsed=3., overrides=None, directory="policy_graph/gpu0",
               graph=True):
        parent = self.study / directory
        query_dir = parent / "queries"
        query_dir.mkdir(parents=True, exist_ok=True)
        save(parent / "config.json", self.config)
        solved = status == "solved"
        record = {"row": int(name.split("-")[-1]),
                  "request": {"smiles": smiles, "seed": seed, "calls_limit": cap, "soft_seconds": seconds},
                  "status": status, "solved": solved, "verified_solved": solved,
                  "model_calls": calls, "search_seconds": elapsed,
                  "cost_observation_complete": status in ("solved", "unsolved_within_budget", "execution_error"),
                  "audit_complete": status in ("solved", "unsolved_within_budget"),
                  "audit": {"verified_saved_graph": solved, **({"target": smiles} if solved else {})},
                  "cache_mode": "new"}
        if overrides:
            record.update(overrides)
        path = query_dir / f"{name}.result.json"
        save(path, record)
        if graph:
            path.with_name(path.name.replace(".result.json", ".graph.pkl")).write_bytes(pickle.dumps({"witness": smiles}))
        return path, record

    def store(self, suffix="current"):
        return ExactStudyStore(self.root / suffix, self.config, self.study,
                               lambda: self.backend, lambda: self.auditor)

    def test_exact6_reloads_original_graph_and_never_constructs_search_backend(self):
        old, original = self.modern()
        def no_search_backend():
            self.fail("Search model/backend was constructed for exact reuse")
        store = ExactStudyStore(self.root / "current", self.config, self.study,
                                no_search_backend, lambda: self.auditor)
        q = store("candidate", 99, 6, 45.)
        self.assertEqual(q["cache_mode"], "explicit_study_exact")
        self.assertEqual(q["adopted_seed"], 7)
        self.assertEqual(q["actual_new_calls"], 0)
        self.assertEqual(q["old_work"]["model_calls"], 1)
        self.assertEqual(q["new_work"]["model_calls"], 0)
        self.assertTrue(q["verified_solved"])
        self.assertTrue(q["study_saved_graph_reaudited"])
        self.assertEqual(self.auditor.calls, [({"witness": "candidate"}, "candidate")])
        self.assertEqual(q["graph_path"], str(old.with_name(old.name.replace(".result.json", ".graph.pkl"))))
        self.assertEqual(q["original_result"], original)
        self.assertEqual(read(old), original)

    def test_cap12_is_not_a_cap6_prefix(self):
        self.modern(cap=12, calls=1, elapsed=1.)
        self.backend.responses = [("unsolved_within_budget", 6, 2.)]
        q = self.store()("candidate", 99, 6, 45.)
        self.assertEqual(q["cache_mode"], "new")
        self.assertEqual(len(self.backend.requests), 1)
        self.assertEqual(self.backend.requests[0]["calls_limit"], 6)
        self.assertEqual(q["actual_new_calls"], 6)
        self.assertFalse(self.auditor.calls)

    def test_second_branch_time_requires_exact_match(self):
        self.modern(seconds=29.778339915000288)
        self.backend.responses = [("unsolved_within_budget", 6, 2.)]
        q = self.store()("candidate", 99, 6, 29.77833991500029)
        self.assertEqual(q["cache_mode"], "new")
        self.assertEqual(len(self.backend.requests), 1)

    def test_lex_first_negative_is_not_replaced_by_solved(self):
        first, original = self.modern(status="unsolved_within_budget", calls=6, elapsed=17.)
        self.modern(name="query-00001", status="solved", seed=8)
        q = self.store()("candidate", 99, 6, 45.)
        self.assertEqual(q["status"], "unsolved_within_budget")
        self.assertFalse(q["verified_solved"])
        self.assertEqual(q["study_query_path"], str(first))
        self.assertEqual(q["adopted_seed"], 7)
        self.assertEqual(q["search_seconds"], 17.)
        self.assertFalse(self.backend.requests)

    def test_lex_first_unknown_is_not_skipped_for_normal_terminal(self):
        first, original = self.modern(status="execution_error", calls=2, elapsed=4.)
        self.modern(name="query-00001", status="solved", seed=8)
        store = self.store()
        q = store("candidate", 99, 6, 45.)
        self.assertEqual(q["status"], "execution_error")
        self.assertEqual(q["study_query_path"], str(first))
        self.assertFalse(q["verified_solved"])
        self.assertEqual(q["model_calls"], 2)
        self.assertEqual(q["actual_new_calls"], 0)
        again = store("candidate", 100, 6, 45.)
        self.assertEqual(again["status"], "execution_error")
        self.assertFalse(self.backend.requests)

    def test_running_cache_retains_unknown_spend_and_is_not_researched(self):
        self.modern(status="running", calls=0, elapsed=0., graph=False)
        self.modern(name="query-00001", status="solved", seed=8)
        q = self.store()("candidate", 99, 6, 45.)
        self.assertEqual(q["status"], "saved_incomplete_query")
        self.assertFalse(q["cost_observation_complete"])
        resumed = self.store()("candidate", 100, 6, 45.)
        self.assertEqual(resumed["status"], "saved_incomplete_query")
        self.assertFalse(self.backend.requests)

    def test_current_reaudit_failure_stays_unknown_without_alternate_or_search(self):
        first, original = self.modern()
        self.modern(name="query-00001", seed=8)
        self.auditor.fail = True
        store = self.store()
        q = store("candidate", 99, 6, 45.)
        self.assertEqual(q["status"], "study_reuse_limitation")
        self.assertFalse(q["verified_solved"])
        self.assertEqual(q["old_work"]["model_calls"], 1)
        self.assertEqual(q["study_query_path"], str(first))
        self.assertEqual(store("candidate", 100, 6, 45.)["status"], "study_reuse_limitation")
        self.assertEqual(len(self.auditor.calls), 1)
        self.assertFalse(self.backend.requests)

    def test_equivalent_alias_resolves_to_included_physical_origin(self):
        origin, original = self.modern(directory="z_physical/gpu0")
        alias, record = self.modern(directory="a_alias/gpu0", graph=False,
            overrides={"cache_mode": "joint512_exact", "old512_query_path": str(origin),
                       "original_result": copy.deepcopy(original)})
        q = self.store()("candidate", 99, 6, 45.)
        self.assertEqual(q["study_query_path"], str(origin))
        self.assertTrue(q["verified_solved"])
        self.assertEqual(len(self.auditor.calls), 1)

    def test_alias_with_different_state_is_not_replaced_by_favourable_origin(self):
        origin, original = self.modern(directory="z_physical/gpu0")
        alias, record = self.modern(directory="a_alias/gpu0", status="execution_error", graph=False,
            overrides={"cache_mode": "joint512_exact", "old512_query_path": str(origin),
                       "original_result": copy.deepcopy(original)})
        q = self.store()("candidate", 99, 6, 45.)
        self.assertEqual(q["study_query_path"], str(alias))
        self.assertEqual(q["status"], "execution_error")
        self.assertFalse(q["verified_solved"])
        self.assertFalse(self.auditor.calls)
        self.assertFalse(self.backend.requests)

    def test_alias_cannot_relabel_cap12_origin_as_cap6(self):
        origin, original = self.modern(directory="z_physical/gpu0", cap=12)
        alias, record = self.modern(directory="a_alias/gpu0", cap=6, graph=False,
            overrides={"cache_mode": "joint512_exact", "old512_query_path": str(origin),
                       "graph_path": str(origin.with_name(origin.name.replace(".result.json", ".graph.pkl"))),
                       "original_result": copy.deepcopy(original)})
        q = self.store()("candidate", 99, 6, 45.)
        self.assertEqual(q["study_query_path"], str(alias))
        self.assertEqual(q["status"], "study_reuse_limitation")
        self.assertFalse(q["verified_solved"])
        self.assertFalse(self.auditor.calls)
        self.assertFalse(self.backend.requests)

    def test_outside_origin_is_not_followed_and_saved_unknown_not_researched(self):
        external = self.root / "unapproved" / "query-00000.result.json"
        external.parent.mkdir()
        external.write_text("not JSON; must never be read", encoding="utf-8")
        self.modern(directory="a_alias/gpu0", graph=False,
            overrides={"cache_mode": "joint512_exact", "old512_query_path": str(external),
                       "graph_path": str(external.with_suffix(".graph.pkl"))})
        q = self.store()("candidate", 99, 6, 45.)
        self.assertEqual(q["status"], "study_reuse_limitation")
        self.assertFalse(q["verified_solved"])
        self.assertFalse(self.auditor.calls)
        self.assertFalse(self.backend.requests)

    def test_authorized_cross_directory_origin_is_read_without_scanning_peers(self):
        from unittest.mock import patch
        from lfr_method import exact_study_store as exact
        legacy = self.root / "authorized_project_development"
        origin, original = self.modern(directory=str(legacy / "prior_policy/gpu0"))
        (origin.parent / "query-00001.result.json").write_text("not JSON; must not be enumerated", encoding="utf-8")
        self.modern(directory="a_alias/gpu0", graph=False,
            overrides={"cache_mode": "joint512_exact", "old512_query_path": str(origin),
                       "original_result": copy.deepcopy(original)})
        with patch.object(exact, "PROJECT_DEV", legacy):
            store = self.store()
            q = store("candidate", 99, 6, 45.)
        self.assertEqual(q["physical_origin_path"], str(origin))
        self.assertEqual(q["study_query_path"], str(origin))
        self.assertTrue(q["verified_solved"])
        self.assertEqual(len(self.auditor.calls), 1)
        self.assertFalse(self.backend.requests)
        report = read(store.output / "explicit_study_exact_index.json")
        self.assertEqual(report["diagnostics"]["declared_external_origin_files_read"], 1)

    def test_authorized_historical12_alias_resolves_to_original_trace_parent(self):
        from unittest.mock import patch
        from lfr_method import exact_study_store as exact
        project_user = self.root / "authorized_project_user"
        parent = project_user / "trace_lfr_previous/parent_graph"
        parent.mkdir(parents=True)
        save(parent / "config.json", self.config["compatible12_settings"])
        original = {"row": 0, "canonical_smiles": "parent", "seed": 7,
                    "status": "unsolved_within_budget", "solved": False,
                    "model_calls": 12, "search_seconds": 20.}
        origin = parent / "row-000.result.json"
        save(origin, original)
        origin.with_name("row-000.graph.pkl").write_bytes(pickle.dumps({"witness": "parent"}))
        self.modern(directory="a_alias/gpu0", smiles="parent", cap=12, status="unsolved_within_budget",
            calls=12, elapsed=20., graph=False, overrides={"cache_mode": "historical12",
                "history_path": str(origin), "original_result": copy.deepcopy(original)})
        with patch.object(exact, "PROJECT_USER", project_user):
            q = self.store()("parent", 99, 12, 45.)
        self.assertEqual(q["physical_origin_path"], str(origin))
        self.assertEqual(q["status"], "unsolved_within_budget")
        self.assertEqual(q["model_calls"], 12)
        self.assertEqual(q["search_seconds"], 20.)
        self.assertEqual(q["actual_new_calls"], 0)
        self.assertFalse(self.backend.requests)

    def test_changed_preselected_configuration_does_not_trigger_new_search(self):
        path, record = self.modern()
        store = self.store()
        cfg = read(path.parent.parent / "config.json")
        cfg["shared_soft_seconds"] = 90
        save(path.parent.parent / "config.json", cfg)
        q = store("candidate", 99, 6, 45.)
        self.assertEqual(q["status"], "study_reuse_limitation")
        self.assertFalse(self.backend.requests)

    def test_auditor_does_not_call_resident_model_constructor(self):
        from unittest.mock import patch
        from lfr_method.exact_study_store import SavedGraphAuditor
        self.assertIs(SavedGraphAuditor.audit, worker.ResidentGraph.audit)
        with patch.object(worker.ResidentGraph, "__init__", side_effect=AssertionError("model constructor forbidden")):
            self.modern()
            q = self.store()("candidate", 99, 6, 45.)
        self.assertTrue(q["verified_solved"])


class ExistingExamplesRoundTripTests(unittest.TestCase):
    def test_both_included_software_examples_exact_graph_reuse_without_model_calls(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output, study = root / "run", root / "explicit_study"
            output.mkdir()
            study.mkdir()
            config = dict(scope="software_test_only", package_version="0.1.0", rf_model_dir="models",
                stock_keys="stock.txt", graph_model_dir="graph-model", device="cuda:0", reuse_study=str(study),
                rf_python="test-rf-python", graph_python="test-graph-python",
                parent_calls=12, parent_soft_seconds=45., candidate_calls=[6, 6], candidate_shared_soft_seconds=45.,
                candidate_score_slots_per_arm=8, candidate_forward_attempts_per_arm=16, policies=list(worker.POLICIES))
            saved_candidates = read(ROOT / "examples/cached/candidate_records.json")
            expected = read(ROOT / "examples/cached/policy_outcomes.json")["rows"]
            raw = []
            parent_dir = study / "parent_graph/shard-0"
            parent_dir.mkdir(parents=True)
            save(parent_dir / "config.json", worker.graph_config(config, VERSIONS)["compatible12_settings"])
            for i, entry in enumerate(saved_candidates["sources"]):
                original = entry["source"]
                s = dict(original, graph_verified=False, graph_result=None, graph_result_path=None)
                raw.append(s)
                record = copy.deepcopy(original["graph_result"])
                record.update(row=i, canonical_smiles=s["canonical_smiles"], seed=s["search_seed"])
                path = parent_dir / f"row-{i:03d}.result.json"
                save(path, record)
                path.with_name(path.name.replace(".result.json", ".graph.pkl")).write_bytes(pickle.dumps({"witness": s["canonical_smiles"]}))
            query_dir = study / "policy_graph/gpu0/queries"
            query_dir.mkdir(parents=True)
            save(query_dir.parent / "config.json", worker.graph_config(config, VERSIONS))
            qi = 0
            for outcome in expected:
                for attempt in outcome["attempts"]:
                    record = copy.deepcopy(attempt["query"])
                    record["row"] = qi
                    path = query_dir / f"query-{qi:05d}.result.json"
                    save(path, record)
                    path.with_name(path.name.replace(".result.json", ".graph.pkl")).write_bytes(pickle.dumps({"witness": record["request"]["smiles"]}))
                    qi += 1
            save(output / "method_config.json", config)
            save(output / "selected_inputs.json", {"rows": raw})
            save(output / "sources_scored.json", {"rows": raw, "supplied_inputs": raw,
                 "rf_model_dir": "models", "actual_new_source_predictions": 0})
            by_id = {entry["source"]["sample_id"]: entry for entry in saved_candidates["sources"]}
            def candidates(output, config):
                parents = read(output / "sources_with_parents.json")["rows"]
                save(output / "candidate_stage.json", {"sources": [
                    dict(by_id[s["sample_id"]], source=s) for s in parents],
                    "actual_new_candidate_predictions": 0, "actual_new_forward_attempts": 0})
            def no_search_backend():
                self.fail("Real/model backend constructed during exact software-fixture reuse")
            auditor = FakeAuditor()
            code = worker.run(output, backend_factory=no_search_backend, candidate_runner=candidates,
                              version_provider=lambda: VERSIONS, auditor_factory=lambda: auditor)
            self.assertEqual(code, 0)
            result = read(output / "method_result.json")
            fields = ("final_smiles", "final_score", "final_qed", "solved", "joint", "replaced", "status")
            indexed = {(r["sample_id"], r["policy"]): r for r in result["rows"]}
            for old in expected:
                new = indexed[(old["sample_id"], old["policy"])]
                self.assertEqual({k: new.get(k) for k in fields}, {k: old.get(k) for k in fields})
            self.assertEqual(result["actual_new_work"]["Graph_queries"], 0)
            self.assertEqual(result["actual_new_work"]["Graph_model_calls"], 0)
            self.assertEqual(result["exact_study_reused_queries"], 5)
            self.assertEqual(len(auditor.calls), 1)
            self.assertEqual(result["historical12_queries"], 2)
            # Counts stay old work; they do not enter physical new-work totals.
            self.assertGreater(result["old_graph_work"]["model_calls"], 0)


if __name__ == "__main__":
    unittest.main()
