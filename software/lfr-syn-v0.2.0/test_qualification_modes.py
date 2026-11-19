"""Qualification integration tests, using saved numbers and injected queries only."""
import copy
import unittest

from lfr_method import graph_core, fixedtask_kernel
from lfr_method.qualification import settings, evaluator, report_qualification
from lfr_method import test_graph_worker as worker_test


def source(task="gsk3b"):
    return dict(sample_id="qualification-example", cohort_id=task + "-rep0", task=task,
                replicate=0, position=0, task_index=int(task == "jnk3"),
                canonical_smiles="parent", valid=True, task_score=.9, qed=.6,
                search_seed=123, graph_verified=False,
                graph_result={"status": "unsolved_within_budget"}, graph_result_path=None)


def candidate(smiles="first", score=.2, qed=.6):
    return dict(smiles=smiles, task_score=score, qed=qed, task_pass=False)


class QualificationModes(unittest.TestCase):
    def test_fixed_accepts_above_task_below_source_without_mutating_old_flag(self):
        original = {"scored_candidates": [candidate()]}
        retained = copy.deepcopy(original)
        self.assertEqual(fixedtask_kernel.candidate_pool(source(), original)[0]["smiles"], "first")
        self.assertEqual(graph_core.candidate_pool(source(), original), [])
        self.assertEqual(original, retained)

    def test_both_fixed_task_floors_are_inclusive(self):
        for task, floor in graph_core.THRESHOLDS.items():
            with self.subTest(task=task):
                pool = fixedtask_kernel.candidate_pool(source(task), {"scored_candidates": [
                    candidate("at_floor", floor), candidate("below_floor", floor - .001)]})
                self.assertEqual([c["smiles"] for c in pool], ["at_floor"])

    def test_qed_tolerance_and_source_exclusion_are_unchanged(self):
        pool = fixedtask_kernel.candidate_pool(source(), {"scored_candidates": [
            candidate("boundary", qed=.55), candidate("too_low", qed=.549), candidate("parent")]})
        self.assertEqual([c["smiles"] for c in pool], ["boundary"])

    def test_duplicate_and_score_slot_limit_remain_errors(self):
        for candidates in ([candidate(), candidate()], [candidate(str(i)) for i in range(9)]):
            with self.assertRaises(ValueError):
                fixedtask_kernel.candidate_pool(source(), {"scored_candidates": candidates})

    def test_fixed_route_execution_keeps_six_plus_six_remaining_time(self):
        requests = []
        def query(smiles, seed, calls, seconds):
            requests.append((smiles, seed, calls, seconds))
            solved = smiles == "second"
            return dict(status="solved" if solved else "unsolved_within_budget",
                        verified_solved=solved, model_calls=6 if not solved else 2,
                        search_seconds=12. if not solved else 3., actual_new_calls=0)
        row = evaluator("fixed_task")(source(), "B_fixedtask6plus6",
                {"scored_candidates": [candidate(), candidate("second", .19)]}, query)
        report_qualification(row, "fixed_task")
        self.assertEqual(requests, [("first", 123, 6, 45.), ("second", 123, 6, 33.)])
        self.assertTrue(row["replaced"])
        self.assertAlmostEqual(row["guiding_RF_change"], -.71)
        self.assertFalse(row["guiding_RF_nondegradation_guarantee"])
        self.assertEqual(row["logical_model_calls"], 8)

    def test_strict_guard_remains_and_unknown_rf_change_is_not_zero(self):
        strict = source()
        arm = {"scored_candidates": [dict(candidate(score=.95), task_pass=True)]}
        row = evaluator("source_preserving")(strict, "B_joint6plus6", arm,
                lambda *args: dict(status="error", model_calls=0, search_seconds=0.))
        report_qualification(row, "source_preserving")
        self.assertTrue(row["guiding_RF_nondegradation_guarantee"])
        self.assertIsNone(row["guiding_RF_change"])

    def test_historical_config_defaults_to_strict_and_mismatched_labels_rejected(self):
        self.assertEqual(settings({"policies": list(graph_core.POLICIES)})[0], "source_preserving")
        with self.assertRaisesRegex(ValueError, "conflict"):
            settings({"qualification": "fixed_task", "policies": list(graph_core.POLICIES)})
        with self.assertRaisesRegex(ValueError, "Unknown"):
            settings({"qualification": "other"})


class FixedWorkerIntegration(unittest.TestCase):
    # Use the existing ordinary worker harness with changed qualification, not
    # another copy of its inherited tests.
    setUp = worker_test.WorkerTests.setUp
    setup_rows = worker_test.WorkerTests.setup_rows
    run_backend = worker_test.WorkerTests.run_backend
    def test_fixed_worker_dispatches_the_validated_kernel(self):
        self.config.update(qualification="fixed_task", policies=list(fixedtask_kernel.POLICIES))
        raw = source()
        self.setup_rows([raw])
        from lfr_method.graph_core import save, read
        save(self.output / "method_config.json", self.config)
        def candidates(output, config):
            parents = read(output / "sources_with_parents.json")["rows"]
            save(output / "candidate_stage.json", {"sources": [
                {"source": s, "arms": {name: {"scored_candidates": [candidate()]}
                 for name in ("task_contribution", "permuted_contribution")}} for s in parents],
                 "actual_new_candidate_predictions": 0, "actual_new_forward_attempts": 0})
        backend = worker_test.FakeBackend([("unsolved_within_budget", 12, 1.), ("solved", 1, 1.)])
        self.assertEqual(self.run_backend(backend, candidates), 0)
        result = read(self.output / "method_result.json")
        self.assertEqual(result["qualification"], "fixed_task")
        self.assertEqual({r["policy"] for r in result["rows"]}, set(fixedtask_kernel.POLICIES))
        self.assertTrue(all(r["replaced"] and r["guiding_RF_change"] == -.7 for r in result["rows"]))
        self.assertTrue(all(not r["guiding_RF_nondegradation_guarantee"] for r in result["rows"]))


if __name__ == "__main__":
    unittest.main()
