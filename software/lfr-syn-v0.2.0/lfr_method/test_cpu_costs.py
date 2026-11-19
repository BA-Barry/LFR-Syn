"""Saved physical-cost accounting tests; no chemistry/model imports or calls."""
import copy
import unittest

from lfr_method.cpu_costs import actual_work_components, summarize_cpu_work


def formed_source():
    return {"source": {"sample_id": "source-a"}, "status": "candidate_scoring_complete",
            "costs": {"source_replay_forward_calls": 121, "retro_template_calls": 1071,
                      "source_RF_path_attributions": 1},
            "arms": {"task_contribution": {"costs": {"new_forward_attempts": 6}},
                     "permuted_contribution": {"costs": {"new_forward_attempts": 7}}}}


class CPUCostTests(unittest.TestCase):
    def test_full_native_count_includes_source_replay_and_one_shared_attribution(self):
        report = summarize_cpu_work([formed_source()])
        self.assertTrue(report["cost_observation_complete"])
        self.assertEqual(actual_work_components(report), {
            "source_replay_forward_reactions": 121, "candidate_forward_reactions": 13,
            "forward_reactions_total": 134, "retro_template_calls": 1071,
            "source_RF_path_attributions": 1})

    def test_saved_complete_record_retains_cumulative_counts_and_is_not_modified(self):
        source = formed_source()
        original = copy.deepcopy(source)
        first = summarize_cpu_work([source])
        resumed = summarize_cpu_work([copy.deepcopy(source)])
        self.assertEqual(first, resumed)
        self.assertEqual(source, original)
        self.assertEqual(resumed["forward_reactions_total"]["count"], 134)

    def test_legacy_formed_source_missing_root_fields_is_unknown_not_zero(self):
        source = formed_source()
        del source["costs"]
        report = summarize_cpu_work([source])
        self.assertFalse(report["cost_observation_complete"])
        self.assertIsNone(report["forward_reactions_total"]["count"])
        self.assertEqual(report["forward_reactions_total"]["observed_count"], 13)
        self.assertIsNone(report["components"]["retro_template_calls"]["count"])
        self.assertIsNone(report["components"]["source_RF_path_attributions"]["count"])
        self.assertTrue(report["forward_reactions_total"]["missing_counter_fields"])

    def test_missing_one_candidate_arm_does_not_make_a_complete_total(self):
        source = formed_source()
        del source["arms"]["permuted_contribution"]["costs"]
        report = summarize_cpu_work([source])
        self.assertIsNone(report["components"]["candidate_forward_reactions"]["count"])
        self.assertIsNone(report["forward_reactions_total"]["count"])
        self.assertEqual(report["forward_reactions_total"]["observed_count"], 127)

    def test_measured_zero_reused_root_and_forward_counters_are_preserved(self):
        source = formed_source()
        source["costs"]["source_replay_forward_calls"] = 0
        source["costs"]["retro_template_calls"] = 0
        for arm in source["arms"].values():
            arm["costs"]["new_forward_attempts"] = 0
            arm["costs"]["cached_forward_attempts"] = 16
        report = summarize_cpu_work([source])
        self.assertTrue(report["cost_observation_complete"])
        self.assertEqual(report["forward_reactions_total"]["count"], 0)
        self.assertEqual(report["components"]["source_RF_path_attributions"]["count"], 1)

    def test_terminal_retention_branch_has_known_no_formation_work(self):
        for status in ["retained_invalid", "retained_joint_success", "parent_execution_limitation"]:
            with self.subTest(status=status):
                report = summarize_cpu_work([{"source": {"sample_id": "retained"},
                                             "status": status, "costs": {}, "arms": {}}])
                self.assertTrue(report["cost_observation_complete"])
                self.assertEqual(report["forward_reactions_total"]["count"], 0)
                self.assertEqual(report["components"]["source_RF_path_attributions"]["count"], 0)

    def test_invalid_or_missing_counter_is_reported_without_fabricating_zero(self):
        for invalid in [-1, True, 1.0, None, "121"]:
            with self.subTest(invalid=invalid):
                source = formed_source()
                source["costs"]["source_replay_forward_calls"] = invalid
                report = summarize_cpu_work([source])
                self.assertFalse(report["forward_reactions_total"]["cost_observation_complete"])
                self.assertIsNone(report["forward_reactions_total"]["count"])
                self.assertEqual(report["forward_reactions_total"]["observed_count"], 13)

    def test_retained_label_with_formation_data_does_not_fill_missing_counters(self):
        record = {"source": {"sample_id": "legacy"}, "status": "retained_joint_success",
                  "arms": {"task_contribution": {"planned_candidates": [{"smiles": "CCO"}]}}}
        report = summarize_cpu_work([record])
        self.assertFalse(report["cost_observation_complete"])
        self.assertIsNone(report["forward_reactions_total"]["count"])


if __name__ == "__main__":
    unittest.main()
