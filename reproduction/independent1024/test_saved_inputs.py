"""Ordinary malformed-input checks; no model or molecular calculations."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import reproduce_independent1024 as reproduce


class SavedInputChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = ROOT / 'data'
        cls.design = reproduce.read(data / 'design.json')
        cls.sources = reproduce.read(data / 'sources_with_parents.json')
        cls.rows = reproduce.read(data / 'policy_rows.json')
        cls.predictions = reproduce.read(data / 'predictions.json')

    def malformed_membership(self, change, pattern):
        with tempfile.TemporaryDirectory(prefix='saved input 空格 ') as tmp:
            data = Path(tmp)
            sources, rows = copy.deepcopy(self.sources), copy.deepcopy(self.rows)
            change(sources, rows)
            reproduce.save(data / 'design.json', self.design)
            reproduce.save(data / 'sources_with_parents.json', sources)
            reproduce.save(data / 'policy_rows.json', rows)
            with self.assertRaisesRegex(ValueError, pattern):
                reproduce.load_inputs(data)

    def test_missing_raw_position_is_not_removed_from_denominator(self):
        self.malformed_membership(lambda sources, rows: sources['rows'].pop(), 'denominator')

    def test_duplicate_policy_id_cannot_replace_a_missing_position(self):
        self.malformed_membership(lambda sources, rows: rows['rows'].__setitem__(1, copy.deepcopy(rows['rows'][0])), 'policy membership')

    def test_duplicate_raw_position_with_distinct_id_is_rejected(self):
        def change(sources, rows):
            sources['rows'][1]['position'] = sources['rows'][0]['position']
            sources['rows'][1]['requested_search_seed'] = sources['rows'][0]['requested_search_seed']
        self.malformed_membership(change, 'raw draw position')

    def test_altered_fixed_candidate_RF_cannot_retain_joint_success(self):
        def change(sources, rows):
            row = next(r for r in rows['rows'] if r['policy'] == 'B_fixedtask6plus6' and r['replaced'])
            row['final_score'] = .0
        self.malformed_membership(change, 'Guiding RF endpoint|qualification violated')

    def test_incomplete_checkpoint_vector_is_rejected(self):
        predictions = copy.deepcopy(self.predictions)
        next(iter(predictions['predictions'].values()))['folds'].pop()
        with self.assertRaisesRegex(ValueError, 'ten finite two-task'):
            reproduce.vectors_analysis.vector_records(predictions)

    def test_swapped_task_columns_are_rejected(self):
        predictions = copy.deepcopy(self.predictions)
        predictions['tasks'].reverse()
        with self.assertRaisesRegex(ValueError, 'Wrong or reordered'):
            reproduce.vectors_analysis.vector_records(predictions)

    def test_fabricated_saved_mean_is_rejected(self):
        predictions = copy.deepcopy(self.predictions)
        record = next(iter(predictions['predictions'].values()))
        record['mean'][0] = .0 if record['mean'][0] > .5 else 1.
        with self.assertRaisesRegex(ValueError, 'means conflict'):
            reproduce.vectors_analysis.vector_records(predictions)

    def test_changed_mechanism_source_is_not_a_new_analysis_cohort(self):
        with tempfile.TemporaryDirectory(prefix='saved mechanism 空格 ') as tmp:
            data = Path(tmp)
            stage = reproduce.read(ROOT / 'data/candidate_stage.json')
            stage['sources'][0]['source']['task_score'] = .0
            reproduce.save(data / 'candidate_stage.json', stage)
            reproduce.save(data / 'sources_with_parents.json', self.sources)
            with self.assertRaisesRegex(ValueError, 'saved source identities'):
                reproduce.study_analysis.mechanism_counts(data, self.design)


if __name__ == '__main__':
    unittest.main()
