"""Sampling figures must retain cohort counts without releasing point records."""
import importlib
import json
import unittest

import pandas as pd


class CoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.api = importlib.import_module('analysis.v4.submission_sampling_coverage')

    def test_sparse_cells_are_removed_but_their_total_is_accounted_for(self):
        frame = pd.DataFrame({
            'latitude': [1.] * 5 + [11.] * 4,
            'longitude': [1.] * 5 + [11.] * 4,
            'accepted_key': ['a', 'a', 'a', 'b', 'b'] + ['private-taxon'] * 4,
            'obs_id': ['private-record-' + str(i) for i in range(9)],
            'user_id': ['private-user'] * 9,
        })
        panel = self.api.aggregate_panel(frame, 'test', taxon_column='accepted_key')
        self.assertEqual(panel['total_records'], 9)
        self.assertEqual(panel['shown_records'], 5)
        self.assertEqual(panel['suppressed_records'], 4)
        self.assertEqual(panel['total_taxa'], 3)
        self.assertEqual(panel['cells'], [{'lat_bin': 18, 'lon_bin': 36, 'records': 5, 'taxa': 2}])
        self.assertNotIn('private-', json.dumps(panel))
        self.api.validate_panel(panel)

    def test_world_edges_are_kept_inside_the_grid(self):
        frame = pd.DataFrame({
            'latitude': [-90.] * 5 + [90.] * 5,
            'longitude': [-180.] * 5 + [180.] * 5,
            'accepted_key': ['a'] * 10,
        })
        panel = self.api.aggregate_panel(frame, 'test', taxon_column='accepted_key')
        self.assertEqual([(r['lat_bin'], r['lon_bin']) for r in panel['cells']], [(0, 0), (35, 71)])
        self.assertEqual(panel['shown_records'], 10)

    def test_invalid_source_coordinates_fail_instead_of_silently_disappearing(self):
        for latitude, longitude in [(91., 0.), (0., -181.), (float('nan'), 0.)]:
            frame = pd.DataFrame({'latitude': [latitude] * 5,
                                  'longitude': [longitude] * 5,
                                  'accepted_key': ['a'] * 5})
            with self.subTest(latitude=latitude, longitude=longitude):
                with self.assertRaisesRegex(ValueError, 'COORDINATE'):
                    self.api.aggregate_panel(frame, 'test', taxon_column='accepted_key')

    def test_public_schema_rejects_identifiers_and_count_drift(self):
        frame = pd.DataFrame({'latitude': [1.] * 5, 'longitude': [2.] * 5,
                              'accepted_key': ['a'] * 5})
        original = self.api.aggregate_panel(frame, 'test', taxon_column='accepted_key')
        for mutation in ['identifier', 'extra_coordinate', 'low_count', 'wrong_total']:
            panel = json.loads(json.dumps(original))
            if mutation == 'identifier':
                panel['obs_id'] = '123'
            elif mutation == 'extra_coordinate':
                panel['cells'][0]['latitude'] = 1.
            elif mutation == 'low_count':
                panel['cells'][0]['records'] = 4
            else:
                panel['shown_records'] += 1
            with self.subTest(mutation=mutation):
                with self.assertRaises(ValueError):
                    self.api.validate_panel(panel)

    def test_panel_order_and_cohort_counts_are_frozen(self):
        self.assertEqual(self.api.EXPECTED_PANEL_COUNTS, {
            's1_angle': 168907,
            'arvense_native': 22936,
            'arvense_introduced': 26171,
            'vulgare_native': 17335,
            'vulgare_introduced': 27200,
        })
        self.assertEqual(self.api.CELL_DEGREES, 5)
        self.assertEqual(self.api.MIN_CELL_RECORDS, 5)


if __name__ == '__main__':
    unittest.main()
