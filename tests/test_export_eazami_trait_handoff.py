import unittest

import pandas as pd

from analysis.export_eazami_trait_handoff import (
    AUXILIARY_COLUMNS,
    AUXILIARY_PROVENANCE_ONLY_ENDPOINTS,
    FINAL_AUXILIARY_INFERENTIAL_ENDPOINTS,
    PRIMARY_COLUMNS,
    build_handoff,
)


class TestEAzamiTraitHandoff(unittest.TestCase):
    def fixture(self):
        payload = {}
        for column in PRIMARY_COLUMNS + AUXILIARY_COLUMNS:
            if column == "taxon_name":
                payload[column] = ["Cirsium beta", "Cirsium alpha"]
            else:
                payload[column] = [1.0, 2.0]
        return pd.DataFrame(payload)

    def test_build_handoff_sorts_and_preserves_boundaries(self):
        result = build_handoff(self.fixture())
        self.assertEqual(result["taxon_name"].tolist(), ["Cirsium alpha", "Cirsium beta"])
        self.assertTrue(result["phylogeny_status"].str.contains("unassigned").all())
        self.assertTrue(result["involucre_spine_proxy_status"].str.contains("proxy").all())
        self.assertTrue(result["distribution_scope"].str.contains("not_full").all())
        self.assertTrue(result["distribution_requirement"].str.contains("observation_level").all())
        self.assertTrue(result["auxiliary_inference_scope"].str.contains("peak_count_provenance_only").all())
        self.assertIn("orientation_angle_degrees_species_median", result.columns)
        self.assertIn("spine_relative_length_max_proxy_species_median", result.columns)

    def test_auxiliary_scope_matches_final_manuscript(self):
        self.assertEqual(
            FINAL_AUXILIARY_INFERENTIAL_ENDPOINTS,
            [
                "involucre_projection_roughness",
                "involucre_spread_fraction",
                "spine_relative_length_max_proxy",
            ],
        )
        self.assertEqual(AUXILIARY_PROVENANCE_ONLY_ENDPOINTS, ["spine_peak_count_proxy"])

    def test_duplicate_taxa_fail_closed(self):
        frame = self.fixture()
        frame.loc[1, "taxon_name"] = frame.loc[0, "taxon_name"]
        with self.assertRaises(ValueError):
            build_handoff(frame)


if __name__ == "__main__":
    unittest.main()
