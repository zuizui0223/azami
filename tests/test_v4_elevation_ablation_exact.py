"""Synthetic, source-free checks for matched V4 elevation ablation."""
import unittest

import numpy as np
import pandas as pd

from analysis.v4.elevation_ablation_exact import (
    M3, M4, wls_design_response, summary, check_public_payload,
)


class ElevationAblationTests(unittest.TestCase):
    def test_original_model_is_strict_superset(self):
        self.assertEqual(M3, ("chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean"))
        self.assertEqual(M4, (*M3, "topo_elevation"))

    def test_matched_nested_model_reduces_weighted_residuals(self):
        rng = np.random.default_rng(20261031)
        n = 360
        climate = rng.normal(size=(n, 3))
        elevation = rng.normal(size=n)
        y = 0.7 * climate[:, 0] - 0.2 * climate[:, 2] + 0.45 * elevation
        y += rng.normal(scale=0.8, size=n)
        weight = np.ones(n)
        group = np.repeat(np.arange(36), 10)
        m3 = wls_design_response(np.column_stack([np.ones(n), climate]), y, weight, group)
        m4 = wls_design_response(np.column_stack([np.ones(n), climate, elevation]), y, weight, group)
        self.assertLessEqual(np.sum(m4.resid**2), np.sum(m3.resid**2) + 1e-8)
        self.assertGreater(m3.aic, m4.aic)
        self.assertEqual(set(summary(m3, M3)), set(M3))

    def test_public_payload_rejects_private_columns_recursively(self):
        skeleton = {
            "S1": {"status": "EXACT_FROZEN_S1_REPRODUCED_NEW_EXPLORATORY_M3"},
            "S2": {"strata": [{"species": "Cirsium vulgare"}] * 6},
            "comparison": {"M3": list(M3), "M4": list(M4)},
        }
        check_public_payload(skeleton)
        skeleton["S2"]["strata"][0] = {"user_id": "hidden"}
        with self.assertRaisesRegex(ValueError, "PRIVATE_SOURCE_KEY"):
            check_public_payload(skeleton)

    def test_output_is_not_recast_as_original_holm_family(self):
        import inspect
        from analysis.v4 import elevation_ablation_exact as mod
        source = inspect.getsource(mod.run)
        self.assertIn("NEW_EXPLORATORY_SENSITIVITY_NOT_ORIGINAL_HOLM_INFERENCE", source)
        self.assertIn("original_20261024_frozen_estimates_unchanged", source)
        self.assertIn("no_elevation_exposure_coefficient_in_M3", source)


if __name__ == "__main__":
    unittest.main()
