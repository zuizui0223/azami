"""Source-free numerical and privacy tests for original-cohort BIO18 geography."""
from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from analysis.v4.s1_bio18_geographic_replication import (
    GEO_COLS, SECTORS, FOUR_WARM, geography, fit_taxon, summarise_geography,
)


class NativeBio18GeographyTests(unittest.TestCase):
    def synthetic_taxon_means(self, n=120):
        rng = np.random.default_rng(20261029)
        sector = np.arange(n) % 3
        x = rng.normal(size=(n, 4))
        y = 0.18*x[:, 0] + 0.06*x[:, 1] - 0.12*x[:, 2] + 0.03*x[:, 3] \
            + 0.30*(sector == 2) + rng.normal(scale=0.13, size=n)
        data = pd.DataFrame({"response": y,
                             **{key: x[:, i] for i, key in enumerate(FOUR_WARM)}})
        for j, key in enumerate(SECTORS):
            data[f"sector_{key}"] = (sector == j).astype(float)
        for j, key in enumerate(GEO_COLS):
            data[key] = rng.normal(size=n) if j < 7 else rng.normal(size=n)
        data.index = [f"private_taxon_{i}" for i in range(n)]
        return data

    def test_geography_is_fixed_and_response_blind(self):
        f, sec = geography([0, 0, 0, 45], [-31, -30, 60, 130])
        self.assertEqual(f.shape, (4, 8))
        self.assertEqual(sec.tolist(), [0, 1, 2, 2])
        self.assertAlmostEqual(f[0, 0], np.cos(np.deg2rad(-31)), places=10)
        with self.assertRaises(ValueError):
            geography([100], [0])
        with self.assertRaises(ValueError):
            geography([float("nan")], [0])

    def test_ols_fits_expected_sign_without_double_interpreting_p(self):
        frame = self.synthetic_taxon_means()
        unadjusted = fit_taxon(frame, "none")
        geography_adjusted = fit_taxon(frame, "longitude_sector")
        sphere = fit_taxon(frame, "spherical_degree2")
        for model in (unadjusted, geography_adjusted, sphere):
            self.assertEqual(model["status"], "estimated")
            self.assertGreater(model["standardised_beta"], 0.0)
            self.assertEqual(len(model["HC3_pointwise_95CI"]), 2)
            self.assertTrue(model["no_holm_correction_for_new_sensitivity_family"])
            # Neither observation identifiers nor individual taxon labels leak.
            self.assertNotIn("private_taxon", str(model))
        self.assertEqual(fit_taxon(frame.iloc[:12])["status"], "too_few_taxa")

    def test_sector_leave_one_out_and_invalid_design(self):
        frame = self.synthetic_taxon_means()
        output = summarise_geography({
            "photo_weighted": frame,
            "equal_5degree_cell": frame.copy(),
        }, verify_original_baseline=False)
        self.assertEqual(set(output), {"photo_weighted", "equal_5degree_cell"})
        for grouping in output.values():
            self.assertEqual(set(grouping["within_predeclared_longitude_sectors"]), set(SECTORS))
            self.assertEqual(set(grouping["leave_one_predeclared_sector_out"]), set(SECTORS))
            for sector in SECTORS:
                self.assertEqual(grouping["within_predeclared_longitude_sectors"][sector]["n_taxa"], 40)
                self.assertEqual(grouping["leave_one_predeclared_sector_out"][sector]["heldout_n_taxa"], 40)

    def test_invalid_geographic_design_not_silently_replaced(self):
        frame = self.synthetic_taxon_means(30)
        frame.loc[:, "sector_W180_W030"] = 0.0
        frame.loc[:, "sector_W030_E060"] = 0.0
        invalid = fit_taxon(frame, "longitude_sector")
        self.assertEqual(invalid["status"], "not_identifiable")
        self.assertEqual(fit_taxon(frame, "bad")["status"], "not_identifiable")


if __name__ == "__main__":
    unittest.main()
