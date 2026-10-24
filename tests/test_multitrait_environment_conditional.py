"""Lightweight synthetic unit tests; no frozen-result claim without source archives."""
import unittest

import numpy as np
import pandas as pd

from analysis.v3 import run_multitrait_environment_conditional as mt
from analysis.v3 import run_construct_scale_upgrade as up


class MultitraitEnvironmentConditionalTests(unittest.TestCase):
    def test_weighted_residual_removes_linear_pc_effect(self):
        rng = np.random.default_rng(4)
        pcs = rng.normal(size=(150, 3))
        weights = np.linspace(1, 2, 150)
        y = 5 + 2 * pcs[:, 0] - 3 * pcs[:, 2]
        adjusted = mt.weighted_residual(y, pcs, weights)[:, 0]
        self.assertLess(np.max(np.abs(adjusted)), 1e-10)

    def test_weighted_pcs_are_orthogonal(self):
        rng = np.random.default_rng(100)
        x = rng.normal(size=(60, 9))
        weights = rng.uniform(0.5, 2, size=60)
        scores, explained = mt.weighted_pcs(x, weights, 3)
        weighted_cov = (scores * weights[:, None]).T @ scores / sum(weights)
        self.assertEqual(scores.shape, (60, 3))
        self.assertTrue(np.isfinite(scores).all())
        self.assertLess(np.max(np.abs(weighted_cov -
                                    np.diag(np.diag(weighted_cov)))), 1e-10)
        self.assertTrue(0 < sum(explained) <= 1.0)

    def test_conditioned_cohort_and_pair_dimensions(self):
        rng = np.random.default_rng(812)
        frame_rows = []
        taxa = [f"taxon_{i:02d}" for i in range(12)]
        columns = [f"{construct}__0" for construct in up.CORE]
        fmap = {name: [f"{name}__0"] for name in up.CORE}
        for i, taxon in enumerate(taxa):
            for j in range(6):
                env = rng.normal(size=9) + 0.3 * i
                row = {"taxon_name": taxon, "obs_id": f"{i}_{j}"}
                row.update(dict(zip(mt.axis.PREDICTORS, env)))
                row.update({col: 0.4 * env[0] + 0.2 * i + rng.normal()
                            for col in columns})
                frame_rows.append(row)
        frame = pd.DataFrame(frame_rows)
        frame.attrs["feature_map"] = fmap
        within, among, explained = mt.conditioned_data(frame, mt.axis.PREDICTORS, 3)
        result = mt.conditioned_pairs(within, among, fmap)
        self.assertEqual(len(result), 36)
        self.assertEqual(len(within), 72)
        self.assertEqual(len(among), 12)
        self.assertTrue(np.isfinite(result.conditioned_within_rv).all())
        self.assertTrue(np.isfinite(result.conditioned_among_rv).all())
        self.assertEqual(len(explained["among_taxon_pc_explained"]), 3)
        matched, intervals = mt.matched_n(within, among, fmap, result, 50, 918)
        self.assertEqual(len(matched), 50)
        self.assertEqual(len(intervals), 36)
        self.assertTrue((matched.pairs_stronger_among <= 36).all())

    def test_rank_fail_closed(self):
        x = np.ones((12, 9))
        with self.assertRaises(ValueError):
            mt.weighted_pcs(x, np.ones(12), 3)
        with self.assertRaises(ValueError):
            mt.matched_n(pd.DataFrame(), pd.DataFrame(), {}, pd.DataFrame(), 10, 1)


if __name__ == "__main__":
    unittest.main()
