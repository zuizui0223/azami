from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
sys.path.insert(0, str(ANALYSIS))


def load_module():
    path = ANALYSIS / "run_geb_v2_expanded_niche_sorting.py"
    spec = importlib.util.spec_from_file_location("geb_v2_expanded_niche", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NICHE = load_module()


class ExpandedNicheSortingTests(unittest.TestCase):
    def test_batched_null_metrics_are_finite_and_keep_group_sizes(self):
        env = np.array([
            [-2.0, 0.0, 1.0],
            [-1.0, 1.0, 0.0],
            [0.0, -1.0, 1.0],
            [1.0, 0.0, -1.0],
            [2.0, 1.0, 0.0],
            [3.0, -1.0, 1.0],
        ])
        labels = np.array([-1, -1, 0, 0, 1, 1], dtype=np.int8)
        centroid, overlap = NICHE.batched_null_metrics(env, labels, 50, np.random.default_rng(1))
        self.assertEqual(centroid.shape, (50,))
        self.assertEqual(overlap.shape, (50,))
        self.assertTrue(np.isfinite(centroid).all())
        self.assertTrue(np.isfinite(overlap).all())
        self.assertTrue(((overlap >= 0) & (overlap <= 1)).all())

    def test_common_taxon_scope_detects_strong_linear_environmental_sorting(self):
        rng = np.random.default_rng(9)
        traits = []
        environment = []
        for taxon_index in range(32):
            taxon = f"Cirsium niche{taxon_index:02d}"
            for replicate in range(5):
                obs_id = f"{taxon_index}_{replicate}"
                environment.append({
                    "obs_id": obs_id,
                    "taxon_name": taxon,
                    "chelsa_bio01": taxon_index + rng.normal(0, 0.05),
                    "chelsa_bio04": rng.normal(),
                    "chelsa_bio12": rng.normal(),
                    "chelsa_bio15": rng.normal(),
                })
                traits.append({
                    "obs_id": obs_id,
                    "taxon_name": taxon,
                    "endpoint_id": "linear_test",
                    "module": "shape",
                    "analysis_tier": "candidate",
                    "analysis_eligible_bool": True,
                    "circular_group": "",
                    "value": taxon_index + rng.normal(0, 0.1),
                })
        linear, hue, scores, report = NICHE.run_scope(
            pd.DataFrame(traits),
            pd.DataFrame(environment),
            NICHE.ALIGN.DEFAULT_PREDICTORS,
            min_trait_observations=5,
            minimum_common_taxa=20,
            permutations=999,
            seed=44,
            scope="test_scope",
        )
        self.assertEqual(report["n_common_taxa"], 32)
        self.assertEqual(len(hue), 0)
        self.assertEqual(scores["taxon_name"].nunique(), 32)
        row = linear.iloc[0]
        self.assertLess(row["centroid_permutation_p"], 0.01)
        self.assertLess(row["overlap_permutation_p"], 0.05)
        self.assertTrue(bool(row["both_metrics_fdr_0_05"]))

    def test_joint_hue_uses_no_quartile_rows(self):
        x = np.linspace(-2, 2, 40)
        y = np.column_stack([x, x * 0.5 + 0.1])
        env = np.column_stack([x, np.sin(x), np.cos(x)])
        observed = NICHE.multivariate_r2(y, env)
        self.assertGreater(observed, 0.9)


if __name__ == "__main__":
    unittest.main()
