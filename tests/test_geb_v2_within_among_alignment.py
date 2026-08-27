from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis" / "run_geb_v2_within_among_alignment.py"


def load_module():
    spec = importlib.util.spec_from_file_location("geb_v2_within_among", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ANALYSE = load_module()


class WithinAmongAlignmentTests(unittest.TestCase):
    def test_circular_hue_is_one_inferential_unit(self):
        traits = pd.DataFrame([
            {"endpoint_id": "corolla_hue_sin", "module": "colour", "analysis_tier": "primary", "analysis_eligible_bool": True, "circular_group": "corolla_hue", "value": 0.5},
            {"endpoint_id": "corolla_hue_cos", "module": "colour", "analysis_tier": "primary", "analysis_eligible_bool": True, "circular_group": "corolla_hue", "value": 0.8},
            {"endpoint_id": "linear_test", "module": "shape", "analysis_tier": "candidate", "analysis_eligible_bool": True, "circular_group": "", "value": 1.0},
        ])
        units = ANALYSE.inferential_units(traits)
        self.assertEqual(len(units), 2)
        hue = next(unit for unit in units if unit["endpoint_id"] == "corolla_hue")
        self.assertEqual(hue["inferential_unit"], "circular_joint")
        self.assertEqual(set(hue["members"]), {"corolla_hue_sin", "corolla_hue_cos"})

    def test_primary_among_scope_detects_strong_taxon_level_signal(self):
        rng = np.random.default_rng(7)
        traits = []
        environment = []
        for taxon_index in range(30):
            taxon = f"Cirsium test{taxon_index:02d}"
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
                    "value": 2.0 * taxon_index + rng.normal(0, 0.2),
                })
        result, taxon_values = ANALYSE.run_scope(
            pd.DataFrame(traits),
            pd.DataFrame(environment),
            ANALYSE.DEFAULT_PREDICTORS,
            min_trait_observations=5,
            minimum_taxa=20,
            permutations=999,
            seed=123,
            scope="among_taxon_min5",
        )
        bio1 = result[result["predictor"].eq("chelsa_bio01")].iloc[0]
        self.assertEqual(bio1["status"], "ok")
        self.assertGreater(bio1["beta_std_among"], 0.95)
        self.assertLessEqual(bio1["p_value"], 0.002)
        self.assertTrue(bool(bio1["fdr_significant_0_05"]))
        self.assertEqual(taxon_values["taxon_name"].nunique(), 30)

    def test_alignment_classes_and_linear_sign_concordance(self):
        among = pd.DataFrame([
            {
                "scope": "among_taxon_min5",
                "endpoint_id": "linear_test",
                "module": "shape",
                "analysis_tier": "candidate",
                "inferential_unit": "linear_endpoint",
                "predictor": "chelsa_bio01",
                "status": "ok",
                "beta_std_among": 0.4,
                "p_value": 0.001,
                "q_fdr_bh_within_tier_scope": 0.004,
                "fdr_significant_0_05": True,
                "n_taxa": 30,
            }
        ])
        within = pd.DataFrame([
            {
                "endpoint_id": "linear_test",
                "module": "shape",
                "analysis_tier": "candidate",
                "inferential_unit": "linear_endpoint",
                "predictor": "chelsa_bio01",
                "status": "ok",
                "beta_std_within": 0.2,
                "p_value": 0.002,
                "q_fdr_bh_within_tier": 0.01,
                "fdr_significant_0_05": True,
            }
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "within.csv"
            within.to_csv(path, index=False)
            aligned = ANALYSE.build_alignment(among, path, "among_taxon_min5")
        row = aligned.iloc[0]
        self.assertEqual(row["cross_scale_class"], "both_scales")
        self.assertTrue(bool(row["linear_sign_concordant"]))


if __name__ == "__main__":
    unittest.main()
