from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


VALIDATE = load_module("69_validate_submission_outputs.py", "ch1_submission_validate")
MANIFEST = load_module("70_build_submission_manifest.py", "ch1_submission_manifest")


class TestSubmissionContract(unittest.TestCase):
    def build_bundle(self, root: Path) -> tuple[Path, dict]:
        bundle = root / "bundle"
        bundle.mkdir()
        config = json.loads((V2 / "submission_config.json").read_text(encoding="utf-8"))

        endpoint_rows = []
        for index in range(9):
            endpoint_rows.append({
                "analysis_tier": "main",
                "trait_group": "primary",
                "observation_variable": f"obs_endpoint_{index}",
                "species_variable": f"species_endpoint_{index}",
            })
        for index in range(3):
            endpoint_rows.append({
                "analysis_tier": "auxiliary",
                "trait_group": "auxiliary",
                "observation_variable": f"obs_aux_{index}",
                "species_variable": f"species_aux_{index}",
            })
        registry = pd.DataFrame(endpoint_rows)
        registry.to_csv(bundle / "integrated_continuous_endpoint_registry.csv", index=False)

        observation = pd.DataFrame({
            "obs_id": [f"obs_{index:04d}" for index in range(3725)],
            "taxon_name": [f"Cirsium species_{index % 216:03d}" for index in range(3725)],
        })
        species = pd.DataFrame({
            "taxon_name": [f"Cirsium species_{index:03d}" for index in range(216)]
        })
        for index in range(9):
            observation[f"obs_endpoint_{index}"] = index + 0.1
            species[f"species_endpoint_{index}"] = index + 0.2
        for index in range(3):
            observation[f"obs_aux_{index}"] = index + 0.3
            species[f"species_aux_{index}"] = index + 0.4
        observation.to_csv(bundle / "integrated_primary_auxiliary_observation.csv", index=False)
        species.to_csv(bundle / "integrated_primary_auxiliary_species.csv", index=False)

        (bundle / "qc_bias_and_join_report.json").write_text(
            json.dumps({"n_heads": 6626}), encoding="utf-8"
        )

        pre_rows = []
        for endpoint in range(9):
            for predictor in range(4):
                pre_rows.append({
                    "analysis_tier": "main",
                    "model_scale": "observation_within_species",
                    "cohort": "positional_accuracy_le_10km",
                    "endpoint": f"obs_endpoint_{endpoint}",
                    "predictor": f"env_{predictor}",
                    "estimate_standardized": 0.0,
                    "p_fdr_within_tier_scale_cohort": 1.0,
                })
        pd.DataFrame(pre_rows).to_csv(
            bundle / "continuous_environment_model_coefficients.csv", index=False
        )

        (bundle / "historical_tree_audit_report.json").write_text(json.dumps({
            "n_input_taxa": 216,
            "gbotb_lcvp": {
                "n_direct_backbone_tips": 54,
                "n_randomized_scenario2_trees": 50,
            },
        }), encoding="utf-8")
        (bundle / "historical_constraint_model_report.json").write_text(json.dumps({
            "n_main_endpoints": 9,
            "n_auxiliary_endpoints": 3,
            "n_successful_model_fits": 522,
            "n_failed_model_fits": 0,
        }), encoding="utf-8")
        pd.DataFrame([{
            "covariance_model": "PGLS_Pagel_lambda",
            "pagel_lambda": 0.5,
        }]).to_csv(
            bundle / "historical_constraint_model_coefficients_clean.csv", index=False
        )
        pd.DataFrame([{
            "n_expected_random_trees": 50,
            "n_successful_random_trees": 50,
            "complete_random_tree_support": True,
        } for _ in range(36)]).to_csv(
            bundle / "historical_constraint_random_tree_summary.csv", index=False
        )
        pd.DataFrame({
            "accepted_taxon": species["taxon_name"],
            "ploidy_class": "unknown",
            "hybrid_status": "unknown",
        }).to_csv(bundle / "ploidy_hybrid_evidence_queue.csv", index=False)
        return bundle, config

    def test_validation_and_manifest_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, config = self.build_bundle(root)
            report = VALIDATE.validate_bundle(bundle, config)
            self.assertEqual(report["status"], "pass")
            self.assertGreater(report["n_checks"], 30)

            validation_path = root / "validation.json"
            validation_path.write_text(json.dumps(report), encoding="utf-8")
            output = root / "provenance"
            manifest = MANIFEST.build_manifest(
                bundle_root=bundle,
                validation_report=validation_path,
                config_path=V2 / "submission_config.json",
                out_dir=output,
                git_sha="fixture-sha",
            )
            self.assertEqual(manifest["analysis_version"], config["analysis_version"])
            self.assertEqual(manifest["n_files"], 10)
            self.assertTrue((output / "submission_file_manifest.csv").is_file())
            self.assertTrue((output / "submission_environment.json").is_file())

    def test_model_weight_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, config = self.build_bundle(root)
            (bundle / "accidental_best.pt").write_bytes(b"not a result table")
            with self.assertRaisesRegex(ValueError, "no_model_weights"):
                VALIDATE.validate_bundle(bundle, config)


if __name__ == "__main__":
    unittest.main()
