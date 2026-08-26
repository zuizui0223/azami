from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

try:
    import cv2
except ModuleNotFoundError:  # Minimal CI environments may install submission extras later.
    cv2 = None

try:
    import statsmodels  # noqa: F401
except ModuleNotFoundError:
    statsmodels = None


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"
CONTRACT = V2 / "ontology" / "ch1_continuous_trait_contract.csv"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


VALIDATE = load("validate_continuous_contract", V2 / "87_validate_continuous_trait_contract.py")
BUILD = load("build_continuous_universe", V2 / "88_build_continuous_trait_universe.py")
MEASURE = load("measure_extended_continuous", V2 / "89_measure_extended_continuous_traits.py") if cv2 is not None else None
ANALYSE = load("analyse_continuous_universe", V2 / "90_run_continuous_trait_universe_climate.py") if statsmodels is not None else None
PREPARE = load("prepare_extended_cohort", V2 / "91_prepare_exhaustive_extended_cohort.py")


class ContinuousTraitContractTests(unittest.TestCase):
    def test_contract_is_category_free_and_formula_unique(self):
        contract = pd.read_csv(CONTRACT, dtype=str, keep_default_na=False)
        report = VALIDATE.validate_contract(contract)
        self.assertEqual(report["status"], "valid")
        self.assertEqual(report["n_endpoints"], len(contract))
        self.assertEqual(report["category_columns_allowed_in_analysis"], 0)
        self.assertEqual(contract["formula_id"].nunique(), len(contract))

    def test_duplicate_formula_is_rejected(self):
        contract = pd.read_csv(CONTRACT, dtype=str, keep_default_na=False)
        contract.loc[1, "formula_id"] = contract.loc[0, "formula_id"]
        with self.assertRaisesRegex(ValueError, "formula_id must be unique"):
            VALIDATE.validate_contract(contract)

    def test_universe_builder_uses_every_executed_continuous_endpoint(self):
        contract = pd.read_csv(CONTRACT, dtype=str, keep_default_na=False)
        primary_modules = {"orientation", "colour", "display", "shape"}
        primary_contract = contract[contract["module"].isin(primary_modules)]
        extended_contract = contract[~contract["module"].isin(primary_modules)]

        primary_observation = {"obs_id": ["o1"], "taxon_name": ["Cirsium test"]}
        primary_species = {"taxon_name": ["Cirsium test"]}
        extended_observation = {"obs_id": ["o1"], "taxon_name": ["Cirsium test"]}
        extended_species = {"taxon_name": ["Cirsium test"]}
        for index, row in primary_contract.reset_index(drop=True).iterrows():
            primary_observation[row.observation_variable] = [0.2 + index * 0.01]
            primary_species[row.species_variable] = [0.2 + index * 0.01]
            primary_observation[row.source_status_column] = [1]
            species_status = row.source_status_column.replace("n_usable_heads_", "n_observations_usable_")
            primary_species[species_status] = [1]
        for index, row in extended_contract.reset_index(drop=True).iterrows():
            extended_observation[row.observation_variable] = [0.3 + index * 0.01]
            extended_species[row.species_variable] = [0.3 + index * 0.01]
            extended_observation[row.source_status_column] = [1]
            extended_species[row.source_status_column] = [1]

        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            eo = directory / "extended_observation.csv"
            es = directory / "extended_species.csv"
            pd.DataFrame(extended_observation).to_csv(eo, index=False)
            pd.DataFrame(extended_species).to_csv(es, index=False)
            observation, species, report = BUILD.build_universe(
                contract,
                pd.DataFrame(primary_observation),
                pd.DataFrame(primary_species),
                str(eo),
                str(es),
            )
        self.assertEqual(report["status"], "ready")
        self.assertEqual(report["n_available_endpoints"], len(contract))
        self.assertEqual(len(observation), len(contract))
        self.assertEqual(len(species), len(contract))
        self.assertNotIn("trait_state", observation.columns)
        self.assertTrue(observation["measurement_available"].all())


class ExhaustiveCohortTests(unittest.TestCase):
    @staticmethod
    def metadata() -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "queue_id": "q1", "obs_id": "o1", "photo_id": "p1",
                    "taxon_name": "Cirsium alpha", "det_index": 0,
                    "yolo_conf": 0.91, "bbox_x1": 0, "bbox_y1": 0,
                    "bbox_x2": 190, "bbox_y2": 180,
                    "medium_image_url": "https://example.test/p1.jpg",
                },
                {
                    "queue_id": "q2", "obs_id": "o1", "photo_id": "p2",
                    "taxon_name": "Cirsium alpha", "det_index": 0,
                    "yolo_conf": 0.85, "bbox_x1": 0, "bbox_y1": 0,
                    "bbox_x2": 260, "bbox_y2": 230,
                    "medium_image_url": "https://example.test/p2.jpg",
                },
                {
                    "queue_id": "q3", "obs_id": "o2", "photo_id": "p3",
                    "taxon_name": "Cirsium beta", "det_index": 0,
                    "yolo_conf": 0.99, "bbox_x1": 0, "bbox_y1": 0,
                    "bbox_x2": 120, "bbox_y2": 220,
                    "medium_image_url": "https://example.test/p3.jpg",
                },
                {
                    "queue_id": "q4", "obs_id": "o3", "photo_id": "p4",
                    "taxon_name": "Cirsium gamma", "det_index": 0,
                    "yolo_conf": 0.98, "bbox_x1": 0, "bbox_y1": 0,
                    "bbox_x2": 300, "bbox_y2": 300,
                    "medium_image_url": "https://example.test/p4.jpg",
                },
            ]
        )

    def test_selection_is_trait_blind_and_deterministic(self):
        strict = pd.DataFrame({"obs_id": ["o1", "o2"]})
        selected, report = PREPARE.prepare_cohort(
            self.metadata(), strict, min_bbox_dimension=150, min_yolo_confidence=0.0
        )
        self.assertEqual(selected["obs_id"].tolist(), ["o1"])
        self.assertEqual(selected.iloc[0]["photo_id"], "p2")
        self.assertTrue(report["trait_blind"])
        self.assertEqual(report["n_selected_heads"], 1)
        self.assertEqual(report["n_selected_observations"], 1)

    def test_selection_rejects_trait_derived_fields(self):
        metadata = self.metadata()
        metadata["trait_state"] = "upright"
        with self.assertRaisesRegex(ValueError, "Trait-derived fields are forbidden"):
            PREPARE.prepare_cohort(
                metadata,
                pd.DataFrame({"obs_id": ["o1", "o2"]}),
                min_bbox_dimension=150,
                min_yolo_confidence=0.0,
            )


@unittest.skipIf(cv2 is None, "opencv submission dependency is not installed")
class ExtendedMeasurementTests(unittest.TestCase):
    @staticmethod
    def synthetic_head() -> np.ndarray:
        image = np.full((360, 360, 3), (45, 120, 45), dtype=np.uint8)
        cv2.ellipse(image, (180, 210), (78, 108), 0, 0, 360, (70, 105, 135), -1)
        cv2.ellipse(image, (180, 108), (65, 46), 0, 0, 360, (175, 45, 180), -1)
        for x in range(125, 236, 18):
            cv2.line(image, (x, 165), (x - 10, 255), (205, 205, 205), 2)
        for y in range(180, 280, 20):
            cv2.line(image, (130, y), (230, y), (65, 75, 95), 1)
        return image

    def test_extended_metrics_are_continuous_and_mirror_repeatable(self):
        args = argparse.Namespace(
            min_architecture_dimension=150,
            min_surface_dimension=300,
            min_architecture_sharpness=20.0,
            min_surface_sharpness=20.0,
        )
        image = self.synthetic_head()
        combined = MEASURE.combine(
            MEASURE.measure_once(image, args),
            MEASURE.measure_once(cv2.flip(image, 1), args),
        )
        self.assertEqual(combined["architecture_status"], "usable")
        for metric in MEASURE.ARCHITECTURE_METRICS:
            self.assertTrue(np.isfinite(combined[metric]), metric)
        self.assertGreater(combined["involucre_length_width_ratio"], 0)
        self.assertGreaterEqual(combined["involucre_spread_fraction"], 0)
        self.assertLessEqual(combined["involucre_spread_fraction"], 1)

    def test_flip_failure_excludes_only_the_affected_endpoint(self):
        args = argparse.Namespace(
            min_architecture_dimension=150,
            min_surface_dimension=300,
            min_architecture_sharpness=20.0,
            min_surface_sharpness=20.0,
        )
        original = MEASURE.measure_once(self.synthetic_head(), args)
        mirrored = MEASURE.measure_once(cv2.flip(self.synthetic_head(), 1), args)
        mirrored["involucre_projection_max"] = (
            original["involucre_projection_max"]
            + MEASURE.FLIP_TOLERANCE["involucre_projection_max"]
            + 0.01
        )
        combined = MEASURE.combine(original, mirrored)
        self.assertIn(
            "unstable_under_horizontal_flip",
            combined["involucre_projection_max_status"],
        )
        self.assertEqual(combined["involucre_length_width_ratio_status"], "usable")
        self.assertEqual(combined["architecture_status"], "one_or_more_endpoint_failures")


@unittest.skipIf(statsmodels is None, "statsmodels submission dependency is not installed")
class ContinuousAnalysisTests(unittest.TestCase):
    def test_linear_and_circular_traits_use_distinct_inferential_units(self):
        rng = np.random.default_rng(4)
        records = []
        environment_rows = []
        for taxon_index in range(12):
            taxon = f"Cirsium t{taxon_index:02d}"
            offset = rng.normal()
            for observation_index in range(10):
                obs_id = f"{taxon_index}-{observation_index}"
                bio1 = rng.normal()
                angle = (0.9 * bio1 + offset + rng.normal(scale=0.25))
                hue = 0.6 * bio1 + rng.normal(scale=0.35)
                environment_rows.append({"obs_id": obs_id, "bio1": bio1})
                for endpoint, value, circular in (
                    ("orientation_image_vertical_angle", angle, ""),
                    ("corolla_hue_sin", np.sin(hue), "corolla_hue"),
                    ("corolla_hue_cos", np.cos(hue), "corolla_hue"),
                ):
                    records.append({
                        "obs_id": obs_id,
                        "taxon_name": taxon,
                        "endpoint_id": endpoint,
                        "module": "orientation" if endpoint.startswith("orientation") else "colour",
                        "analysis_tier": "primary",
                        "analysis_eligible": True,
                        "circular_group": circular,
                        "value": value,
                    })
        coefficients, coverage, report = ANALYSE.run_models(
            pd.DataFrame(records),
            pd.DataFrame(environment_rows),
            ["bio1"],
            minimum_observations=100,
            minimum_taxa=10,
            permutations=99,
            seed=9,
        )
        self.assertEqual(report["n_linear_endpoints_modelled"], 1)
        self.assertEqual(report["n_circular_traits_modelled"], 1)
        self.assertEqual(set(coefficients["inferential_unit"]), {"linear_endpoint", "circular_joint"})
        orientation = coefficients[coefficients["endpoint_id"].eq("orientation_image_vertical_angle")].iloc[0]
        self.assertGreater(orientation["beta_std_within"], 0.5)
        self.assertIn("corolla_hue", set(coefficients["endpoint_id"]))
        self.assertEqual(len(coverage), 2)


if __name__ == "__main__":
    unittest.main()
