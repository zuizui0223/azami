from __future__ import annotations

import argparse
import importlib.util
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


QUEUE = load("71_build_within_species_expansion_queue.py", "exhaustive_queue")
COHORT = load("74_build_postprediction_analysis_cohorts.py", "postprediction_cohort")


class TestExhaustivePhotoPipeline(unittest.TestCase):
    def photo_fixture(self) -> pd.DataFrame:
        rows = []
        for species, observations in (("Cirsium alpha", 3), ("Cirsium beta", 2)):
            for obs in range(observations):
                for photo_index in range(2):
                    rows.append({
                        "obs_id": f"{species}_{obs}",
                        "photo_id": f"p_{species}_{obs}_{photo_index}",
                        "photo_index": str(photo_index),
                        "taxon_name": species,
                        "taxon_rank": "species",
                        "quality_grade": "research",
                        "captive": "false",
                        "latitude": str(30 + obs),
                        "longitude": str(130 + obs),
                        "coordinate_usable_for_environment": "true",
                        "positional_accuracy": "1000",
                        "geoprivacy": "open",
                        "obscured": "false",
                        "observed_on": "2026-01-01",
                        "small_image_url": "s",
                        "medium_image_url": "m",
                        "large_image_url": "l",
                        "source_file_stem": f"stem_{species}_{obs}_{photo_index}",
                        "photo_license_code": "",
                        "user_id": f"u_{obs}",
                        "user_login": f"observer_{obs}",
                        "photo_attribution": "fixture",
                    })
        rows.append({
            **rows[0],
            "photo_id": "captive_photo",
            "source_file_stem": "captive_photo",
            "captive": "true",
        })
        return pd.DataFrame(rows)

    def queue_args(self) -> argparse.Namespace:
        return argparse.Namespace(
            quality_grades="research",
            photo_licenses="any",
            require_open_coordinate=False,
            max_positional_accuracy_m=0,
            target_species="",
        )

    def test_all_non_captive_photos_are_retained_before_prediction(self) -> None:
        queue, audit = QUEUE.build_queue(self.photo_fixture(), self.queue_args())
        self.assertEqual(len(queue), 10)
        self.assertEqual(queue["obs_id"].nunique(), 5)
        self.assertEqual(queue["photo_id"].nunique(), 10)
        self.assertEqual(set(queue["sampling_role"]), {"exhaustive_pre_prediction"})
        self.assertEqual(int(audit["n_photos"].sum()), 10)

    def test_multiple_photos_from_one_observation_are_not_collapsed(self) -> None:
        queue, _ = QUEUE.build_queue(self.photo_fixture(), self.queue_args())
        counts = queue.groupby("obs_id")["photo_id"].nunique()
        self.assertTrue(counts.eq(2).all())

    def test_trait_columns_never_control_queue_selection(self) -> None:
        data = self.photo_fixture()
        data["orientation_angle_degrees"] = range(len(data))
        queue, _ = QUEUE.build_queue(data, self.queue_args())
        self.assertEqual(len(queue), 10)
        self.assertNotIn("orientation_angle_degrees", queue.columns)

    def observation_fixture(self) -> pd.DataFrame:
        rows = []
        for species in ("Cirsium alpha", "Cirsium beta"):
            for index in range(55):
                rows.append({
                    "obs_id": f"{species}_{index}",
                    "taxon_name": species,
                    "latitude": 30 + index * 0.1,
                    "longitude": 130 + index * 0.1,
                    "positional_accuracy": 1000 if index < 50 else 20000,
                    "coordinate_usable_for_environment": True,
                    "orientation_angle_degrees_median": float(index),
                })
        return pd.DataFrame(rows)

    def cohort_args(self) -> argparse.Namespace:
        return argparse.Namespace(
            strict_positional_accuracy_m=10_000,
            cell_deg=0.25,
            between_species_max_per_taxon=40,
            seed="fixture",
        )

    def test_thinning_happens_only_in_derived_cohorts(self) -> None:
        source = self.observation_fixture()
        cohorts = COHORT.build_cohorts(source, self.cohort_args())
        self.assertEqual(len(cohorts["all_detected"]), len(source))
        self.assertEqual(len(cohorts["strict_10km"]), 100)
        self.assertLess(len(cohorts["strict_spatial_thin"]), 100)
        self.assertLessEqual(
            int(cohorts["between_species_balanced"].groupby("taxon_name").size().max()),
            40,
        )

    def test_postprediction_selection_is_trait_blind(self) -> None:
        source = self.observation_fixture()
        first = COHORT.build_cohorts(source, self.cohort_args())
        source["orientation_angle_degrees_median"] = source["orientation_angle_degrees_median"] * -999
        second = COHORT.build_cohorts(source, self.cohort_args())
        self.assertEqual(
            list(first["strict_spatial_thin"]["obs_id"]),
            list(second["strict_spatial_thin"]["obs_id"]),
        )


if __name__ == "__main__":
    unittest.main()
