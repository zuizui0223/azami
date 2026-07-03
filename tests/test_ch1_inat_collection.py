from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"


def load_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COLLECT = load_module("04_collect_inat_cirsium_metadata.py", "ch1_collect_inat")
QUEUE = load_module("05_build_image_screening_queue.py", "ch1_screen_queue")
DOWNLOAD = load_module("06_download_inat_screening_queue.py", "ch1_queue_download")
SCREEN = load_module("07_screen_downloaded_images_with_yolo.py", "ch1_yolo_screen")


class TestINaturalistCollectionFoundation(unittest.TestCase):
    def observation(self, obs_id: int, photo_ids: list[int], species: str, lat: float, lon: float, *, obscured: bool = False) -> dict:
        return {
            "id": obs_id,
            "taxon": {"id": obs_id + 1000, "name": species, "rank": "species"},
            "species_guess": species,
            "quality_grade": "research",
            "geojson": {"coordinates": [lon, lat]},
            "obscured": obscured,
            "annotations": [
                {
                    "controlled_attribute": {"id": 12, "label": "Flowers and Fruits"},
                    "controlled_value": {"id": 13, "label": "Flowers"},
                }
            ],
            "photos": [{"id": photo_id, "url": f"https://example.org/photos/{photo_id}/small.jpg"} for photo_id in photo_ids],
        }

    def test_photo_rows_preserve_provenance_and_flag_obscured_coordinates(self) -> None:
        rows = COLLECT.observation_to_photo_rows(self.observation(1, [11, 12], "Cirsium alpha", 35.0, 135.0))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["photo_key"], "obs_1_photo_11")
        self.assertTrue(rows[0]["inat_flowers_annotation"])
        self.assertTrue(rows[0]["coordinate_usable_for_environment"])
        self.assertEqual(rows[0]["source_file_stem"], "inat_photo_11")
        self.assertIn("/medium.", rows[0]["medium_image_url"])

        obscured = COLLECT.observation_to_photo_rows(self.observation(2, [21], "Cirsium alpha", 35.0, 135.0, obscured=True))[0]
        self.assertFalse(obscured["coordinate_usable_for_environment"])

    def test_queue_deduplicates_observations_and_retains_species_balance(self) -> None:
        observations = [
            self.observation(1, [11, 12], "Cirsium alpha", 35.0, 135.0),
            self.observation(2, [13], "Cirsium alpha", 36.0, 136.0),
            self.observation(3, [14], "Cirsium beta", 40.0, 140.0),
            self.observation(4, [15], "Cirsium beta", 41.0, 141.0),
        ]
        rows = [row for obs in observations for row in COLLECT.observation_to_photo_rows(obs)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "photo_metadata.csv"
            pd.DataFrame(rows).to_csv(source, index=False)
            out_dir = root / "queue"
            with patch.object(sys, "argv", [
                "queue",
                "--input", str(source),
                "--out-dir", str(out_dir),
                "--max-photos", "10",
                "--min-photos-per-species", "2",
                "--max-photos-per-species", "3",
                "--seed", "123",
            ]):
                QUEUE.main()
            queue = pd.read_csv(out_dir / "image_screening_queue.csv", dtype=str)
            self.assertEqual(len(queue), 4)
            self.assertTrue(queue["obs_id"].is_unique)
            self.assertSetEqual(set(queue["taxon_name"]), {"Cirsium alpha", "Cirsium beta"})
            self.assertFalse(queue["screen_download_filename"].str.contains("Cirsium", regex=False).any())

    def test_url_choice_and_crop_bounds_are_safe(self) -> None:
        row = pd.Series({"small_image_url": "small", "medium_image_url": "medium", "large_image_url": "large"})
        self.assertEqual(DOWNLOAD.selected_url(row, "medium"), "medium")
        self.assertEqual(SCREEN.safe_crop_box(10, 10, 20, 30, 100, 100, 0.1), (9, 8, 21, 32))
        with self.assertRaises(ValueError):
            SCREEN.safe_crop_box(10, 10, 10, 10, 100, 100, 0.1)


if __name__ == "__main__":
    unittest.main()
