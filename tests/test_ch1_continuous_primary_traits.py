from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import cv2
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MEASURE = load_module(
    "52_measure_primary_traits_continuous.py", "continuous_primary_measurement"
)
AGGREGATE = load_module(
    "53_aggregate_primary_trait_continuous.py", "continuous_primary_aggregation"
)
REBUILD = load_module(
    "54_rebuild_global_continuous_trait_packet.py", "continuous_primary_rebuild"
)


class TestContinuousPrimaryTraits(unittest.TestCase):
    def test_colour_metrics_are_flip_invariant(self) -> None:
        image = np.full((128, 128, 3), (40, 120, 40), dtype=np.uint8)
        cv2.ellipse(
            image, (64, 60), (35, 28), 0, 0, 360, (180, 60, 220), -1
        )
        left = MEASURE.colour_measurement(image)
        right = MEASURE.colour_measurement(cv2.flip(image, 1))
        self.assertGreater(left["floral_pixel_fraction"], 0.5)
        self.assertAlmostEqual(
            left["median_lab_chroma"], right["median_lab_chroma"], places=5
        )
        self.assertLess(
            MEASURE.circular_difference_degrees(
                left["mean_hue_degrees"], right["mean_hue_degrees"]
            ),
            1.0,
        )

    def test_shape_metrics_distinguish_circle_and_elongated_outline(self) -> None:
        circle = np.full((160, 160, 3), (40, 120, 40), dtype=np.uint8)
        cv2.ellipse(
            circle, (80, 80), (42, 42), 0, 0, 360, (80, 70, 120), -1
        )
        elongated = np.full((180, 140, 3), (40, 120, 40), dtype=np.uint8)
        cv2.ellipse(
            elongated, (70, 90), (32, 60), 0, 0, 360, (80, 70, 120), -1
        )
        round_result = MEASURE.shape_measurement(circle)
        elongated_result = MEASURE.shape_measurement(elongated)
        self.assertLess(round_result["aspect_ratio"], 1.15)
        self.assertGreater(elongated_result["aspect_ratio"], 1.5)
        self.assertGreater(
            round_result["circularity"], elongated_result["circularity"]
        )

    def test_orientation_angle_is_flip_stable(self) -> None:
        context = np.full((240, 240, 3), (225, 225, 225), dtype=np.uint8)
        cv2.line(context, (120, 220), (120, 120), (35, 100, 35), 8)
        cv2.ellipse(
            context, (120, 88), (36, 32), 0, 0, 360, (55, 100, 55), -1
        )
        cv2.circle(context, (120, 68), 22, (180, 60, 220), -1)
        head = context[45:122, 82:159].copy()

        original = MEASURE.orientation_measurement(head, context)
        flipped = MEASURE.orientation_measurement(
            cv2.flip(head, 1), cv2.flip(context, 1)
        )
        self.assertLess(original["angle_degrees"], 10)
        self.assertLess(
            abs(original["angle_degrees"] - flipped["angle_degrees"]), 3
        )

    def test_rebuild_helpers_preserve_box_and_annotation_identity(self) -> None:
        from PIL import Image

        image = Image.new("RGB", (100, 80), "white")
        row = {
            "audit_id": "ch1atlas_0000001",
            "queue_id": "ch1atlas_0000001",
            "det_index": "0",
            "bbox_x1": "20.2",
            "bbox_y1": "15.8",
            "bbox_x2": "60.1",
            "bbox_y2": "55.2",
        }
        self.assertEqual(
            REBUILD.annotation_unit_id(row), "ch1atlas_0000001_head_01"
        )
        head, context = REBUILD.crop_pair(image, row, 1.5)
        self.assertEqual(head.size, (41, 41))
        self.assertGreaterEqual(context.size[0], head.size[0])
        self.assertGreaterEqual(context.size[1], head.size[1])

    def test_aggregation_outputs_observation_and_species_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            measurements = pd.DataFrame([
                {
                    "annotation_unit_id": f"head_{index}",
                    "colour_status": "usable",
                    "shape_status": "usable",
                    "orientation_status": "usable",
                    "corolla_visible_fraction": 0.2,
                    "corolla_lab_lightness": 50 + index,
                    "corolla_lab_chroma": 30 + index,
                    "corolla_hue_sin": 0.5,
                    "corolla_hue_cos": 0.5,
                    "corolla_white_fraction": 0.1,
                    "corolla_redmagenta_fraction": 0.6,
                    "corolla_purple_fraction": 0.2,
                    "corolla_yellow_fraction": 0.1,
                    "shape_aspect_ratio": 1.0 + index * 0.1,
                    "shape_circularity": 0.8,
                    "shape_solidity": 0.9,
                    "shape_width_cv": 0.1,
                    "shape_narrow_end_width": 0.5,
                    "shape_middle_width": 1.0,
                    "orientation_angle_degrees": 20 + index,
                }
                for index in range(4)
            ])
            metadata = pd.DataFrame([
                {
                    "annotation_unit_id": f"head_{index}",
                    "obs_id": f"obs_{index // 2}",
                    "photo_id": f"photo_{index // 2}",
                    "taxon_name": "Cirsium_test",
                    "latitude": "35.0",
                    "longitude": "135.0",
                }
                for index in range(4)
            ])
            measurements.to_csv(root / "measurements.csv", index=False)
            metadata.to_csv(root / "metadata.csv", index=False)
            output = root / "output"
            argv = [
                "aggregate",
                "--head-measurements", str(root / "measurements.csv"),
                "--head-metadata", str(root / "metadata.csv"),
                "--out-dir", str(output),
                "--min-observations-per-species", "2",
            ]
            with patch.object(sys, "argv", argv):
                AGGREGATE.main()
            observation = pd.read_csv(
                output / "primary_trait_continuous_observation_level.csv"
            )
            species = pd.read_csv(
                output / "primary_trait_continuous_species_level.csv"
            )
            self.assertEqual(len(observation), 2)
            self.assertEqual(len(species), 1)
            self.assertEqual(
                int(species.iloc[0]["n_observations_per_species"]), 2
            )


if __name__ == "__main__":
    unittest.main()
