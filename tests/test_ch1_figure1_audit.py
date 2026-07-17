from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "ch1_global" / "v2" / "79_evaluate_capitulum_detector.py"
SPEC = importlib.util.spec_from_file_location("detector_evaluation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class TestFigure1DetectorAudit(unittest.TestCase):
    def test_iou_identity_and_disjoint(self) -> None:
        self.assertAlmostEqual(
            MODULE.intersection_over_union((0, 0, 10, 10), (0, 0, 10, 10)),
            1.0,
        )
        self.assertAlmostEqual(
            MODULE.intersection_over_union((0, 0, 10, 10), (20, 20, 30, 30)),
            0.0,
        )

    def test_iou_partial_overlap(self) -> None:
        observed = MODULE.intersection_over_union((0, 0, 10, 10), (5, 5, 15, 15))
        self.assertAlmostEqual(observed, 25 / 175)

    def test_match_image_is_one_to_one(self) -> None:
        import pandas as pd

        truth = pd.DataFrame([
            {"object_id": "a", "x1": 0, "y1": 0, "x2": 10, "y2": 10},
            {"object_id": "b", "x1": 20, "y1": 20, "x2": 30, "y2": 30},
        ])
        predictions = pd.DataFrame([
            {"prediction_id": "p1", "bbox_x1": 0, "bbox_y1": 0, "bbox_x2": 10, "bbox_y2": 10},
            {"prediction_id": "p2", "bbox_x1": 0, "bbox_y1": 0, "bbox_x2": 10, "bbox_y2": 10},
        ])
        rows = MODULE.match_image("image", truth, predictions, 0.5)
        outcomes = [row["outcome"] for row in rows]
        self.assertEqual(outcomes.count("true_positive"), 1)
        self.assertEqual(outcomes.count("false_positive"), 1)
        self.assertEqual(outcomes.count("false_negative"), 1)


if __name__ == "__main__":
    unittest.main()
