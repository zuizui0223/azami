from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "ch1_global" / "v2" / "79_evaluate_capitulum_detector.py"
SPEC = importlib.util.spec_from_file_location("detector_evaluation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_iou_identity_and_disjoint() -> None:
    assert MODULE.intersection_over_union((0, 0, 10, 10), (0, 0, 10, 10)) == pytest.approx(1.0)
    assert MODULE.intersection_over_union((0, 0, 10, 10), (20, 20, 30, 30)) == pytest.approx(0.0)


def test_iou_partial_overlap() -> None:
    observed = MODULE.intersection_over_union((0, 0, 10, 10), (5, 5, 15, 15))
    assert observed == pytest.approx(25 / 175)


def test_match_image_is_one_to_one() -> None:
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
    assert outcomes.count("true_positive") == 1
    assert outcomes.count("false_positive") == 1
    assert outcomes.count("false_negative") == 1
