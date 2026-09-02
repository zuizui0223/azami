"""Shared matching and interval utilities for independent detector audits."""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

Box = tuple[float, float, float, float]


def validate_box(box: Sequence[float]) -> Box:
    values = tuple(float(value) for value in box)
    if len(values) != 4 or not np.isfinite(values).all():
        raise ValueError(f"Invalid box: {values}")
    x1, y1, x2, y2 = values
    if x2 <= x1 or y2 <= y1:
        raise ValueError(f"Box has non-positive extent: {values}")
    return x1, y1, x2, y2


def intersection_over_union(first: Box, second: Box) -> float:
    ax1, ay1, ax2, ay2 = validate_box(first)
    bx1, by1, bx2, by2 = validate_box(second)
    ix1, iy1 = max(ax1, bx1), max(ay1, by1)
    ix2, iy2 = min(ax2, bx2), min(ay2, by2)
    intersection = max(0.0, ix2 - ix1) * max(0.0, iy2 - iy1)
    union = (ax2 - ax1) * (ay2 - ay1) + (bx2 - bx1) * (by2 - by1) - intersection
    return float(intersection / union) if union > 0 else 0.0


def greedy_match(
    truth_boxes: Sequence[Box], prediction_boxes: Sequence[Box], threshold: float
) -> tuple[list[tuple[int, int, float]], list[int], list[int]]:
    if not 0 < threshold <= 1:
        raise ValueError("IoU threshold must be in (0,1]")
    candidates: list[tuple[float, int, int]] = []
    for truth_index, truth in enumerate(truth_boxes):
        for prediction_index, prediction in enumerate(prediction_boxes):
            candidates.append(
                (intersection_over_union(truth, prediction), truth_index, prediction_index)
            )
    candidates.sort(reverse=True)
    used_truth: set[int] = set()
    used_predictions: set[int] = set()
    matches: list[tuple[int, int, float]] = []
    for iou, truth_index, prediction_index in candidates:
        if iou < threshold:
            break
        if truth_index in used_truth or prediction_index in used_predictions:
            continue
        used_truth.add(truth_index)
        used_predictions.add(prediction_index)
        matches.append((truth_index, prediction_index, float(iou)))
    unmatched_truth = [index for index in range(len(truth_boxes)) if index not in used_truth]
    unmatched_predictions = [
        index for index in range(len(prediction_boxes)) if index not in used_predictions
    ]
    return matches, unmatched_truth, unmatched_predictions


def precision_recall_f1(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = (
        2 * precision * recall / (precision + recall)
        if np.isfinite(precision) and np.isfinite(recall) and precision + recall
        else float("nan")
    )
    return float(precision), float(recall), float(f1)


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    estimate = successes / total
    denominator = 1 + z * z / total
    centre = estimate + z * z / (2 * total)
    spread = z * math.sqrt(estimate * (1 - estimate) / total + z * z / (4 * total * total))
    return (centre - spread) / denominator, (centre + spread) / denominator


def classify_relative_area(box: Box, width: float, height: float) -> str:
    if width <= 0 or height <= 0:
        return "unknown"
    x1, y1, x2, y2 = validate_box(box)
    fraction = ((x2 - x1) * (y2 - y1)) / (width * height)
    if fraction < 0.01:
        return "small_lt_1pct"
    if fraction < 0.05:
        return "medium_1_to_5pct"
    return "large_ge_5pct"
