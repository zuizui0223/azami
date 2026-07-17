#!/usr/bin/env python3
"""Evaluate the frozen capitulum detector against independent source-image annotations.

Annotation rows represent every independently audited source image, including images
with zero visible capitula. Prediction rows are the frozen ``yolo_crop_metadata.csv``
output. Matching is one-to-one within each image using descending IoU.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ANNOTATION_REQUIRED = {
    "image_id", "object_id", "assessable_for_detection", "x1", "y1", "x2", "y2"
}
PREDICTION_REQUIRED = {"bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate frozen capitulum detector predictions.")
    parser.add_argument("--annotations", required=True, help="Independent detector audit CSV.")
    parser.add_argument("--predictions", required=True, help="Frozen yolo_crop_metadata.csv.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--annotation-image-column", default="image_id")
    parser.add_argument("--prediction-image-column", default="queue_id")
    parser.add_argument("--iou-threshold", type=float, default=0.5)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"1", "true", "yes", "y"}


def numeric_box(row: pd.Series, names: tuple[str, str, str, str]) -> tuple[float, float, float, float]:
    values = tuple(float(row[name]) for name in names)
    x1, y1, x2, y2 = values
    if not np.isfinite(values).all() or x2 <= x1 or y2 <= y1:
        raise ValueError(f"Invalid box: {values}")
    return values


def intersection_over_union(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> float:
    ax1, ay1, ax2, ay2 = a
    bx1, by1, bx2, by2 = b
    ix1, iy1 = max(ax1, bx1), max(ay1, by1)
    ix2, iy2 = min(ax2, bx2), min(ay2, by2)
    intersection = max(0.0, ix2 - ix1) * max(0.0, iy2 - iy1)
    union = (ax2 - ax1) * (ay2 - ay1) + (bx2 - bx1) * (by2 - by1) - intersection
    return float(intersection / union) if union > 0 else 0.0


def validate_annotations(table: pd.DataFrame, image_column: str) -> pd.DataFrame:
    missing = ANNOTATION_REQUIRED.difference(table.columns)
    if image_column not in table.columns:
        missing.add(image_column)
    if missing:
        raise ValueError(f"Annotation table missing columns: {sorted(missing)}")
    work = table.copy()
    work[image_column] = work[image_column].map(text)
    if work[image_column].eq("").any():
        raise ValueError("Annotation image identifiers cannot be empty")
    work["assessable_for_detection_bool"] = work["assessable_for_detection"].map(as_bool)
    object_rows = work["object_id"].map(text).ne("")
    for _, row in work.loc[object_rows].iterrows():
        numeric_box(row, ("x1", "y1", "x2", "y2"))
    duplicated = work.loc[object_rows, [image_column, "object_id"]].duplicated()
    if duplicated.any():
        raise ValueError("Annotation object_id must be unique within image")
    return work


def validate_predictions(table: pd.DataFrame, image_column: str) -> pd.DataFrame:
    missing = PREDICTION_REQUIRED.difference(table.columns)
    if image_column not in table.columns:
        missing.add(image_column)
    if missing:
        raise ValueError(f"Prediction table missing columns: {sorted(missing)}")
    work = table.copy()
    work[image_column] = work[image_column].map(text)
    work = work.loc[work[image_column].ne("")].copy()
    for _, row in work.iterrows():
        numeric_box(row, ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"))
    work["prediction_id"] = [f"prediction_{index:08d}" for index in range(1, len(work) + 1)]
    return work


def match_image(image_id: str, annotations: pd.DataFrame, predictions: pd.DataFrame, threshold: float) -> list[dict[str, Any]]:
    truth_rows = annotations.loc[annotations["object_id"].map(text).ne("")].copy()
    candidates: list[tuple[float, int, int]] = []
    truth_boxes = [numeric_box(row, ("x1", "y1", "x2", "y2")) for _, row in truth_rows.iterrows()]
    pred_boxes = [numeric_box(row, ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")) for _, row in predictions.iterrows()]
    for truth_index, truth_box in enumerate(truth_boxes):
        for pred_index, pred_box in enumerate(pred_boxes):
            candidates.append((intersection_over_union(truth_box, pred_box), truth_index, pred_index))
    candidates.sort(reverse=True)
    used_truth: set[int] = set()
    used_predictions: set[int] = set()
    rows: list[dict[str, Any]] = []
    for iou, truth_index, pred_index in candidates:
        if iou < threshold or truth_index in used_truth or pred_index in used_predictions:
            continue
        used_truth.add(truth_index)
        used_predictions.add(pred_index)
        truth = truth_rows.iloc[truth_index]
        prediction = predictions.iloc[pred_index]
        rows.append({"image_id": image_id, "outcome": "true_positive", "object_id": text(truth["object_id"]), "prediction_id": text(prediction["prediction_id"]), "iou": iou, "error_category": ""})
    for truth_index, (_, truth) in enumerate(truth_rows.iterrows()):
        if truth_index not in used_truth:
            rows.append({"image_id": image_id, "outcome": "false_negative", "object_id": text(truth["object_id"]), "prediction_id": "", "iou": "", "error_category": text(truth.get("missed_head_category", "unclassified_miss")) or "unclassified_miss"})
    for pred_index, (_, prediction) in enumerate(predictions.iterrows()):
        if pred_index not in used_predictions:
            rows.append({"image_id": image_id, "outcome": "false_positive", "object_id": "", "prediction_id": text(prediction["prediction_id"]), "iou": "", "error_category": text(prediction.get("false_positive_category", "unclassified_false_positive")) or "unclassified_false_positive"})
    return rows


def main() -> None:
    args = parse_args()
    if not 0 < args.iou_threshold <= 1:
        raise ValueError("--iou-threshold must be in (0, 1]")
    annotations = validate_annotations(pd.read_csv(args.annotations, dtype=str, keep_default_na=False), args.annotation_image_column)
    predictions = validate_predictions(pd.read_csv(args.predictions, dtype=str, keep_default_na=False), args.prediction_image_column)
    assessable_images = sorted(set(annotations.loc[annotations["assessable_for_detection_bool"], args.annotation_image_column]))
    match_rows: list[dict[str, Any]] = []
    image_rows: list[dict[str, Any]] = []
    for image_id in assessable_images:
        truth = annotations.loc[annotations[args.annotation_image_column].eq(image_id)]
        pred = predictions.loc[predictions[args.prediction_image_column].eq(image_id)]
        rows = match_image(image_id, truth, pred, args.iou_threshold)
        match_rows.extend(rows)
        outcomes = pd.Series([row["outcome"] for row in rows], dtype=str)
        image_rows.append({"image_id": image_id, "n_truth": int(truth["object_id"].map(text).ne("").sum()), "n_predictions": int(len(pred)), "true_positive": int(outcomes.eq("true_positive").sum()), "false_positive": int(outcomes.eq("false_positive").sum()), "false_negative": int(outcomes.eq("false_negative").sum())})
    matches = pd.DataFrame(match_rows, columns=["image_id", "outcome", "object_id", "prediction_id", "iou", "error_category"])
    image_summary = pd.DataFrame(image_rows)
    counts = matches["outcome"].value_counts() if not matches.empty else pd.Series(dtype=int)
    tp, fp, fn = int(counts.get("true_positive", 0)), int(counts.get("false_positive", 0)), int(counts.get("false_negative", 0))
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if np.isfinite(precision + recall) and precision + recall else float("nan")
    metrics = {"iou_threshold": args.iou_threshold, "n_annotated_images": int(annotations[args.annotation_image_column].nunique()), "n_assessable_images": int(len(assessable_images)), "n_unassessable_images": int(annotations.loc[~annotations["assessable_for_detection_bool"], args.annotation_image_column].nunique()), "n_true_objects": int(annotations.loc[annotations["assessable_for_detection_bool"], "object_id"].map(text).ne("").sum()), "n_predictions": int(predictions.loc[predictions[args.prediction_image_column].isin(assessable_images)].shape[0]), "true_positive": tp, "false_positive": fp, "false_negative": fn, "precision": precision, "recall": recall, "f1": f1, "matching_rule": "one-to-one greedy matching by descending IoU within source image"}
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    matches.to_csv(out_dir / "detector_object_matches.csv", index=False, encoding="utf-8-sig")
    image_summary.to_csv(out_dir / "detector_image_summary.csv", index=False, encoding="utf-8-sig")
    errors = matches.loc[matches["outcome"].ne("true_positive")]
    errors.groupby(["outcome", "error_category"], dropna=False).size().rename("n").reset_index().to_csv(out_dir / "detector_error_categories.csv", index=False, encoding="utf-8-sig")
    (out_dir / "detector_metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=True), encoding="utf-8")
    print(json.dumps(metrics, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
