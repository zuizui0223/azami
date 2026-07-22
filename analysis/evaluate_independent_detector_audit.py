#!/usr/bin/env python3
"""Evaluate the frozen visible-capitulum detector on adjudicated independent labels."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from azami_ch1.detection_audit import (
        classify_relative_area, greedy_match, precision_recall_f1, validate_box,
        wilson_interval,
    )
except ImportError:
    module_path = Path(__file__).resolve().parents[1] / "src" / "azami_ch1" / "detection_audit.py"
    spec = importlib.util.spec_from_file_location("azami_detection_audit", module_path)
    if spec is None or spec.loader is None:
        raise
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    classify_relative_area = module.classify_relative_area
    greedy_match = module.greedy_match
    precision_recall_f1 = module.precision_recall_f1
    validate_box = module.validate_box
    wilson_interval = module.wilson_interval

ANNOTATION_REQUIRED = {
    "image_id", "object_id", "assessable_for_detection", "annotation_state",
    "image_quality", "image_width", "image_height", "x1", "y1", "x2", "y2",
    "life_stage", "occlusion", "edge_truncated",
}
PREDICTION_REQUIRED = {
    "queue_id", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "yolo_conf",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--private-key", default="")
    parser.add_argument("--production-confidence", type=float, default=0.25)
    parser.add_argument("--minimum-prediction-confidence", type=float, default=0.01)
    parser.add_argument("--iou-threshold", type=float, default=0.5)
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260722)
    parser.add_argument("--minimum-stratum-images", type=int, default=10)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes", "y"}


def load_annotations(path: str) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = ANNOTATION_REQUIRED.difference(table.columns)
    if missing:
        raise ValueError(f"Annotations missing columns: {sorted(missing)}")
    table["image_id"] = table["image_id"].map(text)
    table["assessable_bool"] = table["assessable_for_detection"].map(as_bool)
    table["image_width_num"] = pd.to_numeric(table["image_width"], errors="raise")
    table["image_height_num"] = pd.to_numeric(table["image_height"], errors="raise")
    object_rows = table["object_id"].map(text).ne("")
    for index, row in table.loc[object_rows].iterrows():
        validate_box((row["x1"], row["y1"], row["x2"], row["y2"]))
        if not row["assessable_bool"]:
            raise ValueError(f"Unassessable image contains target box: row {index}")
    return table


def load_predictions(path: str, minimum_confidence: float) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = PREDICTION_REQUIRED.difference(table.columns)
    if missing:
        raise ValueError(f"Predictions missing columns: {sorted(missing)}")
    table["queue_id"] = table["queue_id"].map(text)
    table["confidence"] = pd.to_numeric(table["yolo_conf"], errors="raise")
    table = table.loc[table["confidence"].ge(minimum_confidence)].copy()
    for _, row in table.iterrows():
        validate_box((row["bbox_x1"], row["bbox_y1"], row["bbox_x2"], row["bbox_y2"]))
    table["prediction_id"] = [f"prediction_{index:08d}" for index in range(1, len(table) + 1)]
    return table


def image_boxes(frame: pd.DataFrame, names: tuple[str, str, str, str]) -> list[tuple[float, float, float, float]]:
    rows = frame.loc[frame[names[0]].map(text).ne("")]
    return [validate_box(tuple(row[name] for name in names)) for _, row in rows.iterrows()]


def evaluate_threshold(
    annotations: pd.DataFrame,
    predictions: pd.DataFrame,
    threshold: float,
    iou_threshold: float,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    assessable = annotations.loc[annotations["assessable_bool"]].copy()
    image_ids = sorted(assessable["image_id"].unique())
    predictions = predictions.loc[predictions["confidence"].ge(threshold)].copy()
    object_rows: list[dict[str, Any]] = []
    image_rows: list[dict[str, Any]] = []
    for image_id in image_ids:
        truth = assessable.loc[assessable["image_id"].eq(image_id)]
        pred = predictions.loc[predictions["queue_id"].eq(image_id)]
        truth_objects = truth.loc[truth["object_id"].map(text).ne("")].reset_index(drop=True)
        pred_objects = pred.reset_index(drop=True)
        truth_boxes = image_boxes(truth_objects, ("x1", "y1", "x2", "y2"))
        pred_boxes = image_boxes(pred_objects, ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"))
        matches, unmatched_truth, unmatched_predictions = greedy_match(
            truth_boxes, pred_boxes, iou_threshold
        )
        matched_truth = {item[0]: item for item in matches}
        matched_predictions = {item[1]: item for item in matches}
        for truth_index, row in truth_objects.iterrows():
            matched = matched_truth.get(truth_index)
            box = truth_boxes[truth_index]
            object_rows.append({
                "image_id": image_id,
                "object_id": text(row["object_id"]),
                "outcome": "true_positive" if matched else "false_negative",
                "prediction_id": text(pred_objects.iloc[matched[1]]["prediction_id"]) if matched else "",
                "iou": matched[2] if matched else "",
                "confidence": float(pred_objects.iloc[matched[1]]["confidence"]) if matched else "",
                "size_class": classify_relative_area(
                    box, float(row["image_width_num"]), float(row["image_height_num"])
                ),
                "life_stage": text(row["life_stage"]) or "unknown",
                "occlusion": text(row["occlusion"]) or "unknown",
                "edge_truncated": as_bool(row["edge_truncated"]),
                "image_quality": text(row["image_quality"]) or "unknown",
            })
        for prediction_index, row in pred_objects.iterrows():
            if prediction_index in matched_predictions:
                continue
            object_rows.append({
                "image_id": image_id, "object_id": "", "outcome": "false_positive",
                "prediction_id": text(row["prediction_id"]), "iou": "",
                "confidence": float(row["confidence"]), "size_class": "prediction_only",
                "life_stage": "", "occlusion": "", "edge_truncated": "",
                "image_quality": text(truth.iloc[0]["image_quality"]) or "unknown",
                "bbox_x1": float(row["bbox_x1"]), "bbox_y1": float(row["bbox_y1"]),
                "bbox_x2": float(row["bbox_x2"]), "bbox_y2": float(row["bbox_y2"]),
            })
        image_rows.append({
            "image_id": image_id,
            "annotation_state": text(truth.iloc[0]["annotation_state"]),
            "image_quality": text(truth.iloc[0]["image_quality"]) or "unknown",
            "n_truth": len(truth_objects), "n_predictions": len(pred_objects),
            "true_positive": len(matches), "false_positive": len(unmatched_predictions),
            "false_negative": len(unmatched_truth),
        })
    objects = pd.DataFrame(object_rows)
    images = pd.DataFrame(image_rows)
    counts = objects["outcome"].value_counts() if not objects.empty else pd.Series(dtype=int)
    tp = int(counts.get("true_positive", 0))
    fp = int(counts.get("false_positive", 0))
    fn = int(counts.get("false_negative", 0))
    precision, recall, f1 = precision_recall_f1(tp, fp, fn)
    no_target = images["n_truth"].eq(0)
    positive = images["n_truth"].gt(0)
    specificity = float(images.loc[no_target, "n_predictions"].eq(0).mean()) if no_target.any() else float("nan")
    image_sensitivity = float(images.loc[positive, "true_positive"].gt(0).mean()) if positive.any() else float("nan")
    metrics = {
        "confidence_threshold": threshold,
        "iou_threshold": iou_threshold,
        "n_assessable_images": len(images),
        "n_true_objects": tp + fn,
        "n_predictions": tp + fp,
        "true_positive": tp, "false_positive": fp, "false_negative": fn,
        "precision": precision, "recall": recall, "f1": f1,
        "precision_wilson_ci95": wilson_interval(tp, tp + fp),
        "recall_wilson_ci95": wilson_interval(tp, tp + fn),
        "no_target_image_specificity": specificity,
        "positive_image_sensitivity": image_sensitivity,
    }
    return metrics, objects, images


def bootstrap_metrics(images: pd.DataFrame, repeats: int, seed: int) -> dict[str, tuple[float, float]]:
    if repeats < 1 or images.empty:
        return {}
    rng = np.random.default_rng(seed)
    values = {"precision": [], "recall": [], "f1": [], "specificity": [], "image_sensitivity": []}
    n = len(images)
    for _ in range(repeats):
        sample = images.iloc[rng.integers(0, n, size=n)]
        tp = int(sample["true_positive"].sum())
        fp = int(sample["false_positive"].sum())
        fn = int(sample["false_negative"].sum())
        precision, recall, f1 = precision_recall_f1(tp, fp, fn)
        no_target = sample["n_truth"].eq(0)
        positive = sample["n_truth"].gt(0)
        specificity = float(sample.loc[no_target, "n_predictions"].eq(0).mean()) if no_target.any() else np.nan
        sensitivity = float(sample.loc[positive, "true_positive"].gt(0).mean()) if positive.any() else np.nan
        for name, value in zip(values, (precision, recall, f1, specificity, sensitivity)):
            values[name].append(value)
    return {
        name: tuple(np.nanquantile(series, [0.025, 0.975]).astype(float))
        for name, series in values.items()
    }


def recall_strata(objects: pd.DataFrame, column: str) -> pd.DataFrame:
    truth = objects.loc[objects["outcome"].isin(["true_positive", "false_negative"])].copy()
    rows = []
    for value, part in truth.groupby(column, dropna=False):
        tp = int(part["outcome"].eq("true_positive").sum())
        total = len(part)
        low, high = wilson_interval(tp, total)
        rows.append({
            "stratum": column, "level": text(value) or "missing",
            "n_true_objects": total, "true_positive": tp,
            "recall": tp / total if total else np.nan,
            "ci95_low": low, "ci95_high": high,
        })
    return pd.DataFrame(rows)


def image_strata(images: pd.DataFrame, column: str, minimum_images: int) -> pd.DataFrame:
    rows = []
    for value, part in images.groupby(column, dropna=False):
        if len(part) < minimum_images:
            continue
        tp = int(part["true_positive"].sum())
        fp = int(part["false_positive"].sum())
        fn = int(part["false_negative"].sum())
        precision, recall, f1 = precision_recall_f1(tp, fp, fn)
        rows.append({
            "stratum": column, "level": text(value) or "missing",
            "n_images": len(part), "n_true_objects": tp + fn,
            "n_predictions": tp + fp, "precision": precision,
            "recall": recall, "f1": f1,
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    if not 0 <= args.minimum_prediction_confidence <= args.production_confidence <= 1:
        raise ValueError("Require 0 <= minimum-prediction-confidence <= production-confidence <= 1")
    if not 0 < args.iou_threshold <= 1:
        raise ValueError("iou-threshold must be in (0,1]")
    annotations = load_annotations(args.annotations)
    predictions = load_predictions(args.predictions, args.minimum_prediction_confidence)
    assessable_ids = set(annotations.loc[annotations["assessable_bool"], "image_id"])
    predictions = predictions.loc[predictions["queue_id"].isin(assessable_ids)].copy()

    key = None
    if text(args.private_key):
        key = pd.read_csv(args.private_key, dtype=str, keep_default_na=False)
        if "audit_id" not in key.columns:
            raise ValueError("private key must contain audit_id")

    metrics, objects, images = evaluate_threshold(
        annotations, predictions, args.production_confidence, args.iou_threshold
    )
    metrics["n_annotated_images"] = int(annotations["image_id"].nunique())
    metrics["n_unassessable_images"] = int(annotations.loc[~annotations["assessable_bool"], "image_id"].nunique())
    metrics["n_no_target_images"] = int(annotations.loc[annotations["annotation_state"].eq("no_target"), "image_id"].nunique())
    metrics["image_bootstrap_ci95"] = bootstrap_metrics(images, args.bootstrap_repeats, args.bootstrap_seed)
    metrics["matching_rule"] = "one-to-one greedy matching by descending IoU within source image"
    metrics["inference_boundary"] = "independent source-image audit; no training/proposal photo or observation is included"

    if key is not None:
        images = images.merge(key, left_on="image_id", right_on="audit_id", how="left", validate="one_to_one")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    objects.to_csv(out / "detector_object_outcomes.csv", index=False)
    images.to_csv(out / "detector_image_outcomes.csv", index=False)
    objects.loc[objects["outcome"].eq("false_positive")].to_csv(out / "detector_false_positive_review_queue.csv", index=False)
    objects.loc[objects["outcome"].eq("false_negative")].to_csv(out / "detector_false_negative_review_queue.csv", index=False)

    curves = []
    for threshold in np.round(np.arange(0.05, 0.951, 0.05), 2):
        row, _, _ = evaluate_threshold(annotations, predictions, float(threshold), args.iou_threshold)
        curves.append(row)
    pd.DataFrame(curves).to_csv(out / "detector_confidence_threshold_curve.csv", index=False)

    object_strata = pd.concat(
        [recall_strata(objects, column) for column in (
            "size_class", "life_stage", "occlusion", "edge_truncated", "image_quality"
        )], ignore_index=True,
    )
    object_strata.to_csv(out / "detector_object_stratified_recall.csv", index=False)
    image_tables = [image_strata(images, "image_quality", args.minimum_stratum_images)]
    if key is not None:
        for column in ("latitude_band", "spatial_block", "taxon_name"):
            if column in images.columns:
                image_tables.append(image_strata(images, column, args.minimum_stratum_images))
    pd.concat(image_tables, ignore_index=True).to_csv(
        out / "detector_image_stratified_metrics.csv", index=False
    )
    (out / "detector_independent_metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
