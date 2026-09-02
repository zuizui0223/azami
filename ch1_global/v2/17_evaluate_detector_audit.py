#!/usr/bin/env python3
"""Evaluate detector predictions against adjudicated visible-capitulum boxes.

Annotations are source-image boxes, not trait labels. Images marked
unassessable are excluded from accuracy denominators rather than counted as
negative. Detection matching uses confidence-ranked greedy IoU matching.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ANNOTATION_REQUIRED = {"audit_id", "annotation_status", "assessability", "gt_index", "gt_label", "gt_x1", "gt_y1", "gt_x2", "gt_y2"}
PREDICTION_REQUIRED = {"queue_id", "yolo_conf", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate the Ch.1 visible-capitulum detector audit.")
    parser.add_argument("--annotations", required=True, help="Final adjudicated detector_audit annotations CSV")
    parser.add_argument("--predictions", required=True, help="YOLO yolo_crop_metadata.csv from the detector audit queue")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--iou-threshold", type=float, default=0.5)
    parser.add_argument("--min-yolo-conf", type=float, default=0.25)
    parser.add_argument("--audit-key", default=None, help="Optional private key for diagnostic taxon summaries")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def iou(box_a: tuple[float, float, float, float], box_b: tuple[float, float, float, float]) -> float:
    ax1, ay1, ax2, ay2 = box_a
    bx1, by1, bx2, by2 = box_b
    left, top = max(ax1, bx1), max(ay1, by1)
    right, bottom = min(ax2, bx2), min(ay2, by2)
    intersection = max(0.0, right - left) * max(0.0, bottom - top)
    area_a = max(0.0, ax2 - ax1) * max(0.0, ay2 - ay1)
    area_b = max(0.0, bx2 - bx1) * max(0.0, by2 - by1)
    denominator = area_a + area_b - intersection
    return intersection / denominator if denominator > 0 else 0.0


def metric(numerator: int, denominator: int) -> float | None:
    return None if denominator == 0 else numerator / denominator


def valid_box(row: pd.Series, prefix: str) -> tuple[float, float, float, float] | None:
    values = []
    for suffix in ("x1", "y1", "x2", "y2"):
        try:
            values.append(float(row[f"{prefix}_{suffix}"]))
        except (KeyError, TypeError, ValueError):
            return None
    x1, y1, x2, y2 = values
    return (x1, y1, x2, y2) if x2 > x1 and y2 > y1 else None


def greedy_match(gt_boxes: list[tuple[float, float, float, float]], pred_boxes: list[tuple[float, float, float, float, float]], threshold: float) -> tuple[int, int, int]:
    matched_gt: set[int] = set()
    true_positive = 0
    for *coordinates, _confidence in sorted(pred_boxes, key=lambda item: item[-1], reverse=True):
        best_index, best_iou = None, 0.0
        pred_box = tuple(coordinates)
        for index, gt_box in enumerate(gt_boxes):
            if index in matched_gt:
                continue
            value = iou(pred_box, gt_box)
            if value > best_iou:
                best_index, best_iou = index, value
        if best_index is not None and best_iou >= threshold:
            matched_gt.add(best_index)
            true_positive += 1
    false_positive = len(pred_boxes) - true_positive
    false_negative = len(gt_boxes) - true_positive
    return true_positive, false_positive, false_negative


def main() -> None:
    args = parse_args()
    if not 0 < args.iou_threshold <= 1 or not 0 <= args.min_yolo_conf <= 1:
        raise ValueError("iou-threshold must be in (0,1] and min-yolo-conf in [0,1]")
    annotations = pd.read_csv(args.annotations, dtype=str, keep_default_na=False)
    predictions = pd.read_csv(args.predictions, dtype=str, keep_default_na=False)
    missing = ANNOTATION_REQUIRED.difference(annotations.columns)
    if missing:
        raise ValueError(f"Annotations missing required columns: {sorted(missing)}")
    missing = PREDICTION_REQUIRED.difference(predictions.columns)
    if missing:
        raise ValueError(f"Predictions missing required columns: {sorted(missing)}")
    annotations["audit_id"] = annotations["audit_id"].map(text)
    annotations["annotation_status"] = annotations["annotation_status"].map(text).str.lower()
    annotations["assessability"] = annotations["assessability"].map(text).str.lower()
    gold = annotations.loc[
        annotations["annotation_status"].isin({"accepted", "adjudicated"})
        & annotations["assessability"].eq("assessable")
    ].copy()
    if gold.empty:
        raise ValueError("No accepted/adjudicated assessable annotation rows available for evaluation")
    predictions["queue_id"] = predictions["queue_id"].map(text)
    predictions["yolo_conf"] = pd.to_numeric(predictions["yolo_conf"], errors="coerce")
    predictions = predictions.loc[predictions["yolo_conf"].ge(args.min_yolo_conf)].copy()

    audit_ids = sorted(gold["audit_id"].unique())
    image_rows: list[dict[str, object]] = []
    total_tp = total_fp = total_fn = 0
    image_tp = image_fp = image_tn = image_fn = 0
    for audit_id in audit_ids:
        gt_part = gold.loc[gold["audit_id"].eq(audit_id)]
        gt_boxes = [
            box for _, row in gt_part.loc[gt_part["gt_label"].map(text).eq("visible_capitulum")].iterrows()
            if (box := valid_box(row, "gt")) is not None
        ]
        pred_part = predictions.loc[predictions["queue_id"].eq(audit_id)]
        pred_boxes = [
            (*box, float(row["yolo_conf"])) for _, row in pred_part.iterrows()
            if (box := valid_box(row, "bbox")) is not None
        ]
        tp, fp, fn = greedy_match(gt_boxes, pred_boxes, args.iou_threshold)
        total_tp += tp
        total_fp += fp
        total_fn += fn
        gt_present, pred_present = bool(gt_boxes), bool(pred_boxes)
        if gt_present and pred_present:
            image_tp += 1
        elif not gt_present and pred_present:
            image_fp += 1
        elif not gt_present and not pred_present:
            image_tn += 1
        else:
            image_fn += 1
        image_rows.append({
            "audit_id": audit_id,
            "n_gt_visible_capitula": len(gt_boxes),
            "n_yolo_predictions": len(pred_boxes),
            "detector_tp": tp,
            "detector_fp": fp,
            "detector_fn": fn,
            "gt_head_present": gt_present,
            "yolo_head_present": pred_present,
        })
    per_image = pd.DataFrame(image_rows)
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_assessable_images": int(len(per_image)),
        "iou_threshold": args.iou_threshold,
        "min_yolo_conf": args.min_yolo_conf,
        "box_true_positive": int(total_tp),
        "box_false_positive": int(total_fp),
        "box_false_negative": int(total_fn),
        "box_precision": metric(total_tp, total_tp + total_fp),
        "box_recall": metric(total_tp, total_tp + total_fn),
        "image_presence_true_positive": int(image_tp),
        "image_presence_false_positive": int(image_fp),
        "image_presence_true_negative": int(image_tn),
        "image_presence_false_negative": int(image_fn),
        "image_presence_precision": metric(image_tp, image_tp + image_fp),
        "image_presence_recall": metric(image_tp, image_tp + image_fn),
        "note": "Metrics apply only to adjudicated assessable source images. This audit evaluates visible-capitulum detection, not individual-level flower number or morphology traits.",
    }
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.audit_key:
        key = pd.read_csv(args.audit_key, dtype=str, keep_default_na=False)
        if {"audit_id", "taxon_name"}.issubset(key.columns):
            keyed = per_image.merge(key[["audit_id", "taxon_name", "spatial_block"]], on="audit_id", how="left")
            keyed.groupby("taxon_name", dropna=False).agg(
                n_images=("audit_id", "size"),
                mean_gt_count=("n_gt_visible_capitula", "mean"),
                mean_yolo_count=("n_yolo_predictions", "mean"),
                detector_fn=("detector_fn", "sum"),
            ).reset_index().to_csv(out_dir / "detector_audit_taxon_diagnostics.csv", index=False, encoding="utf-8-sig")
    per_image.to_csv(out_dir / "detector_audit_per_image_metrics.csv", index=False, encoding="utf-8-sig")
    (out_dir / "detector_audit_metrics.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
