#!/usr/bin/env python3
"""Validate double annotations and create an adjudicated detector ground truth.

Single-labelled tasks use annotator_1. Double-labelled tasks are accepted when
status and one-to-one boxes agree at the predeclared IoU threshold. Disagreements
are exported for blinded adjudication; canonical ground truth is written only when
all assigned tasks are resolved.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from azami_ch1.detection_audit import greedy_match, validate_box
except ImportError:
    module_path = Path(__file__).resolve().parents[1] / "src" / "azami_ch1" / "detection_audit.py"
    spec = importlib.util.spec_from_file_location("azami_detection_audit", module_path)
    if spec is None or spec.loader is None:
        raise
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    greedy_match = module.greedy_match
    validate_box = module.validate_box

REQUIRED = {
    "audit_id", "annotator_id", "annotation_status", "assessability",
    "image_quality", "image_width", "image_height", "gt_index", "gt_label",
    "life_stage", "occlusion", "edge_truncated", "gt_x1", "gt_y1", "gt_x2",
    "gt_y2", "notes",
}
VALID_ASSESSABILITY = {"assessable", "no_target", "unassessable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--assignments", required=True)
    parser.add_argument("--annotations", action="append", required=True)
    parser.add_argument("--adjudication", default="")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--iou-threshold", type=float, default=0.5)
    parser.add_argument("--allow-unresolved", action="store_true")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes", "y"}


def read_annotation_files(paths: list[str]) -> pd.DataFrame:
    tables = [pd.read_csv(path, dtype=str, keep_default_na=False) for path in paths]
    if not tables:
        raise ValueError("At least one annotation file is required")
    combined = pd.concat(tables, ignore_index=True)
    missing = REQUIRED.difference(combined.columns)
    if missing:
        raise ValueError(f"Annotation table missing columns: {sorted(missing)}")
    combined = combined.drop_duplicates().copy()
    combined["audit_id"] = combined["audit_id"].map(text)
    combined["annotator_id"] = combined["annotator_id"].map(text)
    if combined[["audit_id", "annotator_id"]].eq("").any().any():
        raise ValueError("audit_id and annotator_id cannot be empty")
    return combined


def validate_group(group: pd.DataFrame) -> dict[str, Any]:
    audit_id = text(group.iloc[0]["audit_id"])
    annotator = text(group.iloc[0]["annotator_id"])
    if not group["annotation_status"].map(text).eq("complete").all():
        raise ValueError(f"Incomplete annotation: {audit_id}/{annotator}")
    states = set(group["assessability"].map(text))
    if len(states) != 1 or next(iter(states)) not in VALID_ASSESSABILITY:
        raise ValueError(f"Invalid assessability state: {audit_id}/{annotator}: {states}")
    state = next(iter(states))
    dimensions = group[["image_width", "image_height"]].drop_duplicates()
    if len(dimensions) != 1:
        raise ValueError(f"Image dimensions differ within task: {audit_id}/{annotator}")
    width = int(float(dimensions.iloc[0]["image_width"]))
    height = int(float(dimensions.iloc[0]["image_height"]))
    if width < 1 or height < 1:
        raise ValueError(f"Invalid image dimensions: {audit_id}/{annotator}")
    object_rows = group.loc[group["gt_index"].map(text).ne("")].copy()
    if state == "assessable" and object_rows.empty:
        raise ValueError(f"assessable task has no boxes: {audit_id}/{annotator}")
    if state != "assessable" and not object_rows.empty:
        raise ValueError(f"{state} task cannot have boxes: {audit_id}/{annotator}")
    if object_rows["gt_index"].duplicated().any():
        raise ValueError(f"Duplicate gt_index: {audit_id}/{annotator}")
    boxes = []
    for _, row in object_rows.sort_values("gt_index").iterrows():
        box = validate_box((row["gt_x1"], row["gt_y1"], row["gt_x2"], row["gt_y2"]))
        x1, y1, x2, y2 = box
        if x1 < 0 or y1 < 0 or x2 > width or y2 > height:
            raise ValueError(f"Box outside image: {audit_id}/{annotator}: {box} vs {width}x{height}")
        boxes.append(box)
    return {
        "audit_id": audit_id,
        "annotator_id": annotator,
        "state": state,
        "width": width,
        "height": height,
        "boxes": boxes,
        "rows": group.copy(),
    }


def compare_double(first: dict[str, Any], second: dict[str, Any], threshold: float) -> dict[str, Any]:
    if first["state"] != second["state"]:
        return {
            "agreement": False, "reason": "assessability_disagreement",
            "n_first": len(first["boxes"]), "n_second": len(second["boxes"]),
            "matched": 0, "mean_iou": float("nan"),
        }
    if first["state"] != "assessable":
        return {
            "agreement": True, "reason": "status_agreement", "n_first": 0,
            "n_second": 0, "matched": 0, "mean_iou": float("nan"),
        }
    matches, unmatched_first, unmatched_second = greedy_match(
        first["boxes"], second["boxes"], threshold
    )
    agreement = not unmatched_first and not unmatched_second
    return {
        "agreement": agreement,
        "reason": "box_agreement" if agreement else "box_count_or_location_disagreement",
        "n_first": len(first["boxes"]), "n_second": len(second["boxes"]),
        "matched": len(matches),
        "mean_iou": float(np.mean([item[2] for item in matches])) if matches else float("nan"),
    }


def normalize_ground_truth(group: pd.DataFrame) -> list[dict[str, Any]]:
    first = group.iloc[0]
    state = text(first["assessability"])
    common = {
        "image_id": text(first["audit_id"]),
        "assessable_for_detection": state != "unassessable",
        "annotation_state": state,
        "image_quality": text(first["image_quality"]),
        "image_width": int(float(first["image_width"])),
        "image_height": int(float(first["image_height"])),
        "audit_notes": text(first["notes"]),
    }
    objects = group.loc[group["gt_index"].map(text).ne("")].sort_values("gt_index")
    if objects.empty:
        return [{
            **common, "object_id": "", "x1": "", "y1": "", "x2": "", "y2": "",
            "life_stage": "", "occlusion": "", "edge_truncated": "",
        }]
    rows: list[dict[str, Any]] = []
    for _, row in objects.iterrows():
        rows.append(
            {
                **common,
                "object_id": f"{common['image_id']}_head_{int(float(row['gt_index'])):02d}",
                "x1": float(row["gt_x1"]), "y1": float(row["gt_y1"]),
                "x2": float(row["gt_x2"]), "y2": float(row["gt_y2"]),
                "life_stage": text(row["life_stage"]),
                "occlusion": text(row["occlusion"]),
                "edge_truncated": as_bool(row["edge_truncated"]),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    if not 0 < args.iou_threshold <= 1:
        raise ValueError("iou-threshold must be in (0,1]")
    manifest = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    assignments = pd.read_csv(args.assignments, dtype=str, keep_default_na=False)
    annotations = read_annotation_files(args.annotations)
    expected_pairs = set(zip(assignments["audit_id"].map(text), assignments["annotator_id"].map(text)))
    found_pairs = set(zip(annotations["audit_id"], annotations["annotator_id"]))
    missing_pairs = sorted(expected_pairs - found_pairs)
    unexpected_pairs = sorted(found_pairs - expected_pairs)
    if unexpected_pairs:
        raise ValueError(f"Unexpected annotation assignments: {unexpected_pairs[:10]}")

    validated: dict[tuple[str, str], dict[str, Any]] = {}
    for key, group in annotations.groupby(["audit_id", "annotator_id"], sort=True):
        validated[(text(key[0]), text(key[1]))] = validate_group(group)

    comparison_rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    final_groups: dict[str, pd.DataFrame] = {}
    manifest_ids = set(manifest["audit_id"].map(text))
    if set(assignments["audit_id"].map(text)) != manifest_ids:
        raise ValueError("Assignments and blinded manifest do not contain identical audit IDs")

    for audit_id in sorted(manifest_ids):
        first = validated.get((audit_id, "annotator_1"))
        second_required = bool(
            assignments.loc[assignments["audit_id"].map(text).eq(audit_id), "double_label"].map(as_bool).any()
        )
        if first is None:
            unresolved.append({"audit_id": audit_id, "reason": "missing_annotator_1"})
            continue
        if not second_required:
            final_groups[audit_id] = first["rows"]
            continue
        second = validated.get((audit_id, "annotator_2"))
        if second is None:
            unresolved.append({"audit_id": audit_id, "reason": "missing_annotator_2"})
            continue
        comparison = compare_double(first, second, args.iou_threshold)
        comparison_rows.append({"audit_id": audit_id, **comparison})
        if comparison["agreement"]:
            final_groups[audit_id] = first["rows"]
        else:
            unresolved.append({"audit_id": audit_id, "reason": comparison["reason"]})

    if text(args.adjudication):
        adjudication = read_annotation_files([args.adjudication])
        for audit_id, group in adjudication.groupby("audit_id", sort=True):
            result = validate_group(group)
            final_groups[text(audit_id)] = result["rows"]
        unresolved = [row for row in unresolved if row["audit_id"] not in final_groups]

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(out / "detector_human_double_annotation_agreement.csv", index=False)
    pd.DataFrame(unresolved, columns=["audit_id", "reason"]).to_csv(
        out / "detector_adjudication_required.csv", index=False
    )

    if comparison.empty:
        human_metrics = {"n_double_label_images": 0, "n_agreed_images": 0, "image_agreement": float("nan")}
    else:
        agreed = int(comparison["agreement"].astype(bool).sum())
        human_metrics = {
            "n_double_label_images": int(len(comparison)),
            "n_agreed_images": agreed,
            "image_agreement": agreed / len(comparison),
            "mean_matched_box_iou": float(pd.to_numeric(comparison["mean_iou"], errors="coerce").mean()),
        }
    human_metrics.update({
        "n_expected_annotation_assignments": len(expected_pairs),
        "n_missing_annotation_assignments": len(missing_pairs),
        "n_unresolved_images": len(unresolved),
        "iou_threshold": args.iou_threshold,
    })
    (out / "detector_human_agreement_metrics.json").write_text(
        json.dumps(human_metrics, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    if missing_pairs:
        pd.DataFrame(missing_pairs, columns=["audit_id", "annotator_id"]).to_csv(
            out / "detector_missing_annotation_assignments.csv", index=False
        )
    if unresolved and not args.allow_unresolved:
        raise SystemExit(
            f"{len(unresolved)} audit images require annotation/adjudication; see "
            f"{out / 'detector_adjudication_required.csv'}"
        )

    canonical_rows: list[dict[str, Any]] = []
    for audit_id in sorted(final_groups):
        canonical_rows.extend(normalize_ground_truth(final_groups[audit_id]))
    pd.DataFrame(canonical_rows).to_csv(
        out / "detector_independent_audit_annotations.csv", index=False
    )
    print(json.dumps(human_metrics, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
