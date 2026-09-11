#!/usr/bin/env python3
"""Prepare a trait-blind exhaustive high-resolution continuous-trait cohort.

The cohort starts from all detector heads, intersects the already frozen strict
taxon-by-0.25-degree-cell observation cohort, applies only image-resolution and
URL availability gates, and retains one deterministic highest-resolution head per
observation. No measured trait value or categorical state is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


FORBIDDEN_SELECTION_COLUMNS = {
    "orientation_angle_degrees_median",
    "corolla_lab_lightness_median",
    "corolla_lab_chroma_median",
    "shape_aspect_ratio_median",
    "involucre_projection_roughness",
    "involucre_spread_fraction",
    "ai_candidate_state",
    "trait_state",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--yolo-metadata", required=True)
    parser.add_argument("--strict-observations", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--min-bbox-dimension", type=float, default=150.0)
    parser.add_argument("--min-yolo-confidence", type=float, default=0.0)
    return parser.parse_args()


def _text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def prepare_cohort(
    metadata: pd.DataFrame,
    strict: pd.DataFrame,
    min_bbox_dimension: float,
    min_yolo_confidence: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {
        "queue_id", "obs_id", "photo_id", "taxon_name", "det_index", "yolo_conf",
        "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "medium_image_url",
    }
    missing = sorted(required.difference(metadata.columns))
    if missing:
        raise ValueError(f"YOLO metadata are missing columns: {missing}")
    if "obs_id" not in strict:
        raise ValueError("Strict observation cohort is missing obs_id")
    if strict["obs_id"].astype(str).duplicated().any():
        raise ValueError("Strict observation cohort must be unique by obs_id")
    leaked = sorted(FORBIDDEN_SELECTION_COLUMNS.intersection(metadata.columns))
    if leaked:
        raise ValueError(f"Trait-derived fields are forbidden in cohort selection: {leaked}")

    table = metadata.copy()
    table["obs_id"] = table["obs_id"].astype(str)
    for column in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "yolo_conf", "det_index"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["bbox_width"] = table["bbox_x2"] - table["bbox_x1"]
    table["bbox_height"] = table["bbox_y2"] - table["bbox_y1"]
    table["bbox_min_dimension"] = table[["bbox_width", "bbox_height"]].min(axis=1)
    table["annotation_unit_id"] = [
        f"{_text(queue)}_head_{int(index) + 1:02d}" if pd.notna(index) else ""
        for queue, index in zip(table["queue_id"], table["det_index"])
    ]
    if table["annotation_unit_id"].eq("").any() or table["annotation_unit_id"].duplicated().any():
        raise ValueError("Derived annotation_unit_id values are blank or duplicated")

    high_resolution = table[
        table["bbox_min_dimension"].ge(min_bbox_dimension)
        & table["yolo_conf"].ge(min_yolo_confidence)
        & table["medium_image_url"].map(_text).ne("")
    ].copy()
    strict_ids = set(strict["obs_id"].astype(str))
    eligible = high_resolution[high_resolution["obs_id"].isin(strict_ids)].copy()
    selected = (
        eligible.sort_values(
            ["obs_id", "bbox_min_dimension", "yolo_conf", "photo_id", "det_index"],
            ascending=[True, False, False, True, True],
            kind="mergesort",
        )
        .drop_duplicates("obs_id", keep="first")
        .sort_values(["taxon_name", "obs_id", "annotation_unit_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    if selected["obs_id"].duplicated().any():
        raise RuntimeError("One-head-per-observation selection failed")
    report = {
        "n_input_heads": int(len(table)),
        "n_input_photos": int(table["photo_id"].nunique()),
        "n_input_observations": int(table["obs_id"].nunique()),
        "n_high_resolution_heads": int(len(high_resolution)),
        "n_high_resolution_observations": int(high_resolution["obs_id"].nunique()),
        "n_strict_observations": int(len(strict)),
        "n_high_resolution_heads_in_strict_cohort": int(len(eligible)),
        "n_selected_heads": int(len(selected)),
        "n_selected_observations": int(selected["obs_id"].nunique()),
        "n_selected_taxa": int(selected["taxon_name"].nunique()),
        "min_bbox_dimension": float(min_bbox_dimension),
        "min_yolo_confidence": float(min_yolo_confidence),
        "selection_rule": "strict spatial cohort intersect high-resolution heads; one largest bbox-min-dimension head per observation; confidence and stable ids break ties",
        "trait_blind": True,
    }
    return selected, report


def main() -> None:
    args = parse_args()
    metadata = pd.read_csv(args.yolo_metadata, dtype=str, keep_default_na=False, low_memory=False)
    strict = pd.read_csv(args.strict_observations, dtype=str, keep_default_na=False, low_memory=False)
    selected, report = prepare_cohort(
        metadata, strict, args.min_bbox_dimension, args.min_yolo_confidence
    )
    output = Path(args.out_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output, index=False, encoding="utf-8-sig")
    Path(args.report).write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
