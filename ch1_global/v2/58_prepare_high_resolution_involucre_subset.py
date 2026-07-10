#!/usr/bin/env python3
"""Select high-resolution detected heads for exploratory involucre measurements."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--yolo-metadata", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--min-bbox-dimension", type=float, default=150.0)
    p.add_argument("--min-yolo-confidence", type=float, default=0.0)
    return p.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def annotation_unit(row: pd.Series) -> str:
    queue = text(row.get("audit_id")) or text(row.get("queue_id"))
    if not queue:
        raise ValueError("YOLO row lacks audit_id/queue_id")
    return f"{queue}_head_{int(float(row['det_index'])) + 1:02d}"


def main() -> None:
    args = parse_args()
    table = pd.read_csv(args.yolo_metadata, dtype=str, keep_default_na=False)
    required = {
        "queue_id", "obs_id", "photo_id", "taxon_name", "det_index", "yolo_conf",
        "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "medium_image_url",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"YOLO metadata missing columns: {sorted(missing)}")
    table = table.copy()
    for column in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "yolo_conf"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["bbox_width"] = table["bbox_x2"] - table["bbox_x1"]
    table["bbox_height"] = table["bbox_y2"] - table["bbox_y1"]
    table["bbox_min_dimension"] = table[["bbox_width", "bbox_height"]].min(axis=1)
    table["annotation_unit_id"] = table.apply(annotation_unit, axis=1)
    if table["annotation_unit_id"].duplicated().any():
        raise ValueError("Derived annotation_unit_id values are not unique")

    selected = table.loc[
        table["bbox_min_dimension"].ge(args.min_bbox_dimension)
        & table["yolo_conf"].ge(args.min_yolo_confidence)
        & table["medium_image_url"].map(text).ne("")
    ].copy()
    selected = selected.sort_values(["photo_id", "det_index"])
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out_csv, index=False, encoding="utf-8-sig")

    report = {
        "n_input_heads": int(len(table)),
        "n_selected_heads": int(len(selected)),
        "n_selected_photos": int(selected["photo_id"].nunique()),
        "n_selected_species": int(selected["taxon_name"].nunique()),
        "min_bbox_dimension": args.min_bbox_dimension,
        "min_yolo_confidence": args.min_yolo_confidence,
        "semantic_status": (
            "Resolution prefilter for exploratory involucre/bract continuous proxies. "
            "Not a claim that excluded images lack spines, hairs or spreading bracts."
        ),
    }
    Path(args.report).write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
