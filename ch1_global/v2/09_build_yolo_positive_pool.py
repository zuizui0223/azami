#!/usr/bin/env python3
"""Build image- and head-level outputs from detector-only YOLO screening.

The YOLO detector is the operational gate for whether a photograph contains a
visible Cirsium head. It also yields a *visible-head count per photograph*.
That count is a photographic observation (affected by framing and occlusion),
not a plant-level flower-number estimate.

Outputs
-------
- yolo_image_summary.csv: every queued image, including `head_present_yolo`
  and `visible_head_count_yolo`.
- yolo_positive_heads.csv: one row per valid detected head, with the source
  image's visible-head count repeated on each crop.
- yolo_positive_summary.json: compact audit counts.

Example
-------
python ch1_global/v2/09_build_yolo_positive_pool.py \
  --queue "queue/image_screening_queue.csv" \
  --screening "yolo/screening_results.csv" \
  --crops "yolo/yolo_crop_metadata.csv" \
  --out-dir "yolo/positive_pool"
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Make explicit image- and head-level YOLO-positive datasets.")
    parser.add_argument("--queue", required=True, help="image_screening_queue.csv")
    parser.add_argument("--screening", required=True, help="screening_results.csv")
    parser.add_argument("--crops", required=True, help="yolo_crop_metadata.csv")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-yolo-conf", type=float, default=0.25)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.min_yolo_conf <= 1.0:
        raise ValueError("--min-yolo-conf must be in [0, 1]")

    queue = pd.read_csv(args.queue, dtype=str, keep_default_na=False)
    screening = pd.read_csv(args.screening, dtype=str, keep_default_na=False)
    crop_path = Path(args.crops)
    crops = pd.read_csv(crop_path, dtype=str, keep_default_na=False) if crop_path.exists() and crop_path.stat().st_size else pd.DataFrame()
    required_queue = {"queue_id", "obs_id", "photo_id", "taxon_name"}
    required_screening = {"queue_id", "screen_status", "n_detections", "max_yolo_conf"}
    missing = required_queue.difference(queue.columns)
    if missing:
        raise ValueError(f"Queue missing columns: {sorted(missing)}")
    missing = required_screening.difference(screening.columns)
    if missing:
        raise ValueError(f"Screening results missing columns: {sorted(missing)}")
    if queue["queue_id"].duplicated().any():
        raise ValueError("Queue has duplicate queue_id values")

    screening["_row_order"] = range(len(screening))
    latest = screening.sort_values("_row_order").groupby("queue_id", as_index=False).tail(1).copy()
    latest["visible_head_count_yolo"] = pd.to_numeric(latest["n_detections"], errors="coerce").fillna(0).astype(int)
    latest["max_yolo_conf"] = pd.to_numeric(latest["max_yolo_conf"], errors="coerce")
    latest["head_present_yolo"] = latest["screen_status"].eq("detected") & latest["visible_head_count_yolo"].gt(0)
    latest["yolo_screen_complete"] = latest["screen_status"].isin(["detected", "no_detection", "missing_image"])
    image_summary = queue.merge(
        latest[["queue_id", "screen_status", "visible_head_count_yolo", "max_yolo_conf", "head_present_yolo", "yolo_screen_complete"]],
        on="queue_id", how="left", validate="one_to_one",
    )
    image_summary["screen_status"] = image_summary["screen_status"].fillna("not_screened")
    image_summary["visible_head_count_yolo"] = image_summary["visible_head_count_yolo"].fillna(0).astype(int)
    image_summary["head_present_yolo"] = image_summary["head_present_yolo"].fillna(False).astype(bool)
    image_summary["yolo_screen_complete"] = image_summary["yolo_screen_complete"].fillna(False).astype(bool)
    image_summary["head_count_bin_yolo"] = pd.cut(
        image_summary["visible_head_count_yolo"],
        bins=[-1, 0, 1, 3, float("inf")],
        labels=["0", "1", "2_3", "4_plus"],
    ).astype(str)

    if crops.empty:
        positive_heads = pd.DataFrame(columns=list(queue.columns) + ["visible_head_count_yolo", "head_present_yolo"])
    else:
        required_crops = {"queue_id", "yolo_conf", "crop_path"}
        missing = required_crops.difference(crops.columns)
        if missing:
            raise ValueError(f"YOLO crop metadata missing columns: {sorted(missing)}")
        crops["yolo_conf"] = pd.to_numeric(crops["yolo_conf"], errors="coerce")
        positive_heads = crops.loc[crops["yolo_conf"].ge(args.min_yolo_conf)].copy()
        positive_heads = positive_heads.merge(
            image_summary[["queue_id", "visible_head_count_yolo", "head_present_yolo", "head_count_bin_yolo", "max_yolo_conf"]],
            on="queue_id", how="inner", validate="many_to_one", suffixes=("", "_source_image"),
        )
        positive_heads = positive_heads.loc[positive_heads["head_present_yolo"]].copy()
        positive_heads["source_image_head_count"] = positive_heads["visible_head_count_yolo"].astype(int)
        positive_heads["source_image_multiple_heads"] = positive_heads["source_image_head_count"].gt(1)
        positive_heads = positive_heads.sort_values(["queue_id", "yolo_conf"], ascending=[True, False]).reset_index(drop=True)
        positive_heads["head_rank_by_yolo_conf"] = positive_heads.groupby("queue_id").cumcount() + 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    image_summary.to_csv(out_dir / "yolo_image_summary.csv", index=False, encoding="utf-8-sig")
    positive_heads.to_csv(out_dir / "yolo_positive_heads.csv", index=False, encoding="utf-8-sig")
    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "min_yolo_conf": args.min_yolo_conf,
        "n_queued_images": int(len(queue)),
        "n_screened_images": int(image_summary["yolo_screen_complete"].sum()),
        "n_yolo_positive_images": int(image_summary["head_present_yolo"].sum()),
        "n_yolo_positive_head_crops": int(len(positive_heads)),
        "n_multi_head_images": int((image_summary["visible_head_count_yolo"] > 1).sum()),
        "note": "visible_head_count_yolo counts detector-positive heads visible in each photograph. It is not a plant-level flower-number estimate without a separately validated whole-plant imaging protocol.",
    }
    (out_dir / "yolo_positive_summary.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
