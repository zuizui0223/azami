#!/usr/bin/env python3
"""Build a blinded, prediction-stratified review packet after detector inference.

This script does not create labels and does not alter the original 1,000-image
source-image audit manifest. It only prioritizes a compact, reproducible subset
for later independent assessment of detector errors and calibration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


STRATA = (
    "no_detection",
    "multiple_detections",
    "single_low_confidence",
    "single_mid_confidence",
    "single_high_confidence",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stratify detector predictions into a blinded review subset.")
    parser.add_argument("--screening", required=True)
    parser.add_argument("--blinded-manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-review-images", type=int, default=200)
    parser.add_argument("--low-confidence", type=float, default=0.40)
    parser.add_argument("--high-confidence", type=float, default=0.75)
    parser.add_argument("--seed", default="20260704")
    return parser.parse_args()


def stable_rank(value: str, seed: str) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def number(value: Any) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return float("nan")


def prediction_stratum(row: pd.Series, low: float, high: float) -> str:
    status = str(row.get("screen_status", ""))
    n = int(number(row.get("n_detections", 0)) or 0)
    confidence = number(row.get("max_yolo_conf", ""))
    if status == "no_detection" or n == 0:
        return "no_detection"
    if n >= 2:
        return "multiple_detections"
    if confidence < low:
        return "single_low_confidence"
    if confidence < high:
        return "single_mid_confidence"
    return "single_high_confidence"


def main() -> None:
    args = parse_args()
    if args.n_review_images < 1:
        raise ValueError("--n-review-images must be >= 1")
    if not 0 <= args.low_confidence < args.high_confidence <= 1:
        raise ValueError("confidence thresholds must satisfy 0 <= low < high <= 1")
    screening = pd.read_csv(args.screening, dtype=str, keep_default_na=False)
    manifest = pd.read_csv(args.blinded_manifest, dtype=str, keep_default_na=False)
    required_screening = {"queue_id", "screen_status", "n_detections", "max_yolo_conf", "source_image"}
    required_manifest = {"audit_id", "queue_id", "photo_id", "source_file_stem", "screen_download_filename"}
    if missing := required_screening.difference(screening.columns):
        raise ValueError(f"Screening results missing columns: {sorted(missing)}")
    if missing := required_manifest.difference(manifest.columns):
        raise ValueError(f"Blinded manifest missing columns: {sorted(missing)}")
    if screening["queue_id"].duplicated().any() or manifest["queue_id"].duplicated().any():
        raise ValueError("queue_id must be unique in screening and blinded manifest")
    if not screening["screen_status"].isin({"detected", "no_detection"}).all():
        bad = screening.loc[~screening["screen_status"].isin({"detected", "no_detection"}), "screen_status"].unique().tolist()
        raise ValueError(f"Cannot triage incomplete inference rows: {bad}")

    merged = manifest.merge(screening, on="queue_id", how="inner", validate="one_to_one")
    if len(merged) != len(manifest):
        raise ValueError("Some blinded audit images have no detector prediction")
    merged["prediction_stratum"] = merged.apply(prediction_stratum, axis=1, low=args.low_confidence, high=args.high_confidence)
    merged["selection_rank"] = merged["queue_id"].map(lambda value: stable_rank(str(value), args.seed))

    base = args.n_review_images // len(STRATA)
    remainder = args.n_review_images % len(STRATA)
    selected_parts: list[pd.DataFrame] = []
    selected_ids: set[str] = set()
    for position, stratum in enumerate(STRATA):
        n_target = base + (1 if position < remainder else 0)
        part = merged.loc[merged["prediction_stratum"].eq(stratum)].sort_values("selection_rank").head(n_target).copy()
        selected_parts.append(part)
        selected_ids.update(part["queue_id"].astype(str))
    selected = pd.concat(selected_parts, ignore_index=True) if selected_parts else merged.iloc[0:0].copy()
    remaining_needed = args.n_review_images - len(selected)
    if remaining_needed > 0:
        fill = merged.loc[~merged["queue_id"].astype(str).isin(selected_ids)].sort_values("selection_rank").head(remaining_needed).copy()
        selected = pd.concat([selected, fill], ignore_index=True)
    selected = selected.sort_values(["prediction_stratum", "selection_rank"]).copy()

    all_columns = [
        "audit_id", "queue_id", "photo_id", "source_file_stem", "screen_download_filename",
        "screen_status", "n_detections", "max_yolo_conf", "prediction_stratum", "source_image",
    ]
    review_columns = all_columns + ["selection_rank"]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    merged[all_columns].to_csv(out_dir / "bootstrap_detector_prediction_strata.csv", index=False, encoding="utf-8-sig")
    selected[review_columns].to_csv(out_dir / "bootstrap_detector_review_subset.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_predicted_images": int(len(merged)),
        "n_selected_for_review": int(len(selected)),
        "selection_counts_by_stratum": {stratum: int((selected["prediction_stratum"] == stratum).sum()) for stratum in STRATA},
        "all_prediction_counts_by_stratum": {stratum: int((merged["prediction_stratum"] == stratum).sum()) for stratum in STRATA},
        "low_confidence": args.low_confidence,
        "high_confidence": args.high_confidence,
        "seed": args.seed,
        "interpretation": "Strata and review selection are model-derived triage only. They are not detector labels, error rates, or accuracy estimates.",
    }
    (out_dir / "bootstrap_detector_triage_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
