#!/usr/bin/env python3
"""Build the exhaustive public-photo queue for high-resolution within-species analysis.

Every eligible licensed photo is retained for automated head detection and
continuous-trait measurement. No species cap, environmental balancing, spatial
thinning or trait-based selection is applied at this stage. Any thinning or
weighting occurs only after prediction, and the unthinned predictions remain the
source-of-truth table.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

REQUIRED = {
    "obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank",
    "quality_grade", "captive", "latitude", "longitude",
    "coordinate_usable_for_environment", "positional_accuracy", "geoprivacy",
    "obscured", "observed_on", "small_image_url", "medium_image_url",
    "large_image_url", "source_file_stem", "photo_license_code",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an exhaustive Cirsium photo-inference queue.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--quality-grades", default="research")
    parser.add_argument("--photo-licenses", default="cc0,cc-by,cc-by-sa,cc-by-nc,cc-by-nc-sa")
    parser.add_argument("--require-open-coordinate", action="store_true")
    parser.add_argument("--max-positional-accuracy-m", type=float, default=0,
                        help="0 keeps every public usable coordinate; positive values impose a limit.")
    parser.add_argument("--target-species", default="", help="Optional comma-separated exact taxa.")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def finite_number(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def spatial_block(latitude: Any, longitude: Any, width: float) -> str:
    lat, lon = finite_number(latitude), finite_number(longitude)
    if not math.isfinite(lat) or not math.isfinite(lon):
        return "missing"
    return f"lat{math.floor((lat + 90.0) / width):03d}_lon{math.floor((lon + 180.0) / width):03d}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_queue(metadata: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = sorted(REQUIRED.difference(metadata.columns))
    if missing:
        raise ValueError(f"Merged metadata missing required columns: {missing}")

    quality = {item.strip().lower() for item in args.quality_grades.split(",") if item.strip()}
    licenses = {item.strip().lower() for item in args.photo_licenses.split(",") if item.strip()}
    targets = {item.strip() for item in args.target_species.split(",") if item.strip()}
    if not quality or not licenses:
        raise ValueError("At least one quality grade and photo licence are required")

    work = metadata.copy()
    work["taxon_rank"] = work["taxon_rank"].map(text).str.lower()
    work["quality_grade"] = work["quality_grade"].map(text).str.lower()
    work["photo_license_normalized"] = work["photo_license_code"].map(text).str.lower()
    work["captive_bool"] = work["captive"].map(as_bool)
    work["coordinate_usable_bool"] = work["coordinate_usable_for_environment"].map(as_bool)
    work["positional_accuracy_numeric"] = pd.to_numeric(work["positional_accuracy"], errors="coerce")

    eligible = work.loc[
        work["obs_id"].map(text).ne("")
        & work["photo_id"].map(text).ne("")
        & work["taxon_name"].map(text).ne("")
        & work["taxon_rank"].eq("species")
        & work["quality_grade"].isin(quality)
        & work["photo_license_normalized"].isin(licenses)
        & ~work["captive_bool"]
        & work["medium_image_url"].map(text).ne("")
    ].copy()
    if args.require_open_coordinate:
        eligible = eligible.loc[eligible["coordinate_usable_bool"]].copy()
    if args.max_positional_accuracy_m > 0:
        eligible = eligible.loc[
            eligible["positional_accuracy_numeric"].notna()
            & eligible["positional_accuracy_numeric"].le(args.max_positional_accuracy_m)
        ].copy()
    if targets:
        eligible = eligible.loc[eligible["taxon_name"].isin(targets)].copy()
    if eligible.empty:
        raise ValueError("No eligible photos after filtering")

    eligible = eligible.sort_values(["photo_id", "obs_id", "photo_index"]).drop_duplicates("photo_id", keep="first")
    for width in (1.0, 2.0, 5.0):
        label = str(int(width))
        eligible[f"spatial_block_{label}deg"] = [
            spatial_block(lat, lon, width) for lat, lon in zip(eligible["latitude"], eligible["longitude"])
        ]

    eligible = eligible.sort_values(["taxon_name", "obs_id", "photo_index", "photo_id"]).reset_index(drop=True)
    eligible.insert(0, "queue_id", [f"ch1allphoto_{index:09d}" for index in range(1, len(eligible) + 1)])
    eligible["audit_id"] = eligible["queue_id"]
    eligible["queue_status"] = "not_downloaded"
    eligible["sampling_role"] = "exhaustive_pre_prediction"
    eligible["screen_download_filename"] = eligible["source_file_stem"].map(lambda value: f"{text(value)}.jpg")

    audit = (
        eligible.groupby("taxon_name", sort=True)
        .agg(
            n_photos=("photo_id", "size"),
            n_observations=("obs_id", "nunique"),
            n_users=("user_id", "nunique") if "user_id" in eligible.columns else ("obs_id", "size"),
            n_blocks_1deg=("spatial_block_1deg", lambda values: sum(pd.Series(values).ne("missing").groupby(values).any())),
            n_blocks_2deg=("spatial_block_2deg", lambda values: pd.Series(values)[pd.Series(values).ne("missing")].nunique()),
            n_blocks_5deg=("spatial_block_5deg", lambda values: pd.Series(values)[pd.Series(values).ne("missing")].nunique()),
            n_coordinate_usable=("coordinate_usable_bool", "sum"),
        )
        .reset_index()
        .sort_values(["n_photos", "taxon_name"], ascending=[False, True])
    )
    return eligible, audit


def main() -> None:
    args = parse_args()
    if args.max_positional_accuracy_m < 0:
        raise ValueError("max-positional-accuracy-m cannot be negative")
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    metadata = pd.read_csv(input_path, dtype=str, keep_default_na=False, low_memory=False)
    queue, audit = build_queue(metadata, args)

    keep = [
        "queue_id", "audit_id", "queue_status", "sampling_role", "obs_id", "photo_id", "photo_index",
        "taxon_name", "taxon_rank", "quality_grade", "latitude", "longitude",
        "coordinate_usable_for_environment", "positional_accuracy", "geoprivacy", "obscured",
        "spatial_block_1deg", "spatial_block_2deg", "spatial_block_5deg", "observed_on",
        "captive", "user_id", "user_login", "photo_license_code", "photo_attribution",
        "small_image_url", "medium_image_url", "large_image_url", "source_file_stem",
        "screen_download_filename",
    ]
    for column in keep:
        if column not in queue.columns:
            queue[column] = ""

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue[keep].to_csv(out_dir / "exhaustive_photo_inference_queue.csv", index=False, encoding="utf-8-sig")
    audit.to_csv(out_dir / "exhaustive_photo_species_audit.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": sha256_file(input_path),
        "n_input_photo_rows": int(len(metadata)),
        "n_queue_photos": int(len(queue)),
        "n_queue_observations": int(queue["obs_id"].nunique()),
        "n_queue_species": int(queue["taxon_name"].nunique()),
        "parameters": vars(args),
        "selection_boundary": {
            "before_prediction": "taxonomic quality, licence, captivity and optional coordinate filters only",
            "not_applied": "no species cap, no spatial thinning, no environmental balancing, no trait filtering",
            "after_prediction": "all predictions are archived; weighting/thinning may create derived analysis cohorts only",
        },
    }
    (out_dir / "exhaustive_photo_queue_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
