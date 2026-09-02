#!/usr/bin/env python3
"""Build a species- and space-balanced image screening queue for Chapter 1.

Input is the metadata-only output of 04_collect_inat_cirsium_metadata.py.
This script does not download images and does not use a trait classifier. It
selects image candidates reproducibly before YOLO screening.

Example
-------
python ch1_global/v2/05_build_image_screening_queue.py \
  --input "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_raw\\photo_metadata.csv" \
  --out-dir "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_screen_queue" \
  --max-photos 5000
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


QUEUE_VERSION = "1.0.0"
REQUIRED = {
    "obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank", "latitude", "longitude",
    "coordinate_usable_for_environment", "captive", "small_image_url", "medium_image_url", "source_file_stem",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an image-screening queue without trait leakage.")
    parser.add_argument("--input", required=True, help="photo_metadata.csv produced by 04")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-photos", type=int, default=5000, help="0 means use all eligible rows after caps")
    parser.add_argument("--min-photos-per-species", type=int, default=10)
    parser.add_argument("--max-photos-per-species", type=int, default=150)
    parser.add_argument("--spatial-block-deg", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=20260703)
    parser.add_argument("--allow-non-species", action="store_true")
    parser.add_argument("--include-captive", action="store_true")
    parser.add_argument("--allow-non-open-coordinates", action="store_true")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return value is True or text(value).lower() in {"true", "1", "yes"}


def spatial_block(lat: Any, lon: Any, width: float) -> str:
    try:
        lat_f, lon_f = float(lat), float(lon)
    except (TypeError, ValueError):
        return "missing"
    if not math.isfinite(lat_f) or not math.isfinite(lon_f):
        return "missing"
    return f"lat{int(math.floor((lat_f + 90) / width)):03d}_lon{int(math.floor((lon_f + 180) / width)):03d}"


def choose_spread(df: pd.DataFrame, n: int, rng: np.random.Generator) -> pd.DataFrame:
    """Cycle across species after taking one candidate per species × block."""
    if n <= 0 or df.empty:
        return df.iloc[0:0].copy()
    work = df.copy()
    work["_random"] = rng.random(len(work))
    first = (
        work.sort_values(["taxon_name", "spatial_block", "_random", "photo_id"])
        .groupby(["taxon_name", "spatial_block"], sort=True, as_index=False)
        .head(1)
    )
    selected: list[int] = []
    per_species = {name: list(part.index) for name, part in first.groupby("taxon_name", sort=True)}
    while len(selected) < n and any(per_species.values()):
        for name in sorted(per_species):
            if len(selected) >= n:
                break
            if per_species[name]:
                selected.append(per_species[name].pop(0))
    remaining = work.drop(index=selected, errors="ignore")
    counts = work.loc[selected, "taxon_name"].value_counts().to_dict() if selected else {}
    while len(selected) < n and not remaining.empty:
        remaining = remaining.assign(_species_count=remaining["taxon_name"].map(counts).fillna(0))
        chosen = remaining.sort_values(["_species_count", "_random", "photo_id"]).index[0]
        name = remaining.at[chosen, "taxon_name"]
        selected.append(chosen)
        counts[name] = int(counts.get(name, 0)) + 1
        remaining = remaining.drop(index=chosen)
    return work.loc[selected].drop(columns=["_random"])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    if args.max_photos < 0 or args.min_photos_per_species < 1 or args.max_photos_per_species < 1:
        raise ValueError("Photo limits must be non-negative (max) or positive (min/cap)")
    if args.min_photos_per_species > args.max_photos_per_species:
        raise ValueError("--min-photos-per-species cannot exceed --max-photos-per-species")
    if args.spatial_block_deg <= 0:
        raise ValueError("--spatial-block-deg must be > 0")

    input_path = Path(args.input)
    df = pd.read_csv(input_path, low_memory=False, dtype={"obs_id": str, "photo_id": str})
    missing = REQUIRED.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required metadata columns: {sorted(missing)}")
    if df.empty:
        raise ValueError("Input metadata is empty")

    df["taxon_name"] = df["taxon_name"].map(text)
    df["taxon_rank"] = df["taxon_rank"].map(text).str.lower()
    df["captive_bool"] = df["captive"].map(as_bool)
    df["coordinate_usable_bool"] = df["coordinate_usable_for_environment"].map(as_bool)
    df["spatial_block"] = [spatial_block(lat, lon, args.spatial_block_deg) for lat, lon in zip(df["latitude"], df["longitude"])]

    eligible = df.loc[df["photo_id"].map(text).ne("") & df["taxon_name"].ne("")].copy()
    if not args.allow_non_species:
        eligible = eligible.loc[eligible["taxon_rank"].eq("species")].copy()
    if not args.include_captive:
        eligible = eligible.loc[~eligible["captive_bool"]].copy()
    if not args.allow_non_open_coordinates:
        eligible = eligible.loc[eligible["coordinate_usable_bool"] & eligible["spatial_block"].ne("missing")].copy()
    if eligible.empty:
        raise ValueError("No eligible rows after rank/captive/coordinate filters")

    # Never give serial images from the same observation extra weight at this stage.
    eligible["photo_index_numeric"] = pd.to_numeric(eligible["photo_index"], errors="coerce").fillna(999999)
    eligible = (
        eligible.sort_values(["obs_id", "photo_index_numeric", "photo_id"])
        .groupby("obs_id", sort=False, as_index=False)
        .head(1)
        .copy()
    )
    species_counts = eligible["taxon_name"].value_counts()
    eligible = eligible.loc[eligible["taxon_name"].map(species_counts).ge(args.min_photos_per_species)].copy()
    if eligible.empty:
        raise ValueError("No species meet --min-photos-per-species after observation deduplication")

    rng = np.random.default_rng(args.seed)
    capped_parts = [choose_spread(group, min(len(group), args.max_photos_per_species), rng) for _, group in eligible.groupby("taxon_name", sort=True)]
    capped = pd.concat(capped_parts, ignore_index=True)
    n_target = len(capped) if args.max_photos == 0 else min(args.max_photos, len(capped))
    queue = choose_spread(capped, n_target, rng).copy().reset_index(drop=True)
    queue.insert(0, "queue_id", [f"ch1q_{index:07d}" for index in range(1, len(queue) + 1)])
    queue["screen_download_filename"] = queue["source_file_stem"].map(lambda value: f"{text(value)}.jpg")
    queue["queue_status"] = "not_downloaded"

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    keep = [
        "queue_id", "queue_status", "obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank",
        "spatial_block", "latitude", "longitude", "coordinate_usable_for_environment", "positional_accuracy",
        "observed_on", "quality_grade", "captive", "geoprivacy", "obscured", "inat_flowers_annotation",
        "photo_license_code", "small_image_url", "medium_image_url", "large_image_url", "source_file_stem",
        "screen_download_filename",
    ]
    for column in keep:
        if column not in queue.columns:
            queue[column] = pd.NA
    queue[keep].to_csv(out_dir / "image_screening_queue.csv", index=False, encoding="utf-8-sig")

    species_qc = (
        queue.groupby("taxon_name")
        .agg(n_queue_images=("queue_id", "size"), n_spatial_blocks=("spatial_block", "nunique"))
        .reset_index()
        .sort_values(["n_queue_images", "taxon_name"], ascending=[False, True])
    )
    spatial_qc = (
        queue.groupby("spatial_block")
        .agg(n_queue_images=("queue_id", "size"), n_species=("taxon_name", "nunique"))
        .reset_index()
        .sort_values(["n_queue_images", "spatial_block"], ascending=[False, True])
    )
    species_qc.to_csv(out_dir / "queue_species_qc.csv", index=False, encoding="utf-8-sig")
    spatial_qc.to_csv(out_dir / "queue_spatial_qc.csv", index=False, encoding="utf-8-sig")
    run_info = {
        "queue_version": QUEUE_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path.resolve()),
        "input_sha256": sha256(input_path),
        "rows_input": int(len(df)),
        "rows_eligible_after_filters": int(len(eligible)),
        "queue_rows": int(len(queue)),
        "species_in_queue": int(queue["taxon_name"].nunique()),
        "parameters": vars(args),
        "note": "Flowering annotations were retained as metadata and were not used to choose candidates. Human flowering QC occurs later.",
    }
    (out_dir / "queue_provenance.json").write_text(json.dumps(run_info, ensure_ascii=False, indent=2), encoding="utf-8")
    print("[OK] wrote image screening queue", out_dir.resolve())
    print(json.dumps({key: run_info[key] for key in ["queue_rows", "species_in_queue", "rows_eligible_after_filters"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
