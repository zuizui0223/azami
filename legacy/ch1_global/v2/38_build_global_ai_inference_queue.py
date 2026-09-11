#!/usr/bin/env python3
"""Build the automated global Cirsium inference queue from merged metadata.

This queue is observation-deduplicated before image download, uses no trait
labels or iNaturalist flower annotation for selection, and balances samples
within species across coarse spatial blocks. It is intended for fully automated
detector + ensemble-trait inference, not human annotation.
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


REQUIRED = {
    "obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank", "quality_grade", "captive",
    "latitude", "longitude", "coordinate_usable_for_environment", "positional_accuracy", "geoprivacy", "obscured",
    "observed_on", "small_image_url", "medium_image_url", "large_image_url", "source_file_stem",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a research-grade, spatially balanced automated AI inference queue.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--quality-grades", default="research", help="Comma-separated permitted iNaturalist quality grades")
    parser.add_argument("--min-observations-per-species", type=int, default=3)
    parser.add_argument("--max-observations-per-species", type=int, default=40)
    parser.add_argument("--max-per-species-spatial-block", type=int, default=3)
    parser.add_argument("--spatial-block-deg", type=float, default=5.0)
    parser.add_argument("--seed", type=int, default=20260705)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def spatial_block(latitude: Any, longitude: Any, width: float) -> str:
    try:
        lat, lon = float(latitude), float(longitude)
    except (TypeError, ValueError):
        return "missing"
    if not math.isfinite(lat) or not math.isfinite(lon):
        return "missing"
    return f"lat{math.floor((lat + 90.0) / width):03d}_lon{math.floor((lon + 180.0) / width):03d}"


def choose_species_spread(group: pd.DataFrame, maximum: int, per_block: int, rng: np.random.Generator) -> pd.DataFrame:
    """Prioritise one photo per spatial block before repeat sampling within blocks."""
    work = group.copy()
    work["_random"] = rng.random(len(work))
    work = work.sort_values(["spatial_block", "_random", "photo_id"]).copy()
    work["_rank_in_block"] = work.groupby("spatial_block", sort=True).cumcount()
    work = work.loc[work["_rank_in_block"] < per_block].copy()
    work = work.sort_values(["_rank_in_block", "_random", "photo_id"]).head(maximum).copy()
    return work.drop(columns=["_random", "_rank_in_block"])


def main() -> None:
    args = parse_args()
    if args.min_observations_per_species < 1 or args.max_observations_per_species < 1:
        raise ValueError("Species sample limits must be positive")
    if args.min_observations_per_species > args.max_observations_per_species:
        raise ValueError("Minimum cannot exceed maximum observations per species")
    if args.max_per_species_spatial_block < 1 or args.spatial_block_deg <= 0:
        raise ValueError("Spatial block parameters must be positive")
    quality_grades = {grade.strip().lower() for grade in args.quality_grades.split(",") if grade.strip()}
    if not quality_grades:
        raise ValueError("At least one quality grade is required")

    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    metadata = pd.read_csv(input_path, dtype=str, keep_default_na=False, low_memory=False)
    missing = sorted(REQUIRED.difference(metadata.columns))
    if missing:
        raise ValueError(f"Merged metadata missing required columns: {missing}")

    metadata["taxon_rank"] = metadata["taxon_rank"].map(text).str.lower()
    metadata["quality_grade"] = metadata["quality_grade"].map(text).str.lower()
    metadata["captive_bool"] = metadata["captive"].map(as_bool)
    metadata["coordinate_usable_bool"] = metadata["coordinate_usable_for_environment"].map(as_bool)
    metadata["spatial_block"] = [spatial_block(lat, lon, args.spatial_block_deg) for lat, lon in zip(metadata["latitude"], metadata["longitude"])]

    eligible = metadata.loc[
        metadata["obs_id"].map(text).ne("")
        & metadata["photo_id"].map(text).ne("")
        & metadata["taxon_name"].map(text).ne("")
        & metadata["taxon_rank"].eq("species")
        & metadata["quality_grade"].isin(quality_grades)
        & ~metadata["captive_bool"]
        & metadata["coordinate_usable_bool"]
        & metadata["spatial_block"].ne("missing")
    ].copy()
    if eligible.empty:
        raise ValueError("No eligible metadata rows after research-grade/taxon/coordinate filtering")

    eligible["photo_index_numeric"] = pd.to_numeric(eligible["photo_index"], errors="coerce").fillna(999999)
    eligible = (
        eligible.sort_values(["obs_id", "photo_index_numeric", "photo_id"])
        .groupby("obs_id", as_index=False, sort=False)
        .head(1)
        .copy()
    )
    species_counts = eligible["taxon_name"].value_counts()
    eligible = eligible.loc[eligible["taxon_name"].map(species_counts).ge(args.min_observations_per_species)].copy()
    if eligible.empty:
        raise ValueError("No species meet the minimum observation threshold")

    rng = np.random.default_rng(args.seed)
    selected_parts = [
        choose_species_spread(group, args.max_observations_per_species, args.max_per_species_spatial_block, rng)
        for _, group in eligible.groupby("taxon_name", sort=True)
    ]
    queue = pd.concat(selected_parts, ignore_index=True)
    queue = queue.sort_values(["taxon_name", "spatial_block", "photo_id"]).reset_index(drop=True)
    queue.insert(0, "queue_id", [f"ch1atlas_{index:07d}" for index in range(1, len(queue) + 1)])
    queue["screen_download_filename"] = queue["source_file_stem"].map(lambda value: f"{text(value)}.jpg")
    queue["queue_status"] = "not_downloaded"

    keep = [
        "queue_id", "queue_status", "obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank", "quality_grade",
        "spatial_block", "latitude", "longitude", "coordinate_usable_for_environment", "positional_accuracy", "observed_on",
        "captive", "geoprivacy", "obscured", "inat_flowers_annotation", "inat_annotation_summary",
        "small_image_url", "medium_image_url", "large_image_url", "source_file_stem", "screen_download_filename",
    ]
    for column in keep:
        if column not in queue.columns:
            queue[column] = ""

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue[keep].to_csv(out_dir / "global_ai_inference_queue.csv", index=False, encoding="utf-8-sig")
    (
        queue.groupby("taxon_name", sort=True)
        .agg(n_queue_observations=("queue_id", "size"), n_spatial_blocks=("spatial_block", "nunique"))
        .reset_index()
        .sort_values(["n_queue_observations", "taxon_name"], ascending=[False, True])
        .to_csv(out_dir / "global_ai_queue_species_qc.csv", index=False, encoding="utf-8-sig")
    )
    (
        queue.groupby("spatial_block", sort=True)
        .agg(n_queue_observations=("queue_id", "size"), n_species=("taxon_name", "nunique"))
        .reset_index()
        .sort_values(["n_queue_observations", "spatial_block"], ascending=[False, True])
        .to_csv(out_dir / "global_ai_queue_spatial_qc.csv", index=False, encoding="utf-8-sig")
    )
    provenance = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": sha256_file(input_path),
        "rows_input": int(len(metadata)),
        "rows_eligible_after_observation_deduplication": int(len(eligible)),
        "queue_observations": int(len(queue)),
        "species_in_queue": int(queue["taxon_name"].nunique()),
        "parameters": vars(args),
        "selection_rules": {
            "taxonomic_quality": "species rank, research-grade observation, non-captive",
            "geography": "public coordinate usable for environmental extraction; spatially balanced within species",
            "observation_unit": "one primary image per iNaturalist observation before sampling",
            "trait_leakage": "No trait label, AI trait result, or iNaturalist flower annotation is used for queue selection",
        },
    }
    (out_dir / "global_ai_queue_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "queue_observations": provenance["queue_observations"],
        "species_in_queue": provenance["species_in_queue"],
        "rows_eligible_after_observation_deduplication": provenance["rows_eligible_after_observation_deduplication"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
