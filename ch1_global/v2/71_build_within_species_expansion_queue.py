#!/usr/bin/env python3
"""Build a high-density queue for within-species trait–climate analysis.

This queue is separate from the 40-observation species-balanced atlas queue.
It keeps one primary photo per observation, requires broad spatial coverage, and
allows substantially more observations per eligible species. Selection never
uses measured capitulum traits.
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
    "obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank",
    "quality_grade", "captive", "latitude", "longitude",
    "coordinate_usable_for_environment", "positional_accuracy", "geoprivacy",
    "obscured", "observed_on", "small_image_url", "medium_image_url",
    "large_image_url", "source_file_stem",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--quality-grades", default="research")
    parser.add_argument("--max-positional-accuracy-m", type=float, default=10_000)
    parser.add_argument("--min-observations-per-species", type=int, default=60)
    parser.add_argument("--min-spatial-blocks-per-species", type=int, default=8)
    parser.add_argument("--max-observations-per-species", type=int, default=250)
    parser.add_argument("--max-per-species-spatial-block", type=int, default=15)
    parser.add_argument("--spatial-block-deg", type=float, default=2.0)
    parser.add_argument("--target-species", default="", help="Optional comma-separated taxa")
    parser.add_argument("--seed", type=int, default=20260711)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def finite_number(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def spatial_block(latitude: Any, longitude: Any, width: float) -> str:
    lat, lon = finite_number(latitude), finite_number(longitude)
    if not math.isfinite(lat) or not math.isfinite(lon):
        return "missing"
    return f"lat{math.floor((lat + 90.0) / width):03d}_lon{math.floor((lon + 180.0) / width):03d}"


def choose_spread(group: pd.DataFrame, maximum: int, per_block: int, rng: np.random.Generator) -> pd.DataFrame:
    """Round-robin spatial blocks before repeat observations within blocks."""
    work = group.copy()
    work["_random"] = rng.random(len(work))
    work = work.sort_values(["spatial_block", "_random", "photo_id"])
    work["_rank_in_block"] = work.groupby("spatial_block", sort=True).cumcount()
    work = work.loc[work["_rank_in_block"] < per_block].copy()
    work = work.sort_values(["_rank_in_block", "_random", "spatial_block", "photo_id"])
    return work.head(maximum).drop(columns=["_random", "_rank_in_block"])


def build_queue(metadata: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = sorted(REQUIRED.difference(metadata.columns))
    if missing:
        raise ValueError(f"Merged metadata missing required columns: {missing}")

    quality_grades = {x.strip().lower() for x in args.quality_grades.split(",") if x.strip()}
    targets = {x.strip() for x in args.target_species.split(",") if x.strip()}
    work = metadata.copy()
    work["taxon_rank"] = work["taxon_rank"].map(text).str.lower()
    work["quality_grade"] = work["quality_grade"].map(text).str.lower()
    work["captive_bool"] = work["captive"].map(as_bool)
    work["coordinate_usable_bool"] = work["coordinate_usable_for_environment"].map(as_bool)
    work["positional_accuracy_numeric"] = pd.to_numeric(work["positional_accuracy"], errors="coerce")
    work["spatial_block"] = [
        spatial_block(lat, lon, args.spatial_block_deg)
        for lat, lon in zip(work["latitude"], work["longitude"])
    ]

    eligible = work.loc[
        work["obs_id"].map(text).ne("")
        & work["photo_id"].map(text).ne("")
        & work["taxon_name"].map(text).ne("")
        & work["taxon_rank"].eq("species")
        & work["quality_grade"].isin(quality_grades)
        & ~work["captive_bool"]
        & work["coordinate_usable_bool"]
        & work["spatial_block"].ne("missing")
        & work["positional_accuracy_numeric"].notna()
        & work["positional_accuracy_numeric"].le(args.max_positional_accuracy_m)
    ].copy()
    if targets:
        eligible = eligible.loc[eligible["taxon_name"].isin(targets)].copy()
    if eligible.empty:
        raise ValueError("No eligible observations after filtering")

    eligible["photo_index_numeric"] = pd.to_numeric(eligible["photo_index"], errors="coerce").fillna(999999)
    eligible = (
        eligible.sort_values(["obs_id", "photo_index_numeric", "photo_id"])
        .groupby("obs_id", as_index=False, sort=False)
        .head(1)
        .copy()
    )

    audit = (
        eligible.groupby("taxon_name", sort=True)
        .agg(
            n_eligible_observations=("obs_id", "size"),
            n_spatial_blocks=("spatial_block", "nunique"),
            latitude_min=("latitude", lambda x: pd.to_numeric(x, errors="coerce").min()),
            latitude_max=("latitude", lambda x: pd.to_numeric(x, errors="coerce").max()),
            longitude_min=("longitude", lambda x: pd.to_numeric(x, errors="coerce").min()),
            longitude_max=("longitude", lambda x: pd.to_numeric(x, errors="coerce").max()),
        )
        .reset_index()
    )
    audit["eligible_for_expansion"] = (
        audit["n_eligible_observations"].ge(args.min_observations_per_species)
        & audit["n_spatial_blocks"].ge(args.min_spatial_blocks_per_species)
    )
    accepted = set(audit.loc[audit["eligible_for_expansion"], "taxon_name"])
    if not accepted:
        raise ValueError("No species meet the high-density expansion criteria")

    rng = np.random.default_rng(args.seed)
    parts = [
        choose_spread(group, args.max_observations_per_species, args.max_per_species_spatial_block, rng)
        for _, group in eligible.loc[eligible["taxon_name"].isin(accepted)].groupby("taxon_name", sort=True)
    ]
    queue = pd.concat(parts, ignore_index=True)
    queue = queue.sort_values(["taxon_name", "spatial_block", "photo_id"]).reset_index(drop=True)
    queue.insert(0, "queue_id", [f"ch1within_{i:08d}" for i in range(1, len(queue) + 1)])
    queue["queue_status"] = "not_downloaded"
    queue["screen_download_filename"] = queue["source_file_stem"].map(lambda x: f"{text(x)}.jpg")
    queue["sampling_role"] = "within_species_high_density"

    selected = queue.groupby("taxon_name").agg(
        n_selected_observations=("queue_id", "size"),
        n_selected_spatial_blocks=("spatial_block", "nunique"),
    ).reset_index()
    audit = audit.merge(selected, on="taxon_name", how="left")
    audit[["n_selected_observations", "n_selected_spatial_blocks"]] = audit[
        ["n_selected_observations", "n_selected_spatial_blocks"]
    ].fillna(0).astype(int)
    return queue, audit


def main() -> None:
    args = parse_args()
    if args.max_positional_accuracy_m <= 0 or args.spatial_block_deg <= 0:
        raise ValueError("Coordinate and spatial-block thresholds must be positive")
    if args.min_observations_per_species < 2 or args.min_spatial_blocks_per_species < 2:
        raise ValueError("Within-species expansion requires at least two observations and blocks")
    if args.max_observations_per_species < args.min_observations_per_species:
        raise ValueError("Maximum observations must be at least the minimum")
    if args.max_per_species_spatial_block < 1:
        raise ValueError("Per-block maximum must be positive")

    input_path = Path(args.input)
    metadata = pd.read_csv(input_path, dtype=str, keep_default_na=False, low_memory=False)
    queue, audit = build_queue(metadata, args)

    keep = [
        "queue_id", "queue_status", "sampling_role", "obs_id", "photo_id", "photo_index",
        "taxon_name", "taxon_rank", "quality_grade", "spatial_block", "latitude", "longitude",
        "coordinate_usable_for_environment", "positional_accuracy", "observed_on", "captive",
        "geoprivacy", "obscured", "small_image_url", "medium_image_url", "large_image_url",
        "source_file_stem", "screen_download_filename",
    ]
    for column in keep:
        if column not in queue.columns:
            queue[column] = ""

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue[keep].to_csv(out_dir / "within_species_expansion_queue.csv", index=False, encoding="utf-8-sig")
    audit.to_csv(out_dir / "within_species_expansion_species_audit.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": sha256_file(input_path),
        "n_input_rows": int(len(metadata)),
        "n_queue_observations": int(len(queue)),
        "n_expansion_species": int(queue["taxon_name"].nunique()),
        "parameters": vars(args),
        "design_boundary": {
            "between_species_dataset": "unchanged 40-observation balanced atlas",
            "within_species_dataset": "separate high-density queue with strict coordinate quality",
            "trait_leakage": "selection uses metadata and spatial coverage only, never measured traits",
            "inference": "expanded data evaluate detectability; they do not make associations causal",
        },
    }
    (out_dir / "within_species_expansion_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
