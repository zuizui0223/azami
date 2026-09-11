#!/usr/bin/env python3
"""Create derived analysis cohorts after exhaustive all-photo prediction.

The source observation table is never overwritten. This script builds explicit
coordinate-quality and spatially thinned cohorts for sensitivity analyses, plus
a separate species-balanced atlas cohort. Selection never uses trait values.
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
    "obs_id", "taxon_name", "latitude", "longitude", "positional_accuracy",
    "coordinate_usable_for_environment",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--strict-positional-accuracy-m", type=float, default=10_000)
    parser.add_argument("--cell-deg", type=float, default=0.25)
    parser.add_argument("--between-species-max-per-taxon", type=int, default=40)
    parser.add_argument("--seed", default="20260711")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def stable_hash(seed: str, *parts: Any) -> str:
    payload = "|".join([seed, *[text(part) for part in parts]])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def spatial_cell(lat: Any, lon: Any, width: float) -> str:
    try:
        latitude, longitude = float(lat), float(lon)
    except (TypeError, ValueError):
        return "missing"
    if not math.isfinite(latitude) or not math.isfinite(longitude):
        return "missing"
    return f"lat{math.floor((latitude + 90) / width):04d}_lon{math.floor((longitude + 180) / width):04d}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_cohorts(frame: pd.DataFrame, args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    missing = sorted(REQUIRED.difference(frame.columns))
    if missing:
        raise ValueError(f"Observation table missing required columns: {missing}")
    if frame.empty or frame["obs_id"].astype(str).duplicated().any():
        raise ValueError("Observation table must contain unique detector-positive observations")

    work = frame.copy()
    work["coordinate_usable_bool"] = work["coordinate_usable_for_environment"].map(as_bool)
    work["positional_accuracy_numeric"] = pd.to_numeric(work["positional_accuracy"], errors="coerce")
    work["analysis_cell"] = [
        spatial_cell(lat, lon, args.cell_deg)
        for lat, lon in zip(work["latitude"], work["longitude"])
    ]
    work["selection_hash"] = [
        stable_hash(args.seed, taxon, cell, obs)
        for taxon, cell, obs in zip(work["taxon_name"], work["analysis_cell"], work["obs_id"])
    ]

    all_detected = work.copy()
    coordinate = work.loc[work["coordinate_usable_bool"] & work["analysis_cell"].ne("missing")].copy()
    strict = coordinate.loc[
        coordinate["positional_accuracy_numeric"].notna()
        & coordinate["positional_accuracy_numeric"].le(args.strict_positional_accuracy_m)
    ].copy()
    spatial = (
        strict.sort_values(["taxon_name", "analysis_cell", "selection_hash", "obs_id"])
        .groupby(["taxon_name", "analysis_cell"], as_index=False, sort=False)
        .head(1)
        .copy()
    )
    balanced = (
        spatial.sort_values(["taxon_name", "selection_hash", "analysis_cell", "obs_id"])
        .groupby("taxon_name", as_index=False, sort=False)
        .head(args.between_species_max_per_taxon)
        .copy()
    )

    for cohort_name, cohort in {
        "all_detected": all_detected,
        "coordinate_usable": coordinate,
        "strict_10km": strict,
        "strict_spatial_thin": spatial,
        "between_species_balanced": balanced,
    }.items():
        cohort["analysis_cohort"] = cohort_name
        cohort["source_is_exhaustive"] = True
    return {
        "all_detected": all_detected,
        "coordinate_usable": coordinate,
        "strict_10km": strict,
        "strict_spatial_thin": spatial,
        "between_species_balanced": balanced,
    }


def main() -> None:
    args = parse_args()
    if args.strict_positional_accuracy_m <= 0 or args.cell_deg <= 0:
        raise ValueError("Coordinate threshold and cell width must be positive")
    if args.between_species_max_per_taxon < 1:
        raise ValueError("between-species maximum must be positive")
    source = Path(args.observation)
    frame = pd.read_csv(source, low_memory=False)
    cohorts = build_cohorts(frame, args)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    filenames = {
        "all_detected": "all_detected_observations.csv",
        "coordinate_usable": "coordinate_usable_observations.csv",
        "strict_10km": "strict_10km_observations.csv",
        "strict_spatial_thin": "strict_spatial_thinned_observations.csv",
        "between_species_balanced": "between_species_balanced_40_observations.csv",
    }
    summary_rows = []
    for name, cohort in cohorts.items():
        cohort.drop(columns=["selection_hash"], errors="ignore").to_csv(
            out_dir / filenames[name], index=False, encoding="utf-8-sig"
        )
        summary_rows.append({
            "cohort": name,
            "filename": filenames[name],
            "n_observations": int(len(cohort)),
            "n_taxa": int(cohort["taxon_name"].nunique()),
            "n_cells": int(cohort["analysis_cell"].loc[cohort["analysis_cell"].ne("missing")].nunique()),
            "max_observations_per_taxon": int(cohort.groupby("taxon_name").size().max()) if not cohort.empty else 0,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "postprediction_cohort_summary.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": sha256_file(source),
        "parameters": vars(args),
        "cohorts": summary_rows,
        "selection_boundary": {
            "source": "exhaustive detector-positive observation table",
            "trait_blind": "cohort selection uses only taxon, coordinate quality, spatial cell and stable hash",
            "source_preservation": "all detected observations remain in all_detected_observations.csv; thinned cohorts are derived sensitivity tables",
        },
    }
    (out_dir / "postprediction_cohort_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
