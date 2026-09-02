#!/usr/bin/env python3
"""Create disclosure-safe Figure 2 inputs from the frozen v2 coordinate cohort.

The observation-level source remains an immutable Actions artifact. This helper
exports only two-degree cell summaries, so the repository can reproduce the map
without publishing observation identifiers or point coordinates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reproducibility" / "figures" / "source"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-cohort", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cell-width-degrees", type=float, default=2.0)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    args = parse_args()
    frame = pd.read_csv(args.source_cohort)
    required = {
        "obs_id",
        "taxon_name",
        "latitude",
        "longitude",
        "chelsa_bio01",
        "chelsa_bio04",
        "chelsa_bio12",
        "chelsa_bio15",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise SystemExit(f"source cohort lacks required columns: {missing}")
    if len(frame) != 46276 or frame["taxon_name"].nunique() != 259:
        raise SystemExit("source cohort is not the frozen 46,276-observation / 259-taxon v2 cohort")
    if frame["obs_id"].duplicated().any():
        raise SystemExit("frozen spatial cohort must contain one row per observation")

    width = float(args.cell_width_degrees)
    if width <= 0 or 180 % width or 360 % width:
        raise SystemExit("cell width must divide both 180 and 360 degrees")
    frame = frame.copy()
    frame["cell_lon_min"] = np.floor((frame["longitude"] + 180.0) / width) * width - 180.0
    frame["cell_lat_min"] = np.floor((frame["latitude"] + 90.0) / width) * width - 90.0
    frame["bio1_degrees_c"] = pd.to_numeric(frame["chelsa_bio01"], errors="coerce") * 0.1 - 273.15
    frame["bio12_mm"] = pd.to_numeric(frame["chelsa_bio12"], errors="coerce")
    grouped = (
        frame.groupby(["cell_lon_min", "cell_lat_min"], as_index=False)
        .agg(
            n_observations=("obs_id", "size"),
            n_taxa=("taxon_name", "nunique"),
            bio1_degrees_c_mean=("bio1_degrees_c", "mean"),
            bio12_mm_mean=("bio12_mm", "mean"),
        )
        .sort_values(["cell_lat_min", "cell_lon_min"])
        .reset_index(drop=True)
    )
    grouped["cell_lon_center"] = grouped["cell_lon_min"] + width / 2.0
    grouped["cell_lat_center"] = grouped["cell_lat_min"] + width / 2.0
    grouped = grouped[
        [
            "cell_lon_min",
            "cell_lat_min",
            "cell_lon_center",
            "cell_lat_center",
            "n_observations",
            "n_taxa",
            "bio1_degrees_c_mean",
            "bio12_mm_mean",
        ]
    ]
    if int(grouped["n_observations"].sum()) != 46276:
        raise SystemExit("aggregated cell counts do not recover the frozen cohort size")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "Figure_2_v2_realized_sampling_cells.csv"
    grouped.to_csv(output, index=False, lineterminator="\n")
    report = {
        "status": "PASS",
        "source_lane": "full27_full_environment_only",
        "source_run_id": 32975451732,
        "source_artifact_id": 9612943217,
        "source_artifact_digest": "sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e",
        "source_member": "environment/strict_spatial_chelsa.csv",
        "source_member_sha256": sha256(args.source_cohort),
        "aggregation": f"{width:g}-degree latitude-longitude cells",
        "n_observations": 46276,
        "n_taxa": 259,
        "n_cells": len(grouped),
        "direct_identifiers_in_output": False,
        "point_coordinates_in_output": False,
        "role": "current-v2 realized sampling-domain presentation only",
    }
    (args.output_dir / "Figure_2_v2_realized_sampling_cells_provenance.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
