#!/usr/bin/env python3
"""Prepare observation-level input for the Chapter 1 spatial robustness audit.

This utility does not fit or alter any model. It joins a frozen observation table,
a frozen fitted-value/residual export, and a manually reviewed broad-region table.
The result is validated before it can be passed to `audit_spatial_robustness.py`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from azami_ch1.provenance import write_json
from azami_ch1.tabular import assert_unique, require_columns, require_complete_text

OUTPUT_COLUMNS = [
    "observation_id",
    "taxon_name",
    "endpoint",
    "latitude",
    "longitude",
    "observed",
    "fitted",
    "residual",
    "broad_region",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--observation-id", default="observation_id")
    parser.add_argument("--endpoint", default="endpoint")
    parser.add_argument("--taxon", default="taxon_name")
    parser.add_argument("--latitude", default="latitude")
    parser.add_argument("--longitude", default="longitude")
    parser.add_argument("--observed", default="observed")
    parser.add_argument("--fitted", default="fitted")
    parser.add_argument("--residual", default="residual")
    args = parser.parse_args()

    observations = pd.read_csv(args.observations)
    predictions = pd.read_csv(args.predictions)
    regions = pd.read_csv(args.regions)

    require_columns(
        observations,
        [args.observation_id, args.taxon, args.latitude, args.longitude],
        "observation table",
    )
    require_columns(predictions, [args.observation_id, args.endpoint], "prediction table")
    require_columns(regions, [args.observation_id, "broad_region"], "region table")

    if args.observed not in predictions.columns and args.observed not in observations.columns:
        raise ValueError("Observed values are absent from both input tables")
    if args.fitted not in predictions.columns and args.residual not in predictions.columns:
        raise ValueError("Prediction table must contain fitted values or residuals")

    assert_unique(observations, [args.observation_id], "observation table")
    assert_unique(predictions, [args.observation_id, args.endpoint], "prediction table")
    assert_unique(regions, [args.observation_id], "region table")
    require_complete_text(regions, "broad_region", "reviewed broad_region")

    obs_keep = [args.observation_id, args.taxon, args.latitude, args.longitude]
    if args.observed in observations.columns:
        obs_keep.append(args.observed)
    obs = observations[obs_keep].rename(
        columns={
            args.taxon: "taxon_name",
            args.latitude: "latitude",
            args.longitude: "longitude",
            args.observed: "observed_from_observations",
        }
    )

    pred = predictions.rename(
        columns={
            args.endpoint: "endpoint",
            args.observed: "observed_from_predictions",
            args.fitted: "fitted",
            args.residual: "residual",
        }
    )

    merged = pred.merge(obs, on=args.observation_id, how="left", validate="many_to_one")
    merged = merged.merge(
        regions[[args.observation_id, "broad_region"]],
        on=args.observation_id,
        how="left",
        validate="many_to_one",
    )

    if merged["taxon_name"].isna().any():
        raise ValueError("At least one prediction row did not match an observation")
    require_complete_text(merged, "broad_region", "reviewed broad_region")

    if "observed_from_predictions" in merged.columns:
        merged["observed"] = merged["observed_from_predictions"]
        if "observed_from_observations" in merged.columns:
            left = pd.to_numeric(merged["observed_from_predictions"], errors="coerce")
            right = pd.to_numeric(merged["observed_from_observations"], errors="coerce")
            if not left.equals(right):
                raise ValueError("Observed values disagree between observation and prediction tables")
    else:
        merged["observed"] = merged["observed_from_observations"]

    if "residual" not in merged.columns:
        merged["residual"] = merged["observed"] - merged["fitted"]
    elif "fitted" not in merged.columns:
        merged["fitted"] = merged["observed"] - merged["residual"]

    merged = merged.rename(columns={args.observation_id: "observation_id"})
    output = merged[OUTPUT_COLUMNS].copy()
    numeric_columns = ["latitude", "longitude", "observed", "fitted", "residual"]
    for column in numeric_columns:
        output[column] = pd.to_numeric(output[column], errors="coerce")
    if output[numeric_columns].isna().any().any():
        counts = output[numeric_columns].isna().sum().to_dict()
        raise ValueError(f"Missing or non-numeric diagnostic values: {counts}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    metadata = {
        "observations": str(args.observations),
        "predictions": str(args.predictions),
        "regions": str(args.regions),
        "output": str(args.output),
        "n_rows": int(len(output)),
        "n_observations": int(output["observation_id"].nunique()),
        "n_taxa": int(output["taxon_name"].nunique()),
        "n_endpoints": int(output["endpoint"].nunique()),
        "regions_present": sorted(output["broad_region"].astype(str).unique().tolist()),
    }
    metadata_path = write_json(
        args.output.with_suffix(".metadata.json"),
        metadata,
        include_generated_utc=True,
    )
    print(metadata_path.read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
