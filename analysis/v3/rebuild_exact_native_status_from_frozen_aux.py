#!/usr/bin/env python3
"""Rebuild the frozen v2 observation native-status table without Git LFS.

The original observation-level CSV is still referenced by the immutable v2 tag,
but its Git LFS object is no longer present on the server.  The same frozen
analysis directory retains the exact non-LFS name-resolution and distribution
records that generated it.  This script combines those frozen auxiliary tables
with the frozen 46,276-row strict-spatial observation table and the pinned TDWG
level-3 geometry, then requires the reconstructed CSV to match the original
recorded SHA-256 byte-for-byte.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from analysis.rebuild_frozen_native_status import (
    EXPECTED_OBSERVATION_SHA256,
    EXPECTED_OUTPUT_SHA256,
    EXPECTED_RESOLVED_TAXA,
    EXPECTED_ROWS,
    EXPECTED_STATUS_COUNTS,
    EXPECTED_TAXA,
    EXPECTED_TDWG_GEOJSON_SHA256,
    OUTPUT_COLUMNS,
    classify_observations,
    load_pinned_level3,
    sha256_file,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", type=Path, required=True)
    parser.add_argument("--resolution", type=Path, required=True)
    parser.add_argument("--distributions", type=Path, required=True)
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("analysis/ch1/native_range_sensitivity_contract.json"),
    )
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    parser.add_argument("--max-retries", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    observation_sha = sha256_file(args.observation)
    if observation_sha != EXPECTED_OBSERVATION_SHA256:
        raise SystemExit(
            f"strict-spatial observation SHA changed: {observation_sha}"
        )

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if contract.get("status") != "locked_before_native_range_outcome_execution":
        raise SystemExit("Native-range contract is not the frozen locked contract")

    source = pd.read_csv(
        args.observation,
        usecols=["obs_id", "taxon_name", "latitude", "longitude"],
        low_memory=False,
    )
    if len(source) != EXPECTED_ROWS or source["obs_id"].nunique() != EXPECTED_ROWS:
        raise SystemExit("Strict-spatial source row identity/count is not frozen")
    if source["taxon_name"].nunique() != EXPECTED_TAXA:
        raise SystemExit("Strict-spatial source taxon count is not frozen")

    resolution = pd.read_csv(args.resolution, low_memory=False)
    distributions = pd.read_csv(args.distributions, low_memory=False)
    expected_names = set(source["taxon_name"].astype(str))
    observed_names = set(resolution["input_name"].astype(str))
    if expected_names != observed_names or len(resolution) != EXPECTED_TAXA:
        raise SystemExit("Frozen name-resolution table does not cover the frozen taxon set")
    resolved_taxa = int(
        resolution["resolution_status"].eq("resolved_unique_accepted_key").sum()
    )
    if resolved_taxa != EXPECTED_RESOLVED_TAXA:
        raise SystemExit(
            f"Expected {EXPECTED_RESOLVED_TAXA} resolved taxa, found {resolved_taxa}"
        )

    geojson, tdwg_sha = load_pinned_level3(
        contract, args.timeout_seconds, args.max_retries
    )
    if tdwg_sha != EXPECTED_TDWG_GEOJSON_SHA256:
        raise SystemExit(
            f"Pinned TDWG geometry SHA changed: {tdwg_sha}; "
            f"expected {EXPECTED_TDWG_GEOJSON_SHA256}"
        )

    classified = classify_observations(source, resolution, distributions, geojson)
    status_counts = {
        str(key): int(value)
        for key, value in classified["native_range_status"].value_counts().items()
    }
    if status_counts != EXPECTED_STATUS_COUNTS:
        raise SystemExit(
            f"Native-status counts differ from frozen v2: {status_counts}"
        )

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    classified[OUTPUT_COLUMNS].to_csv(args.out_csv, index=False)
    output_sha = sha256_file(args.out_csv)
    exact = output_sha == EXPECTED_OUTPUT_SHA256

    report = {
        "status": "PASS" if exact else "FAIL",
        "method": "rebuild_from_frozen_non_lfs_auxiliary_tables",
        "source_observation_sha256": observation_sha,
        "resolution_sha256": sha256_file(args.resolution),
        "distributions_sha256": sha256_file(args.distributions),
        "tdwg_geojson_sha256": tdwg_sha,
        "n_rows": int(len(classified)),
        "n_taxa": int(classified["taxon_name"].nunique()),
        "n_resolved_taxa": resolved_taxa,
        "native_status_counts": status_counts,
        "output_sha256": output_sha,
        "expected_output_sha256": EXPECTED_OUTPUT_SHA256,
        "exact_byte_match": exact,
        "claim_boundary": (
            "This only reconstructs the frozen v2 native-status membership from "
            "its frozen auxiliary tables; it adds no new taxonomic or ecological inference."
        ),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not exact:
        raise SystemExit("Frozen native-status byte identity was not recovered")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
