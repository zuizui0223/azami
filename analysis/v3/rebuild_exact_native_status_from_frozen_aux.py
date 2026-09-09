#!/usr/bin/env python3
"""Regenerate the v2 observation native-status table without its missing LFS object.

The original observation-level CSV is still referenced by the immutable v2 tag,
but its Git LFS object is no longer present on the server. The same frozen
analysis directory retains the exact non-LFS name-resolution and distribution
records used by the v2 classification. This script combines those frozen
auxiliary tables with the frozen 46,276-row strict-spatial observation table and
the pinned TDWG level-3 geometry.

The regenerated table is deliberately NOT called byte-identical to the missing
historical CSV. Instead, it must reproduce every frozen structural invariant and
is given its own deterministic SHA-256 fingerprint for downstream reuse.
"""
from __future__ import annotations

import argparse
import hashlib
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

EXPECTED_REGENERATED_SHA256 = "9686b8f515deef3b3aa9311b5137317a094af195e0f2414f2b1b70d7c72b5021"
EXPECTED_RESOLUTION_SHA256 = "89f95030e3c692845a1c32cc3de261e50fa45b5900925be9817106a56b2b7da9"
EXPECTED_DISTRIBUTIONS_SHA256 = "a85436363cd107adee7d254af59faacfe484b74bfcfc2f281b8ba9489ca935e2"


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


def restore_original_resolution_dtypes(path: Path) -> pd.DataFrame:
    resolution = pd.read_csv(
        path,
        dtype={
            "input_name": "string",
            "resolution_status": "string",
            "accepted_key": "string",
            "accepted_name": "string",
        },
        keep_default_na=False,
        low_memory=False,
    )
    resolution["accepted_key"] = resolution["accepted_key"].map(
        lambda value: int(value) if str(value).strip() else ""
    )
    return resolution


def native_id_fingerprint(classified: pd.DataFrame) -> str:
    native_ids = sorted(
        classified.loc[
            classified["native_range_status"].eq("native"), "obs_id"
        ].astype(str)
    )
    return hashlib.sha256(("\n".join(native_ids) + "\n").encode()).hexdigest()


def main() -> int:
    args = parse_args()
    observation_sha = sha256_file(args.observation)
    resolution_sha = sha256_file(args.resolution)
    distributions_sha = sha256_file(args.distributions)
    if observation_sha != EXPECTED_OBSERVATION_SHA256:
        raise SystemExit(f"strict-spatial observation SHA changed: {observation_sha}")
    if resolution_sha != EXPECTED_RESOLUTION_SHA256:
        raise SystemExit(f"frozen name-resolution SHA changed: {resolution_sha}")
    if distributions_sha != EXPECTED_DISTRIBUTIONS_SHA256:
        raise SystemExit(f"frozen distribution-record SHA changed: {distributions_sha}")

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

    resolution = restore_original_resolution_dtypes(args.resolution)
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
        raise SystemExit(f"Native-status counts differ from frozen v2: {status_counts}")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    classified[OUTPUT_COLUMNS].to_csv(args.out_csv, index=False)
    output_sha = sha256_file(args.out_csv)
    if output_sha != EXPECTED_REGENERATED_SHA256:
        raise SystemExit(
            "Deterministic regenerated native-status fingerprint changed: "
            f"{output_sha}; expected {EXPECTED_REGENERATED_SHA256}"
        )

    old_byte_match = output_sha == EXPECTED_OUTPUT_SHA256
    report = {
        "status": "PASS_REGENERATED_MEMBERSHIP_BYTE_IDENTITY_UNAVAILABLE",
        "method": "regenerate_from_frozen_non_lfs_auxiliary_tables",
        "source_observation_sha256": observation_sha,
        "resolution_sha256": resolution_sha,
        "distributions_sha256": distributions_sha,
        "tdwg_geojson_sha256": tdwg_sha,
        "n_rows": int(len(classified)),
        "n_taxa": int(classified["taxon_name"].nunique()),
        "n_resolved_taxa": resolved_taxa,
        "native_status_counts": status_counts,
        "native_obs_id_fingerprint_sha256": native_id_fingerprint(classified),
        "output_size_bytes": int(args.out_csv.stat().st_size),
        "regenerated_output_sha256": output_sha,
        "missing_historical_lfs_output_sha256": EXPECTED_OUTPUT_SHA256,
        "historical_byte_identity_available": old_byte_match,
        "provenance_note": (
            "The immutable tag retains the original Git LFS pointer, but GitHub's "
            "LFS batch API reports that object as absent. The regenerated cohort "
            "therefore uses the frozen source rows, frozen WCVP name-resolution "
            "table, frozen WCVP distribution table, and pinned TDWG geometry."
        ),
        "claim_boundary": (
            "Use this as a deterministic regenerated frozen-v2 native cohort. "
            "Do not describe it as byte-identical recovery of the missing LFS CSV."
        ),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
