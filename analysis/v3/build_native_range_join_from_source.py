"""Attach WCVP/TDWG taxon and native-range status to source-only observations.

This runner is intentionally independent of v2 trait cohorts and v3 image data.
It consumes only the phenotype-blind source-observation CSV produced from the six
immutable iNaturalist metadata/API archives.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.rebuild_frozen_native_status import (
    fetch_distributions,
    fetch_name_resolution,
    load_pinned_level3,
    normalized,
    request_json,
)
from .build_native_range_join import classify_frame, load_contract, sha256_file


def read_source(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    required = {
        "obs_id",
        "source_taxon_name",
        "source_taxon_rank",
        "analysis_latitude",
        "analysis_longitude",
        "coordinate_status",
        "captive_state",
        "date_status",
        "observation_month",
        "position_accuracy_status",
        "position_accuracy_m",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Source-only observation table is missing: " + ", ".join(missing))
    if frame["obs_id"].isna().any() or frame["obs_id"].astype(str).duplicated().any():
        raise ValueError("Source-only observation IDs must be present and unique")
    return frame


def normalized_state(values: pd.Series) -> pd.Series:
    """Normalize CSV-round-tripped boolean/string state labels without imputing them."""
    return values.astype("string").str.strip().str.casefold()


def run_join(
    source_csv: Path,
    out_dir: Path,
    contract_path: Path,
    timeout: float = 30.0,
    retries: int = 4,
    sleep: float = 0.04,
) -> dict:
    contract = load_contract(contract_path)
    frame = read_source(source_csv)
    names = sorted(
        {normalized(name) for name in frame["source_taxon_name"].dropna() if normalized(name)}
    )
    source = contract["source_taxonomy_and_distribution"]
    metadata = request_json(
        f"https://api.gbif.org/v1/dataset/{source['dataset_key']}", timeout, retries
    )
    if normalized(metadata.get("doi")).casefold() != source["dataset_doi"].casefold():
        raise ValueError("WCVP dataset DOI differs from v3 contract")
    if normalized(metadata.get("modified")) != source["dataset_modified"]:
        raise ValueError("WCVP dataset modification timestamp differs from v3 contract")

    resolution = fetch_name_resolution(
        names, source["dataset_key"], timeout, retries, sleep
    )
    distributions = fetch_distributions(resolution, timeout, retries, sleep)
    legacy_shape = {
        "tdwg_level3": {
            "commit": contract["tdwg_level3"]["commit"],
            "path": contract["tdwg_level3"]["path"],
        }
    }
    geojson, tdwg_sha = load_pinned_level3(legacy_shape, timeout, retries)
    if tdwg_sha != contract["tdwg_level3"]["expected_sha256"]:
        raise ValueError("Pinned TDWG geometry SHA-256 differs from v3 contract")

    joined = classify_frame(frame, resolution, distributions, geojson)
    counts = joined["native_range_status"].value_counts(dropna=False).to_dict()
    taxon_counts = joined.groupby("native_range_status")["source_taxon_name"].nunique().to_dict()
    out_dir.mkdir(parents=True, exist_ok=False)
    joined.to_csv(out_dir / "native_range_join.csv", index=False)
    resolution.to_csv(out_dir / "taxon_resolution.csv", index=False)
    distributions.to_csv(out_dir / "wcvp_distribution_records.csv", index=False)

    primary = joined["native_range_status"].eq("native")
    primary &= normalized_state(joined["captive_state"]).eq("false")
    primary &= joined["date_status"].eq("exact_day")
    primary &= joined["coordinate_status"].eq("public_location_present_precision_not_gated")
    primary &= joined["taxon_resolution_status"].eq("resolved_unique_accepted_key")
    report = {
        "status": "FULL_SOURCE_NATIVE_RANGE_JOIN_COMPLETE_NOT_TRAIT_MODELLED",
        "contract_id": contract["contract_id"],
        "input_source_sha256": sha256_file(source_csv),
        "n_observations": int(len(joined)),
        "n_source_taxa": int(joined["source_taxon_name"].nunique()),
        "status_counts": {str(k): int(v) for k, v in counts.items()},
        "taxa_by_status": {str(k): int(v) for k, v in taxon_counts.items()},
        "primary_native_wild_public_exact_date_resolved_rows": int(primary.sum()),
        "source_rows_deleted": 0,
        "trait_files_read": 0,
        "image_files_read": 0,
        "ecological_models_executed": 0,
        "limits": [
            "Name resolution does not independently validate the photographed plant identification.",
            "Native status is regional at TDWG level 3 and cannot resolve fine range boundaries.",
            "This step does not decide the final environmental representation or trait-support thresholds.",
        ],
    }
    (out_dir / "native_range_join_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("analysis/v3/native_range_join_contract.json"),
    )
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    parser.add_argument("--max-retries", type=int, default=4)
    parser.add_argument("--sleep-seconds", type=float, default=0.04)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run_join(
        args.source_csv,
        args.out_dir,
        args.contract,
        args.timeout_seconds,
        args.max_retries,
        args.sleep_seconds,
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
