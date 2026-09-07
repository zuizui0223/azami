"""Build the v3 full-source taxon/native-range join before ecological fitting.

The input is the private full-source observation-annotation SQLite database. The
output is a local observation-level status table plus aggregate receipts. No
capitulum trait file is read and no source observation is deleted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sqlite3

import pandas as pd

from analysis.rebuild_frozen_native_status import (
    collapse_distribution_status,
    fetch_distributions,
    fetch_name_resolution,
    load_pinned_level3,
    normalized,
    assign_level3,
    request_json,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "analysis" / "v3" / "native_range_join_contract.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_contract(path: Path = DEFAULT_CONTRACT) -> dict:
    contract = json.loads(path.read_text(encoding="utf-8"))
    if contract.get("status") != "source_backed_join_before_ecological_fitting":
        raise ValueError("Native-range join contract is not active")
    return contract


def read_annotation_input(path: Path) -> pd.DataFrame:
    with sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True) as db:
        required = {
            row[1]
            for row in db.execute("PRAGMA table_info(annotations)")
        }
        needed = {
            "obs_id",
            "source_taxon_name",
            "source_taxon_rank",
            "analysis_latitude",
            "analysis_longitude",
            "coordinate_status",
            "captive_state",
            "date_status",
        }
        missing = sorted(needed - required)
        if missing:
            raise ValueError("Observation annotation database is missing: " + ", ".join(missing))
        frame = pd.read_sql_query(
            "SELECT obs_id,source_taxon_name,source_taxon_rank,analysis_latitude,analysis_longitude,coordinate_status,captive_state,date_status FROM annotations ORDER BY obs_id",
            db,
        )
    if frame["obs_id"].duplicated().any():
        raise ValueError("Observation annotation database has duplicate obs_id")
    return frame


def classify_frame(
    frame: pd.DataFrame,
    resolution: pd.DataFrame,
    distributions: pd.DataFrame,
    geojson: dict,
) -> pd.DataFrame:
    output = frame.copy()
    output["taxon_resolution_status"] = "not_attempted"
    output["accepted_key"] = pd.NA
    output["accepted_name"] = ""
    output["tdwg_level3_code"] = ""
    output["native_range_status"] = "not_classifiable"

    names = resolution.rename(columns={"input_name": "source_taxon_name"}).copy()
    name_cols = ["source_taxon_name", "resolution_status", "accepted_key", "accepted_name"]
    output = output.merge(names[name_cols], on="source_taxon_name", how="left", validate="many_to_one", suffixes=("", "_resolved"))
    output["taxon_resolution_status"] = output["resolution_status"].fillna("unresolved_or_missing_name")
    output["accepted_key"] = output["accepted_key_resolved"].where(output["accepted_key_resolved"].notna(), output["accepted_key"])
    output["accepted_name"] = output["accepted_name_resolved"].fillna(output["accepted_name"])
    output = output.drop(columns=["resolution_status", "accepted_key_resolved", "accepted_name_resolved"])

    coordinate_ok = output["coordinate_status"].eq("public_location_present_precision_not_gated")
    coordinate_ok &= pd.to_numeric(output["analysis_latitude"], errors="coerce").notna()
    coordinate_ok &= pd.to_numeric(output["analysis_longitude"], errors="coerce").notna()
    if coordinate_ok.any():
        geo_input = pd.DataFrame(
            {
                "latitude": pd.to_numeric(output.loc[coordinate_ok, "analysis_latitude"], errors="raise"),
                "longitude": pd.to_numeric(output.loc[coordinate_ok, "analysis_longitude"], errors="raise"),
            },
            index=output.index[coordinate_ok],
        )
        output.loc[coordinate_ok, "tdwg_level3_code"] = assign_level3(geo_input, geojson).astype(str)

    lookup = collapse_distribution_status(distributions)
    statuses: list[str] = []
    for row in output[["taxon_resolution_status", "accepted_key", "tdwg_level3_code", "coordinate_status"]].to_dict("records"):
        if row["taxon_resolution_status"] != "resolved_unique_accepted_key":
            statuses.append("unresolved_taxon")
            continue
        if row["coordinate_status"] != "public_location_present_precision_not_gated" or not normalized(row["tdwg_level3_code"]):
            statuses.append("unmapped_or_unusable_location")
            continue
        try:
            key = int(row["accepted_key"])
        except (TypeError, ValueError):
            statuses.append("unresolved_taxon")
            continue
        statuses.append(lookup.get(key, {}).get(str(row["tdwg_level3_code"]), "unlisted"))
    output["native_range_status"] = statuses
    return output


def run_join(
    annotation_db: Path,
    out_dir: Path,
    contract_path: Path = DEFAULT_CONTRACT,
    timeout: float = 30.0,
    retries: int = 4,
    sleep: float = 0.04,
) -> dict:
    contract = load_contract(contract_path)
    frame = read_annotation_input(annotation_db)
    names = sorted({normalized(name) for name in frame["source_taxon_name"].dropna() if normalized(name)})

    source = contract["source_taxonomy_and_distribution"]
    metadata = request_json(
        f"https://api.gbif.org/v1/dataset/{source['dataset_key']}",
        timeout,
        retries,
    )
    if normalized(metadata.get("doi")).casefold() != source["dataset_doi"].casefold():
        raise ValueError("WCVP dataset DOI differs from v3 contract")
    if normalized(metadata.get("modified")) != source["dataset_modified"]:
        raise ValueError("WCVP dataset modification timestamp differs from v3 contract")

    resolution = fetch_name_resolution(names, source["dataset_key"], timeout, retries, sleep)
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
    primary_candidate = joined["native_range_status"].eq("native")
    wild_state = joined["captive_state"].eq("false")
    exact_date = joined["date_status"].eq("exact_day")
    public_coord = joined["coordinate_status"].eq("public_location_present_precision_not_gated")

    out_dir.mkdir(parents=True, exist_ok=False)
    joined.to_csv(out_dir / "native_range_join.csv", index=False)
    resolution.to_csv(out_dir / "taxon_resolution.csv", index=False)
    distributions.to_csv(out_dir / "wcvp_distribution_records.csv", index=False)
    report = {
        "status": "FULL_SOURCE_NATIVE_RANGE_JOIN_COMPLETE_NOT_TRAIT_MODELLED",
        "contract_id": contract["contract_id"],
        "input_annotation_sha256": sha256_file(annotation_db),
        "n_observations": int(len(joined)),
        "n_source_taxa": int(joined["source_taxon_name"].nunique()),
        "status_counts": {str(k): int(v) for k, v in counts.items()},
        "taxa_by_status": {str(k): int(v) for k, v in taxon_counts.items()},
        "primary_native_candidate_rows": int(primary_candidate.sum()),
        "primary_native_wild_public_exact_date_rows": int((primary_candidate & wild_state & public_coord & exact_date).sum()),
        "source_rows_deleted": 0,
        "trait_files_read": 0,
        "ecological_models_executed": 0,
        "limits": [
            "Name resolution does not independently validate image identification.",
            "Native-range classification is regional (TDWG level 3) and cannot represent fine-scale range edges.",
            "This join does not decide final taxon-support, positional-accuracy or measurement-eligibility thresholds.",
        ],
    }
    (out_dir / "native_range_join_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-db", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    parser.add_argument("--max-retries", type=int, default=4)
    parser.add_argument("--sleep-seconds", type=float, default=0.04)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run_join(
        args.annotation_db,
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
