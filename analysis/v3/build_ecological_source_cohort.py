"""Build the frozen observation-level source cohort for Chapter 1 v3 ecology.

This step consumes only the source-backed native-range join. It never reads image
measurements, trait values, environmental values or v2 candidate results. Unlike
the environment-diagnostic cohort it does not cap taxa or thin space: dominance
is handled downstream by the hierarchical model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "analysis" / "v3" / "ecological_source_cohort_contract.json"
ALLOWED_RANKS = {"species", "subspecies", "variety"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalized_state(values: pd.Series) -> pd.Series:
    return values.astype("string").str.strip().str.casefold()


def build(frame: pd.DataFrame, contract: dict) -> tuple[pd.DataFrame, dict]:
    required = {
        "obs_id", "accepted_key", "accepted_name", "source_taxon_name",
        "source_taxon_rank", "taxon_resolution_status", "native_range_status",
        "captive_state", "date_status", "coordinate_status", "analysis_latitude",
        "analysis_longitude", "observation_month", "position_accuracy_status",
        "position_accuracy_m",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Native-range join is missing: " + ", ".join(missing))
    if frame["obs_id"].isna().any() or frame["obs_id"].astype(str).duplicated().any():
        raise ValueError("obs_id must be present and unique")

    x = frame.copy()
    x["source_taxon_rank"] = x["source_taxon_rank"].fillna("").astype(str).str.strip().str.casefold()
    x["analysis_latitude"] = pd.to_numeric(x["analysis_latitude"], errors="coerce")
    x["analysis_longitude"] = pd.to_numeric(x["analysis_longitude"], errors="coerce")
    x["observation_month"] = pd.to_numeric(x["observation_month"], errors="coerce")
    x["accepted_key_numeric"] = pd.to_numeric(x["accepted_key"], errors="coerce").astype("Int64")

    eligibility = (
        x["native_range_status"].astype(str).eq("native")
        & x["taxon_resolution_status"].astype(str).eq("resolved_unique_accepted_key")
        & x["source_taxon_rank"].isin(ALLOWED_RANKS)
        & normalized_state(x["captive_state"]).eq("false")
        & x["date_status"].astype(str).eq("exact_day")
        & x["coordinate_status"].astype(str).eq("public_location_present_precision_not_gated")
        & x["accepted_key_numeric"].notna()
        & x["analysis_latitude"].between(-90, 90)
        & x["analysis_longitude"].between(-180, 180)
        & x["observation_month"].between(1, 12)
    )
    cohort = x.loc[eligibility].copy()
    cohort["accepted_key"] = cohort.pop("accepted_key_numeric").astype(str)
    cohort["quarter_degree_lat"] = (cohort["analysis_latitude"] * 4).floordiv(1).div(4)
    cohort["quarter_degree_lon"] = (cohort["analysis_longitude"] * 4).floordiv(1).div(4)
    cohort["quarter_degree_cell"] = cohort["quarter_degree_lat"].map(lambda v: f"{v:.2f}") + ":" + cohort["quarter_degree_lon"].map(lambda v: f"{v:.2f}")

    cohort = cohort.sort_values(["accepted_key", "obs_id"], kind="mergesort").reset_index(drop=True)
    expected = contract["expected_phenotype_blind_denominator_from_completed_upstream"]
    if len(cohort) != int(expected["rows"]):
        raise ValueError(f"Phenotype-blind ecological denominator changed: {len(cohort)} != {expected['rows']}")
    if cohort["accepted_key"].nunique() != int(expected["accepted_taxa"]):
        raise ValueError("Accepted-taxon denominator changed")

    support = cohort.groupby("accepted_key").agg(
        n_observations=("obs_id", "size"),
        n_quarter_degree_cells=("quarter_degree_cell", "nunique"),
        accepted_name=("accepted_name", "first"),
    ).sort_values(["n_observations", "accepted_name"], ascending=[False, True])
    n = len(cohort)
    support["row_fraction"] = support["n_observations"] / float(n)

    support_counts = support["n_observations"]
    cell_counts = support["n_quarter_degree_cells"]
    report = {
        "status": "FINAL_PHENOTYPE_BLIND_ECOLOGICAL_SOURCE_COHORT_BUILT_BEFORE_TRAIT_EXECUTION",
        "contract_id": contract["contract_id"],
        "rows": int(n),
        "accepted_taxa": int(len(support)),
        "row_cap": None,
        "taxon_cap": None,
        "spatial_thinning": False,
        "top_taxon_deletion": False,
        "trait_columns_read": 0,
        "environment_columns_read": 0,
        "v2_results_read": 0,
        "support": {
            "minimum_observations_per_taxon": int(support_counts.min()),
            "median_observations_per_taxon": float(support_counts.median()),
            "p90_observations_per_taxon": float(support_counts.quantile(.90)),
            "maximum_observations_per_taxon": int(support_counts.max()),
            "minimum_quarter_degree_cells_per_taxon": int(cell_counts.min()),
            "median_quarter_degree_cells_per_taxon": float(cell_counts.median()),
            "maximum_quarter_degree_cells_per_taxon": int(cell_counts.max()),
            "taxa_with_at_least_5_source_observations": int((support_counts >= 5).sum()),
            "taxa_with_at_least_10_source_observations": int((support_counts >= 10).sum()),
            "taxa_with_at_least_10_observations_and_4_cells": int(((support_counts >= 10) & (cell_counts >= 4)).sum()),
            "largest_raw_taxon_fraction": float(support["row_fraction"].max()),
        },
        "top_10_taxa_by_source_rows": [
            {
                "accepted_key": str(key),
                "accepted_name": str(row["accepted_name"]),
                "n_observations": int(row["n_observations"]),
                "n_quarter_degree_cells": int(row["n_quarter_degree_cells"]),
                "row_fraction": float(row["row_fraction"]),
            }
            for key, row in support.head(10).iterrows()
        ],
        "post_measurement_reporting_thresholds": contract["post_measurement_support_rules_fixed_before_coefficients"],
        "limits": [
            "This is the source-support cohort, not the realized endpoint-specific cohort after image availability and QC.",
            "Quarter-degree cells summarize support only and do not replace the required spatial model.",
            "Source taxon resolution is not independent image-identification validation.",
        ],
    }
    return cohort, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-join", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    if args.out_dir.exists():
        raise ValueError("Output exists; preserve prior cohort build")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if contract.get("status") != "fixed_before_original_trait_execution_and_trait_environment_join":
        raise ValueError("Ecological source cohort contract is not active")
    source = pd.read_csv(args.native_join, low_memory=False)
    cohort, report = build(source, contract)
    args.out_dir.mkdir(parents=True)
    keep = [
        "obs_id", "accepted_key", "accepted_name", "source_taxon_name", "source_taxon_rank",
        "analysis_latitude", "analysis_longitude", "observation_month",
        "position_accuracy_status", "position_accuracy_m", "quarter_degree_cell",
        "native_range_status",
    ]
    cohort[keep].to_csv(args.out_dir / "ecological_source_cohort.csv", index=False, lineterminator="\n")
    report["cohort_sha256"] = sha256_file(args.out_dir / "ecological_source_cohort.csv")
    (args.out_dir / "ecological_source_cohort_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
