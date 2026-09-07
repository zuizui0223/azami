"""Build a taxon-balanced, phenotype-blind cohort for environment diagnostics.

This cohort is used only to choose the abiotic representation before any trait
join. Eligibility is based on native-range/taxonomic/date/location metadata. A
stable hash cap prevents a few very common taxa from dominating correlations and
VIF diagnostics. No environmental value or trait value enters the selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


SPECIFICATION = {
    "version": "v3_environment_design_cohort_v1",
    "purpose": "phenotype-blind environmental representation diagnostics",
    "eligibility": [
        "native_range_status == native",
        "taxon_resolution_status == resolved_unique_accepted_key",
        "captive_state == false",
        "date_status == exact_day",
        "coordinate_status == public_location_present_precision_not_gated",
        "observation_month in 1..12",
        "accepted_key present",
    ],
    "selection": "SHA256(v3-env-design:20260907:obs_id) ordered within accepted taxon concept",
    "taxon_cap": 60,
    "minimum_taxon_support_for_diagnostics": 5,
    "meaning": "diagnostic sampling balance only; not the final ecological fitting cohort",
    "trait_values_read": 0,
    "environment_values_read": 0,
}


def stable_hash(obs_id: str) -> str:
    return hashlib.sha256(("v3-env-design:20260907:" + str(obs_id)).encode("utf-8")).hexdigest()


def build_design_cohort(
    frame: pd.DataFrame,
    taxon_cap: int = SPECIFICATION["taxon_cap"],
    minimum_taxon_support: int = SPECIFICATION["minimum_taxon_support_for_diagnostics"],
) -> tuple[pd.DataFrame, dict]:
    required = {
        "obs_id", "accepted_key", "accepted_name", "source_taxon_name",
        "taxon_resolution_status", "native_range_status", "captive_state",
        "date_status", "coordinate_status", "analysis_latitude",
        "analysis_longitude", "observation_month",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Native-range join is missing: " + ", ".join(missing))
    if frame["obs_id"].isna().any() or frame["obs_id"].astype(str).duplicated().any():
        raise ValueError("obs_id must be unique and present")
    x = frame.copy()
    x["observation_month"] = pd.to_numeric(x["observation_month"], errors="coerce")
    x["analysis_latitude"] = pd.to_numeric(x["analysis_latitude"], errors="coerce")
    x["analysis_longitude"] = pd.to_numeric(x["analysis_longitude"], errors="coerce")
    eligible = (
        x["native_range_status"].astype(str).eq("native")
        & x["taxon_resolution_status"].astype(str).eq("resolved_unique_accepted_key")
        & x["captive_state"].astype(str).eq("false")
        & x["date_status"].astype(str).eq("exact_day")
        & x["coordinate_status"].astype(str).eq("public_location_present_precision_not_gated")
        & x["accepted_key"].notna()
        & x["analysis_latitude"].between(-90, 90)
        & x["analysis_longitude"].between(-180, 180)
        & x["observation_month"].between(1, 12)
    )
    x = x.loc[eligible].copy()
    if x.empty:
        raise ValueError("No eligible native-range observations for environment design")
    x["accepted_key"] = x["accepted_key"].astype(str)
    support = x.groupby("accepted_key")["obs_id"].size().rename("eligible_taxon_n")
    retained_keys = support[support >= int(minimum_taxon_support)].index
    x = x.loc[x["accepted_key"].isin(retained_keys)].copy()
    if x.empty:
        raise ValueError("No taxa meet phenotype-blind environment support threshold")
    x["eligible_taxon_n"] = x["accepted_key"].map(support).astype(int)
    x["selection_hash"] = x["obs_id"].astype(str).map(stable_hash)
    x = x.sort_values(["accepted_key", "selection_hash", "obs_id"])
    x["taxon_design_rank"] = x.groupby("accepted_key").cumcount() + 1
    selected = x.loc[x["taxon_design_rank"] <= int(taxon_cap)].copy()
    selected_n = selected.groupby("accepted_key")["obs_id"].size().rename("selected_taxon_n")
    selected["selected_taxon_n"] = selected["accepted_key"].map(selected_n).astype(int)
    selected["equal_taxon_weight"] = 1.0 / selected["selected_taxon_n"].astype(float)

    total_by_taxon = selected.groupby("accepted_key")["equal_taxon_weight"].sum()
    if not ((total_by_taxon - 1.0).abs() < 1e-12).all():
        raise ValueError("Equal-taxon diagnostic weights do not sum to one per taxon")

    report = {
        "status": "PHENOTYPE_BLIND_TAXON_BALANCED_ENVIRONMENT_COHORT_BUILT",
        "specification": {
            **SPECIFICATION,
            "taxon_cap": int(taxon_cap),
            "minimum_taxon_support_for_diagnostics": int(minimum_taxon_support),
        },
        "eligible_rows_before_taxon_support": int(eligible.sum()),
        "eligible_taxa_before_support": int(support.size),
        "taxa_meeting_support": int(len(retained_keys)),
        "selected_rows": int(len(selected)),
        "selected_taxa": int(selected["accepted_key"].nunique()),
        "maximum_rows_per_taxon": int(selected_n.max()),
        "minimum_rows_per_taxon": int(selected_n.min()),
        "largest_selected_taxon_fraction": float(selected_n.max() / len(selected)),
        "trait_columns_read": 0,
        "environment_columns_read": 0,
        "ecological_models_executed": 0,
        "note": "This balances environment diagnostics only. Final ecological models retain their own hierarchical weighting/nesting contract and do not inherit this row cap automatically.",
    }
    return selected.reset_index(drop=True), report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-join", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--taxon-cap", type=int, default=SPECIFICATION["taxon_cap"])
    parser.add_argument("--minimum-taxon-support", type=int, default=SPECIFICATION["minimum_taxon_support_for_diagnostics"])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.out_dir.exists():
        raise ValueError("Output exists; preserve prior environment-design cohort")
    source = pd.read_csv(args.native_join, low_memory=False)
    selected, report = build_design_cohort(
        source, taxon_cap=args.taxon_cap, minimum_taxon_support=args.minimum_taxon_support
    )
    args.out_dir.mkdir(parents=True)
    keep = [
        "obs_id", "accepted_key", "accepted_name", "source_taxon_name",
        "analysis_latitude", "analysis_longitude", "observation_month",
        "eligible_taxon_n", "selected_taxon_n", "equal_taxon_weight",
        "taxon_design_rank", "selection_hash", "native_range_status",
    ]
    selected[keep].rename(
        columns={"analysis_latitude": "latitude", "analysis_longitude": "longitude"}
    ).to_csv(args.out_dir / "environment_design_cohort.csv", index=False)
    (args.out_dir / "environment_design_cohort_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
