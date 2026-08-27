#!/usr/bin/env python3
"""Materialize the frozen complete-18 observation cohort for process-environment extraction.

The whole-capitulum environmental analysis is defined on observations with all
18 measured primary/candidate endpoints available.  This utility applies that
measurement-completeness rule *without inspecting trait magnitudes* and subsets
the already frozen strict-spatial environment table to exactly those obs_id.
It changes I/O volume only; it does not change the estimand used by
run_capitulum_environment_blocks.py.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--out-csv", required=True, type=Path)
    p.add_argument("--report", required=True, type=Path)
    return p.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def complete18_ids(traits: pd.DataFrame) -> tuple[pd.Index, list[str]]:
    required = {
        "obs_id", "endpoint_id", "analysis_tier", "analysis_eligible", "value"
    }
    missing = required.difference(traits.columns)
    if missing:
        raise ValueError(f"Trait table missing columns: {sorted(missing)}")

    x = traits[
        traits["analysis_tier"].isin(["primary", "candidate"])
        & as_bool(traits["analysis_eligible"])
    ].copy()
    x["value"] = pd.to_numeric(x["value"], errors="coerce")
    x = x[x["value"].notna()]
    endpoints = sorted(x["endpoint_id"].dropna().astype(str).unique())
    if len(endpoints) != 18:
        raise ValueError(f"Expected 18 measured inferential endpoints, found {len(endpoints)}")
    if x.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Trait table must be unique by obs_id/endpoint_id")

    counts = x.groupby("obs_id")["endpoint_id"].nunique()
    ids = counts[counts.eq(len(endpoints))].index
    return ids, endpoints


def main() -> int:
    args = parse_args()
    traits = pd.read_csv(
        args.traits_long,
        usecols=["obs_id", "endpoint_id", "analysis_tier", "analysis_eligible", "value"],
        low_memory=False,
    )
    environment = pd.read_csv(args.environment, low_memory=False)
    if "obs_id" not in environment or "taxon_name" not in environment:
        raise ValueError("Environment table must contain obs_id and taxon_name")
    if environment["obs_id"].astype(str).duplicated().any():
        raise ValueError("Environment table must be unique by obs_id")

    ids, endpoints = complete18_ids(traits)
    cohort = environment[environment["obs_id"].isin(ids)].copy()
    if len(cohort) != len(ids):
        missing_ids = len(ids) - len(cohort)
        raise ValueError(f"Strict-spatial environment is missing {missing_ids} complete-18 obs_id")
    if len(cohort) == 0:
        raise ValueError("Complete-18 environment cohort is empty")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    cohort.to_csv(args.out_csv, index=False)
    report = {
        "selection_rule": "all 18 measured primary/candidate endpoints non-missing and analysis-eligible; no trait-value threshold",
        "n_measured_endpoints": len(endpoints),
        "endpoint_ids": endpoints,
        "n_observations": int(len(cohort)),
        "n_taxa": int(cohort["taxon_name"].nunique()),
        "source_environment_rows": int(len(environment)),
        "purpose": "I/O restriction only for the predeclared complete-18 multivariate environmental estimand",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
