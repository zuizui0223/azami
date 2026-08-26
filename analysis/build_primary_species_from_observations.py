#!/usr/bin/env python3
"""Derive a contract-compatible species table from strict observation-level traits.

The frozen exhaustive artifact stores several continuous endpoints directly at
observation level. GEB v2 must preserve every contract endpoint already present
there (not only the nine legacy primary endpoints), so candidate display and
other source-derived continuous measurements are not silently dropped before
building the unified trait universe.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--min-observations-per-taxon", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    observation = pd.read_csv(args.observation, low_memory=False)
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    if "obs_id" not in observation or "taxon_name" not in observation:
        raise ValueError("Observation table must contain obs_id and taxon_name")
    if observation["obs_id"].astype(str).duplicated().any():
        raise ValueError("Observation table must be unique by obs_id")

    present_contract = contract.loc[
        contract["observation_variable"].isin(observation.columns)
    ].copy()
    rows: list[dict[str, object]] = []
    for taxon_name, group in observation.groupby("taxon_name", sort=True):
        row: dict[str, object] = {
            "taxon_name": taxon_name,
            "n_observations_per_species": int(group["obs_id"].nunique()),
        }
        for endpoint in present_contract.to_dict("records"):
            obs_var = endpoint["observation_variable"]
            species_var = endpoint["species_variable"]
            status = endpoint["source_status_column"]
            values = pd.to_numeric(group[obs_var], errors="coerce")
            row[species_var] = float(values.median()) if values.notna().any() else np.nan
            if status in group:
                usable = pd.to_numeric(group[status], errors="coerce").fillna(0)
                row[status] = max(int(row.get(status, 0)), int((usable > 0).sum()))
        rows.append(row)

    species = pd.DataFrame(rows)
    species = species.loc[
        species["n_observations_per_species"].ge(args.min_observations_per_taxon)
    ].copy()
    if species["taxon_name"].astype(str).duplicated().any():
        raise RuntimeError("Derived species table is not unique by taxon_name")

    output = Path(args.out_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    species.to_csv(output, index=False, encoding="utf-8-sig")
    available_by_tier = (
        present_contract.groupby("analysis_tier")["endpoint_id"].nunique().astype(int).to_dict()
        if not present_contract.empty else {}
    )
    report = {
        "source_observations": int(len(observation)),
        "source_taxa": int(observation["taxon_name"].nunique()),
        "derived_species_rows": int(len(species)),
        "contract_endpoints_present_in_source": int(len(present_contract)),
        "source_endpoint_counts_by_tier": available_by_tier,
        "source_endpoint_ids": present_contract["endpoint_id"].tolist(),
        "min_observations_per_taxon": int(args.min_observations_per_taxon),
        "semantic_status": (
            "Species medians derived from every contract continuous endpoint already present in the same strict spatial "
            "observation cohort; no phenotype-based selection is introduced."
        ),
    }
    Path(args.report).write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
