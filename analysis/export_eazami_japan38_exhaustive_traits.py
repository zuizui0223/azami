#!/usr/bin/env python3
"""Extract Japan-38 primary image-trait summaries from the frozen exhaustive layer.

Japan-38 paper concepts are reduced only to their genus+species binomial for a
coverage/sensitivity export. This intentionally does not assign broad image
records to a named variety or subspecies. The output is downstream evidence for
EAzami, not a new Chapter 1 result.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ENDPOINTS = (
    "orientation_angle_degrees_median",
    "corolla_lab_lightness_median",
    "corolla_lab_chroma_median",
    "corolla_hue_sin_median",
    "corolla_hue_cos_median",
    "shape_aspect_ratio_median",
    "shape_circularity_median",
    "shape_solidity_median",
    "shape_width_cv_median",
)


def binomial(value: object) -> str:
    text = re.sub(r"\s+", " ", str(value or "").strip())
    text = re.sub(r"^C\.\s+", "Cirsium ", text)
    tokens = text.split()
    if len(tokens) < 2:
        return text
    return " ".join([tokens[0], tokens[1].strip(",;()")])


def mad(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return float("nan")
    median = float(numeric.median())
    return float((numeric - median).abs().median())


def build(observations: pd.DataFrame, membership: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    required_obs = {"taxon_name", "obs_id", *ENDPOINTS}
    missing = sorted(required_obs.difference(observations.columns))
    if missing:
        raise ValueError(f"Observation table missing required columns: {missing}")
    if "paper_japan_member_id" not in membership.columns or "paper_taxon_concept" not in membership.columns:
        raise ValueError("Japan-38 membership table is incomplete")

    concepts = membership[["paper_japan_member_id", "paper_taxon_concept"]].copy()
    concepts["species_binomial"] = concepts["paper_taxon_concept"].map(binomial)
    concept_map = (
        concepts.groupby("species_binomial", sort=True)
        .agg(
            japan38_member_ids=("paper_japan_member_id", lambda x: "|".join(sorted(set(map(str, x))))),
            japan38_paper_concepts=("paper_taxon_concept", lambda x: "|".join(sorted(set(map(str, x))))),
            n_paper_concepts=("paper_japan_member_id", "nunique"),
        )
        .reset_index()
    )

    target_bins = set(concept_map["species_binomial"])
    focal = observations.loc[observations["taxon_name"].isin(target_bins)].copy()

    rows = []
    for taxon, part in focal.groupby("taxon_name", sort=True):
        row = {
            "taxon_name": taxon,
            "n_observations_detector_positive": int(len(part)),
            "n_spatial_blocks_5deg": int(part["spatial_block_5deg"].nunique()) if "spatial_block_5deg" in part.columns else 0,
            "coordinate_usable_n": int(part["coordinate_usable_bool"].fillna(False).astype(bool).sum()) if "coordinate_usable_bool" in part.columns else 0,
        }
        for endpoint in ENDPOINTS:
            values = pd.to_numeric(part[endpoint], errors="coerce")
            usable = values.dropna()
            row[f"{endpoint}_n"] = int(len(usable))
            row[f"{endpoint}_taxon_median"] = float(usable.median()) if len(usable) else np.nan
            row[f"{endpoint}_taxon_mad"] = mad(usable)
        rows.append(row)

    summary_table = pd.DataFrame(rows)
    if summary_table.empty:
        raise RuntimeError("No Japan-38 binomial taxa were found in exhaustive observations")
    summary_table = summary_table.merge(
        concept_map, left_on="taxon_name", right_on="species_binomial", how="left", validate="one_to_one"
    ).drop(columns=["species_binomial"])
    summary_table["trait_identity_scope"] = np.where(
        summary_table["n_paper_concepts"].gt(1),
        "binomial_only_multiple_paper_concepts_no_infraspecific_assignment",
        "binomial_coverage_paper_concept_may_still_include_infraspecific_name",
    )
    summary_table["claim_boundary"] = (
        "public-image trait summary; no variety assignment, genetic variance, ancestral state, evolutionary rate or adaptation inferred"
    )

    coverage = concept_map.copy()
    coverage["azami_detector_positive_trait_taxon_present"] = coverage["species_binomial"].isin(
        set(summary_table["taxon_name"])
    )
    coverage = coverage.merge(
        summary_table[["taxon_name", "n_observations_detector_positive"]],
        left_on="species_binomial",
        right_on="taxon_name",
        how="left",
    ).drop(columns=["taxon_name"])
    coverage["n_observations_detector_positive"] = coverage["n_observations_detector_positive"].fillna(0).astype(int)

    endpoint_n_columns = [f"{endpoint}_n" for endpoint in ENDPOINTS]
    summary = {
        "contract_version": "azami_eazami_japan38_exhaustive_traits_v1",
        "n_japan38_paper_concepts": int(concepts["paper_japan_member_id"].nunique()),
        "n_distinct_japan38_binomials": int(concepts["species_binomial"].nunique()),
        "n_distinct_binomial_trait_taxa_present": int(len(summary_table)),
        "n_paper_concepts_represented_at_binomial_level": int(
            coverage.loc[coverage["azami_detector_positive_trait_taxon_present"], "n_paper_concepts"].sum()
        ),
        "n_trait_taxa_with_10plus_detector_positive_observations": int(summary_table["n_observations_detector_positive"].ge(10).sum()),
        "n_trait_taxa_with_20plus_detector_positive_observations": int(summary_table["n_observations_detector_positive"].ge(20).sum()),
        "n_trait_taxa_with_50plus_detector_positive_observations": int(summary_table["n_observations_detector_positive"].ge(50).sum()),
        "endpoint_nonmissing_taxa": {
            endpoint: int(summary_table[f"{endpoint}_n"].gt(0).sum()) for endpoint in ENDPOINTS
        },
        "claim_boundary": (
            "Coverage and species-binomial image-trait summaries only. Broad image taxon concepts are not assigned to "
            "paper varieties/subspecies; outputs do not estimate evolutionary rate or adaptation."
        ),
    }
    return summary_table, coverage, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", required=True)
    parser.add_argument("--japan38", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    observations = pd.read_csv(args.observations, low_memory=False)
    membership = pd.read_csv(args.japan38, dtype=str, keep_default_na=False)
    traits, coverage, summary = build(observations, membership)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    traits.to_csv(out / "japan38_exhaustive_primary_trait_summary_v1.csv", index=False, encoding="utf-8")
    coverage.to_csv(out / "japan38_exhaustive_trait_coverage_v1.csv", index=False, encoding="utf-8")
    (out / "japan38_exhaustive_trait_summary_v1.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
