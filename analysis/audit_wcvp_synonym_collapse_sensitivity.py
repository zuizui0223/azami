#!/usr/bin/env python3
"""Audit headline Chapter 1 results after collapsing WCVP synonym conflicts.

This is a sensitivity analysis, not a silent taxonomic rewrite. WCVP synonym
candidates are collapsed only when their accepted-name candidate resolves to an
already active frozen source name. The audit then recomputes the nested point
variance fractions and the 36 primary within-species climate models from frozen
trait/environment values.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

HEAD_TRAITS = {
    "orientation_angle_degrees": "orientation_status",
    "corolla_lab_lightness": "colour_status",
    "corolla_lab_chroma": "colour_status",
    "corolla_hue_sin": "colour_status",
    "corolla_hue_cos": "colour_status",
    "shape_aspect_ratio": "shape_status",
    "shape_circularity": "shape_status",
    "shape_solidity": "shape_status",
    "shape_width_cv": "shape_status",
}
OBS_TRAITS = [
    "orientation_angle_degrees_median",
    "corolla_lab_lightness_median",
    "corolla_lab_chroma_median",
    "corolla_hue_sin_median",
    "corolla_hue_cos_median",
    "shape_aspect_ratio_median",
    "shape_circularity_median",
    "shape_solidity_median",
    "shape_width_cv_median",
]
PREDICTORS = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--wcvp-review", required=True, type=Path)
    p.add_argument("--head-table", required=True, type=Path)
    p.add_argument("--environment-table", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--min-species-model-n", type=int, default=5)
    return p.parse_args()


def active_target(candidate: str, active_names: set[str]) -> str | None:
    candidate = str(candidate).strip()
    matches = [
        name for name in active_names
        if candidate == name or candidate.startswith(name + " ")
    ]
    return max(matches, key=len) if matches else None


def build_mapping(review: pd.DataFrame) -> dict[str, str]:
    active = set(review["input_name"].astype(str).str.strip())
    mapping: dict[str, str] = {}
    synonyms = review.loc[review["recommended_decision"].eq("synonym")]
    for _, row in synonyms.iterrows():
        source = str(row["input_name"]).strip()
        target = active_target(str(row["wcvp_accepted_name_candidate"]), active)
        if not target:
            raise ValueError(f"WCVP synonym target for {source!r} is not an active source unit")
        if source == target:
            raise ValueError(f"Synonym mapping is self-referential: {source}")
        mapping[source] = target
    if not mapping:
        raise ValueError("No WCVP synonym candidates found")
    return mapping


def nested_fractions(frame: pd.DataFrame, trait: str, status: str) -> dict[str, float | int]:
    work = frame.loc[
        frame[status].eq("usable"), ["taxon_name", "photo_id", trait]
    ].copy()
    work[trait] = pd.to_numeric(work[trait], errors="coerce")
    work = work.dropna(subset=[trait])
    y = work[trait].to_numpy(float)
    grand = float(np.mean(y))
    total = float(np.sum((y - grand) ** 2))
    photo = work.groupby(["taxon_name", "photo_id"])[trait].agg(["size", "mean"])
    species = work.groupby("taxon_name")[trait].agg(["size", "mean"])
    species_mean = species["mean"]
    species_ss = float(np.sum(species["size"] * (species["mean"] - grand) ** 2))
    pmeans = photo.index.get_level_values("taxon_name").map(species_mean)
    photo_ss = float(np.sum(photo["size"].to_numpy(float) * (photo["mean"].to_numpy(float) - np.asarray(pmeans, float)) ** 2))
    indexed = work.set_index(["taxon_name", "photo_id"])
    residual = indexed[trait] - indexed.index.map(photo["mean"])
    head_ss = float(np.sum(np.square(residual.to_numpy(float))))
    if not np.isclose(total, species_ss + photo_ss + head_ss, rtol=1e-8, atol=1e-8 * max(1.0, total)):
        raise RuntimeError(f"Nested sums of squares do not close for {trait}")
    return {
        "n_heads": int(len(work)),
        "n_taxa": int(work["taxon_name"].nunique()),
        "among_species_fraction": species_ss / total,
        "among_photos_within_species_fraction": photo_ss / total,
        "among_heads_within_photo_fraction": head_ss / total,
        "within_assigned_species_fraction": (photo_ss + head_ss) / total,
    }


def fit_one(df: pd.DataFrame, trait: str, predictor: str, min_n: int) -> dict[str, object]:
    work = df[["taxon_name", trait, predictor]].copy()
    work[trait] = pd.to_numeric(work[trait], errors="coerce")
    work[predictor] = pd.to_numeric(work[predictor], errors="coerce")
    work = work.dropna()
    keep = []
    for taxon, group in work.groupby("taxon_name"):
        if len(group) >= min_n and group[trait].nunique() >= 2 and group[predictor].nunique() >= 2:
            keep.append(taxon)
    work = work.loc[work["taxon_name"].isin(keep)].copy()
    work["y_dm"] = work[trait] - work.groupby("taxon_name")[trait].transform("mean")
    work["x_dm"] = work[predictor] - work.groupby("taxon_name")[predictor].transform("mean")
    y_sd = float(work["y_dm"].std(ddof=1))
    x_sd = float(work["x_dm"].std(ddof=1))
    if not (math.isfinite(y_sd) and math.isfinite(x_sd) and y_sd > 0 and x_sd > 0):
        raise ValueError(f"Zero variance for {trait} x {predictor}")
    fit = sm.OLS(work["y_dm"] / y_sd, work[["x_dm"]] / x_sd).fit(
        cov_type="cluster", cov_kwds={"groups": work["taxon_name"]}
    )
    return {
        "trait": trait,
        "predictor": predictor,
        "beta_std_within": float(fit.params["x_dm"]),
        "se_cluster_species": float(fit.bse["x_dm"]),
        "p_value": float(fit.pvalues["x_dm"]),
        "n_observations": int(len(work)),
        "n_species": int(work["taxon_name"].nunique()),
    }


def fit_all(df: pd.DataFrame, min_n: int) -> pd.DataFrame:
    out = pd.DataFrame([
        fit_one(df, trait, predictor, min_n)
        for trait in OBS_TRAITS
        for predictor in PREDICTORS
    ])
    reject, q, _, _ = multipletests(out["p_value"].to_numpy(float), alpha=0.05, method="fdr_bh")
    out["q_fdr_bh"] = q
    out["fdr_significant_0_05"] = reject
    return out


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    review = pd.read_csv(args.wcvp_review, low_memory=False)
    head = pd.read_csv(args.head_table, low_memory=False)
    env = pd.read_csv(args.environment_table, low_memory=False)
    mapping = build_mapping(review)

    mapping_rows = [
        {"source_unit": source, "wcvp_accepted_active_unit": target}
        for source, target in sorted(mapping.items())
    ]
    pd.DataFrame(mapping_rows).to_csv(args.out_dir / "wcvp_synonym_collapse_mapping.csv", index=False)

    nested_rows = []
    for trait, status in HEAD_TRAITS.items():
        original = nested_fractions(head, trait, status)
        collapsed_head = head.copy()
        collapsed_head["taxon_name"] = collapsed_head["taxon_name"].replace(mapping)
        collapsed = nested_fractions(collapsed_head, trait, status)
        nested_rows.append({
            "trait": trait,
            "original_n_taxa": original["n_taxa"],
            "collapsed_n_taxa": collapsed["n_taxa"],
            "original_within_fraction": original["within_assigned_species_fraction"],
            "collapsed_within_fraction": collapsed["within_assigned_species_fraction"],
            "delta_within_fraction": collapsed["within_assigned_species_fraction"] - original["within_assigned_species_fraction"],
        })
    nested = pd.DataFrame(nested_rows)
    nested.to_csv(args.out_dir / "wcvp_collapse_nested_variance_sensitivity.csv", index=False)

    original_models = fit_all(env, args.min_species_model_n)
    collapsed_env = env.copy()
    collapsed_env["taxon_name"] = collapsed_env["taxon_name"].replace(mapping)
    before_dedup = len(collapsed_env)
    if "analysis_cell" in collapsed_env.columns:
        collapsed_env["_obs_sort"] = pd.to_numeric(collapsed_env["obs_id"], errors="coerce")
        collapsed_env = (
            collapsed_env.sort_values(["taxon_name", "analysis_cell", "_obs_sort"])
            .drop_duplicates(["taxon_name", "analysis_cell"], keep="first")
            .drop(columns="_obs_sort")
        )
    collapsed_models = fit_all(collapsed_env, args.min_species_model_n)

    comparison = original_models.merge(
        collapsed_models,
        on=["trait", "predictor"],
        suffixes=("_original", "_collapsed"),
        validate="one_to_one",
    )
    comparison["delta_beta"] = comparison["beta_std_within_collapsed"] - comparison["beta_std_within_original"]
    comparison["sign_same"] = np.sign(comparison["beta_std_within_collapsed"]) == np.sign(comparison["beta_std_within_original"])
    comparison["fdr_decision_same"] = comparison["fdr_significant_0_05_collapsed"] == comparison["fdr_significant_0_05_original"]
    comparison.to_csv(args.out_dir / "wcvp_collapse_primary_climate_sensitivity.csv", index=False)

    report = {
        "n_wcvp_synonym_conflicts_collapsed": len(mapping),
        "mapping": mapping,
        "balanced_atlas_original_taxa": int(head["taxon_name"].nunique()),
        "balanced_atlas_collapsed_taxa": int(head["taxon_name"].replace(mapping).nunique()),
        "exhaustive_original_rows": int(len(env)),
        "exhaustive_collapsed_rows_after_cell_dedup": int(len(collapsed_env)),
        "exhaustive_rows_removed_by_postcollapse_cell_dedup": int(before_dedup - len(collapsed_env)),
        "exhaustive_original_taxa": int(env["taxon_name"].nunique()),
        "exhaustive_collapsed_taxa": int(collapsed_env["taxon_name"].nunique()),
        "original_fdr_supported_component_rows": int(original_models["fdr_significant_0_05"].sum()),
        "collapsed_fdr_supported_component_rows": int(collapsed_models["fdr_significant_0_05"].sum()),
        "n_fdr_decisions_changed": int((~comparison["fdr_decision_same"]).sum()),
        "n_coefficient_signs_changed": int((~comparison["sign_same"]).sum()),
        "max_absolute_beta_change": float(comparison["delta_beta"].abs().max()),
        "nested_original_within_fraction_min": float(nested["original_within_fraction"].min()),
        "nested_collapsed_within_fraction_min": float(nested["collapsed_within_fraction"].min()),
        "nested_max_absolute_within_fraction_change": float(nested["delta_within_fraction"].abs().max()),
        "interpretation": "sensitivity only; primary source-platform taxonomic units are not silently relabelled by this audit",
    }
    (args.out_dir / "wcvp_synonym_collapse_sensitivity_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
