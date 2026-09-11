#!/usr/bin/env python3
"""Audit QC-retention bias and join continuous primary traits to environment.

QC retention is an assessability outcome, not measurement accuracy. This script
keeps failed measurements explicit, tests whether retention varies with CHELSA
predictors while clustering by species, and writes observation/species tables
for downstream ecological models.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.stats.multitest import multipletests

TRAITS = ("colour", "shape", "orientation")
ENVIRONMENT = (
    "env_chelsa_bio01_native",
    "env_chelsa_bio04_native",
    "env_chelsa_bio12_native",
    "env_chelsa_bio15_native",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--head", required=True)
    p.add_argument("--observation", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--environment", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def numeric(frame: pd.DataFrame, columns: tuple[str, ...]) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    head = pd.read_csv(args.head)
    observation = pd.read_csv(args.observation)
    species = pd.read_csv(args.species)
    environment = pd.read_csv(args.environment)
    for frame in (head, observation, environment):
        frame["obs_id"] = frame["obs_id"].astype(str)
    if environment["obs_id"].duplicated().any() or observation["obs_id"].duplicated().any():
        raise ValueError("Observation/environment tables must be unique by obs_id")

    identity = observation[["obs_id", "taxon_name"]].merge(
        environment[["obs_id", "taxon_name"]], on="obs_id", suffixes=("_continuous", "_environment")
    )
    if len(identity) != len(observation) or not identity["taxon_name_continuous"].eq(identity["taxon_name_environment"]).all():
        raise ValueError("Observation IDs or taxon identities disagree between inputs")

    payload = environment.drop(columns=["taxon_name", "latitude", "longitude"], errors="ignore")
    head_env = head.merge(payload, on="obs_id", how="left", validate="many_to_one")
    numeric(head_env, ENVIRONMENT)
    for trait in TRAITS:
        required = f"{trait}_status"
        if required not in head_env:
            raise ValueError(f"Missing {required}")
        head_env[f"{trait}_usable"] = head_env[required].eq("usable").astype(int)

    overall_rows: list[dict[str, object]] = []
    for trait in TRAITS:
        counts = head_env[f"{trait}_status"].value_counts()
        failures = counts.drop(labels=["usable"], errors="ignore")
        usable = head_env[f"{trait}_usable"]
        overall_rows.append({
            "trait": trait,
            "n_heads": len(head_env),
            "n_usable": int(usable.sum()),
            "usable_fraction": float(usable.mean()),
            "n_failed_qc": int((1 - usable).sum()),
            "n_species_total": int(head_env["taxon_name"].nunique()),
            "n_species_any_usable": int(head_env.loc[usable.eq(1), "taxon_name"].nunique()),
            "n_observations_any_usable": int(head_env.loc[usable.eq(1), "obs_id"].nunique()),
            "most_common_failure_status": str(failures.idxmax()) if len(failures) else "",
            "most_common_failure_n": int(failures.max()) if len(failures) else 0,
        })
    overall = pd.DataFrame(overall_rows)
    overall.to_csv(out / "qc_overall_summary.csv", index=False)

    species_qc = head_env.groupby("taxon_name").agg(
        n_heads=("annotation_unit_id", "nunique"),
        n_observations=("obs_id", "nunique"),
        latitude_median=("latitude_numeric", "median"),
        longitude_median=("longitude_numeric", "median"),
        n_spatial_blocks=("spatial_block_10deg", "nunique"),
    ).reset_index()
    for trait in TRAITS:
        values = head_env.groupby("taxon_name")[f"{trait}_usable"].agg(["sum", "mean"]).reset_index()
        values.columns = ["taxon_name", f"n_usable_{trait}", f"usable_fraction_{trait}"]
        species_qc = species_qc.merge(values, on="taxon_name", validate="one_to_one")
        species_qc[f"no_usable_{trait}"] = species_qc[f"n_usable_{trait}"].eq(0)
        species_qc[f"low_retention_{trait}"] = species_qc["n_heads"].ge(5) & species_qc[f"usable_fraction_{trait}"].lt(0.5)
    species_qc.to_csv(out / "qc_species_retention.csv", index=False)

    blocks = head_env.groupby("spatial_block_10deg").agg(
        n_heads=("annotation_unit_id", "nunique"),
        n_observations=("obs_id", "nunique"),
        n_species=("taxon_name", "nunique"),
        latitude_median=("latitude_numeric", "median"),
        longitude_median=("longitude_numeric", "median"),
    ).reset_index()
    for trait in TRAITS:
        values = head_env.groupby("spatial_block_10deg")[f"{trait}_usable"].agg(["sum", "mean"]).reset_index()
        values.columns = ["spatial_block_10deg", f"n_usable_{trait}", f"usable_fraction_{trait}"]
        blocks = blocks.merge(values, on="spatial_block_10deg", validate="one_to_one")
        blocks[f"low_retention_{trait}"] = blocks["n_heads"].ge(10) & blocks[f"usable_fraction_{trait}"].lt(0.5)
    blocks.to_csv(out / "qc_spatial_block_retention.csv", index=False)

    bands = np.arange(-90, 91, 20)
    labels = [f"{bands[i]}_to_{bands[i + 1]}" for i in range(len(bands) - 1)]
    head_env["latitude_band_20deg"] = pd.cut(
        pd.to_numeric(head_env["latitude_numeric"], errors="coerce"),
        bins=bands, labels=labels, include_lowest=True, right=False,
    )
    latitude = head_env.groupby("latitude_band_20deg", observed=True).agg(
        n_heads=("annotation_unit_id", "nunique"),
        n_observations=("obs_id", "nunique"),
        n_species=("taxon_name", "nunique"),
    ).reset_index()
    for trait in TRAITS:
        values = head_env.groupby("latitude_band_20deg", observed=True)[f"{trait}_usable"].agg(["sum", "mean"]).reset_index()
        values.columns = ["latitude_band_20deg", f"n_usable_{trait}", f"usable_fraction_{trait}"]
        latitude = latitude.merge(values, on="latitude_band_20deg", validate="one_to_one")
    latitude.to_csv(out / "qc_latitude_band_retention.csv", index=False)

    smd_rows: list[dict[str, object]] = []
    gee_rows: list[dict[str, object]] = []
    for trait in TRAITS:
        usable = head_env[f"{trait}_usable"].astype(bool)
        response = usable.astype(float)
        for predictor in ENVIRONMENT:
            kept = head_env.loc[usable, predictor].dropna()
            failed = head_env.loc[~usable, predictor].dropna()
            pooled = np.sqrt(
                ((len(kept) - 1) * kept.var(ddof=1) + (len(failed) - 1) * failed.var(ddof=1))
                / max(len(kept) + len(failed) - 2, 1)
            )
            smd_rows.append({
                "trait": trait,
                "predictor": predictor,
                "n_usable": len(kept),
                "n_failed": len(failed),
                "usable_mean": kept.mean(),
                "failed_mean": failed.mean(),
                "usable_median": kept.median(),
                "failed_median": failed.median(),
                "standardized_mean_difference": (kept.mean() - failed.mean()) / pooled if pooled > 0 else np.nan,
            })
            x = head_env[predictor].astype(float)
            z = (x - x.mean()) / x.std(ddof=0)
            design = sm.add_constant(z.rename("z"))
            fitted = GEE(
                response, design, groups=head_env["taxon_name"],
                family=Binomial(), cov_struct=Exchangeable(),
            ).fit()
            coefficient = float(fitted.params["z"])
            standard_error = float(fitted.bse["z"])
            gee_rows.append({
                "trait": trait,
                "predictor": predictor,
                "model": "univariate_species_clustered_GEE",
                "n_heads": len(head_env),
                "n_species": head_env["taxon_name"].nunique(),
                "coefficient_per_predictor_sd": coefficient,
                "standard_error": standard_error,
                "odds_ratio_per_predictor_sd": np.exp(coefficient),
                "ci95_low_or": np.exp(coefficient - 1.96 * standard_error),
                "ci95_high_or": np.exp(coefficient + 1.96 * standard_error),
                "p_value": float(fitted.pvalues["z"]),
            })
    pd.DataFrame(smd_rows).to_csv(out / "qc_environment_standardized_differences.csv", index=False)
    gee = pd.DataFrame(gee_rows)
    gee["p_holm_across_12_tests"] = multipletests(gee["p_value"], method="holm")[1]
    gee["qc_environment_bias_flag"] = gee["p_holm_across_12_tests"].lt(0.05)
    gee.to_csv(out / "qc_environment_species_clustered_gee.csv", index=False)

    observation_env = observation.merge(payload, on="obs_id", how="left", validate="one_to_one")
    if observation_env["environment_primary_complete"].isna().any():
        raise ValueError("Some continuous observations lack environmental metadata")
    observation_env.to_csv(out / "primary_continuous_observation_environment.csv", index=False)

    numeric(environment, ENVIRONMENT + ("latitude_numeric", "longitude_numeric"))
    environment_species = environment.groupby("taxon_name").agg(
        n_environment_observations=("obs_id", "nunique"),
        latitude_species_median=("latitude_numeric", "median"),
        longitude_species_median=("longitude_numeric", "median"),
        n_spatial_blocks_10deg=("spatial_block_10deg", "nunique"),
        env_chelsa_bio01_species_median=("env_chelsa_bio01_native", "median"),
        env_chelsa_bio04_species_median=("env_chelsa_bio04_native", "median"),
        env_chelsa_bio12_species_median=("env_chelsa_bio12_native", "median"),
        env_chelsa_bio15_species_median=("env_chelsa_bio15_native", "median"),
    ).reset_index()
    species_env = species_qc.merge(species, on="taxon_name", how="left", validate="one_to_one").merge(
        environment_species, on="taxon_name", how="left", validate="one_to_one"
    )
    if len(species_env) != head_env["taxon_name"].nunique():
        raise ValueError("Species table lost taxa during joining")
    species_env.to_csv(out / "primary_continuous_species_environment.csv", index=False)

    endpoints = pd.DataFrame([
        ("orientation", "orientation_angle_degrees_median", "degrees"),
        ("colour", "corolla_lab_lightness_median", "Lab L*"),
        ("colour", "corolla_lab_chroma_median", "Lab chroma"),
        ("colour", "corolla_hue_sin_median", "unitless"),
        ("colour", "corolla_hue_cos_median", "unitless"),
        ("shape", "shape_aspect_ratio_median", "ratio"),
        ("shape", "shape_circularity_median", "unitless"),
        ("shape", "shape_solidity_median", "unitless"),
        ("shape", "shape_width_cv_median", "unitless"),
    ], columns=["trait_group", "observation_variable", "unit"])
    endpoints["analysis_tier"] = "main"
    endpoints.to_csv(out / "primary_continuous_endpoint_registry.csv", index=False)

    report = {
        "n_heads": int(len(head_env)),
        "n_observations": int(observation_env["obs_id"].nunique()),
        "n_species": int(head_env["taxon_name"].nunique()),
        "overall_qc": overall.to_dict("records"),
        "species_without_any_usable": {
            trait: sorted(species_qc.loc[species_qc[f"no_usable_{trait}"], "taxon_name"].tolist())
            for trait in TRAITS
        },
        "low_retention_species_n_heads_ge_5": {
            trait: int(species_qc[f"low_retention_{trait}"].sum()) for trait in TRAITS
        },
        "low_retention_spatial_blocks_n_heads_ge_10": {
            trait: int(blocks[f"low_retention_{trait}"].sum()) for trait in TRAITS
        },
        "environment_qc_bias_after_holm": gee.loc[
            gee["qc_environment_bias_flag"],
            ["trait", "predictor", "odds_ratio_per_predictor_sd", "p_holm_across_12_tests"],
        ].to_dict("records"),
        "semantic_status": (
            "QC-retention bias audit and environment join. Retention is assessability, "
            "not accuracy; failed measurements remain missing rather than biological absence."
        ),
    }
    (out / "qc_bias_and_join_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
