#!/usr/bin/env python3
"""Run staged continuous trait-environment models before phylogenetic sensitivity.

Observation models estimate within-species associations after demeaning outcome
and predictors within species, with species-clustered standard errors. Species
models are explicitly labelled pre-phylogenetic screens and must later be
rechecked with PGLS/alternative historical constraints.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

OBS_ENV = (
    "env_chelsa_bio01_native",
    "env_chelsa_bio04_native",
    "env_chelsa_bio12_native",
    "env_chelsa_bio15_native",
)
SPECIES_ENV = (
    "env_chelsa_bio01_species_median",
    "env_chelsa_bio04_species_median",
    "env_chelsa_bio12_species_median",
    "env_chelsa_bio15_species_median",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--observation", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--registry", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-observations", type=int, default=100)
    p.add_argument("--min-species", type=int, default=30)
    return p.parse_args()


def zscore(values: pd.Series) -> pd.Series:
    standard_deviation = values.std(ddof=0)
    if not np.isfinite(standard_deviation) or standard_deviation <= 0:
        raise ValueError("Variable has no usable variance")
    return (values - values.mean()) / standard_deviation


def fit_within_species(
    frame: pd.DataFrame,
    outcome: str,
    minimum: int,
) -> tuple[Any, pd.DataFrame]:
    columns = ["taxon_name", outcome, *OBS_ENV]
    data = frame[columns].copy()
    for column in [outcome, *OBS_ENV]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna()
    counts = data.groupby("taxon_name").size()
    data = data.loc[data["taxon_name"].isin(counts.loc[counts.ge(2)].index)].copy()
    if len(data) < minimum or data["taxon_name"].nunique() < 10:
        raise ValueError("Insufficient repeated observations for within-species model")

    response = data[outcome] - data.groupby("taxon_name")[outcome].transform("mean")
    response = zscore(response)
    design = pd.DataFrame(index=data.index)
    for predictor in OBS_ENV:
        centred = data[predictor] - data.groupby("taxon_name")[predictor].transform("mean")
        design[predictor] = zscore(centred)
    design = sm.add_constant(design)
    result = sm.OLS(response, design).fit(
        cov_type="cluster", cov_kwds={"groups": data["taxon_name"]}
    )
    return result, data


def fit_between_species(
    frame: pd.DataFrame,
    outcome: str,
    minimum: int,
) -> tuple[Any, pd.DataFrame]:
    data = frame[["taxon_name", outcome, *SPECIES_ENV]].copy()
    for column in [outcome, *SPECIES_ENV]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna()
    if len(data) < minimum:
        raise ValueError("Insufficient species for between-species model")
    response = zscore(data[outcome])
    design = pd.DataFrame({predictor: zscore(data[predictor]) for predictor in SPECIES_ENV})
    design = sm.add_constant(design)
    result = sm.OLS(response, design).fit(cov_type="HC3")
    return result, data


def coefficient_rows(
    result: Any,
    predictors: tuple[str, ...],
    metadata: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    for predictor in predictors:
        estimate = float(result.params[predictor])
        standard_error = float(result.bse[predictor])
        rows.append({
            **metadata,
            "predictor": predictor,
            "estimate_standardized": estimate,
            "standard_error": standard_error,
            "ci95_low": estimate - 1.96 * standard_error,
            "ci95_high": estimate + 1.96 * standard_error,
            "p_value": float(result.pvalues[predictor]),
        })
    return rows


def main() -> None:
    args = parse_args()
    observation = pd.read_csv(args.observation)
    species = pd.read_csv(args.species)
    registry = pd.read_csv(args.registry, dtype=str, keep_default_na=False)
    required = {
        "analysis_tier", "trait_group", "observation_variable", "species_variable",
    }
    if required.difference(registry.columns):
        raise ValueError("Endpoint registry is incomplete")

    cohorts = {
        "all_coordinate_eligible": observation,
        "positional_accuracy_le_10km": observation.loc[
            observation["coordinate_precision_tier"].isin(["high_le_1km", "moderate_1_to_10km"])
        ].copy(),
    }
    coefficients: list[dict[str, Any]] = []
    statuses: list[dict[str, Any]] = []

    for endpoint in registry.to_dict("records"):
        base = {
            "analysis_tier": endpoint["analysis_tier"],
            "trait_group": endpoint["trait_group"],
        }
        for cohort_name, cohort in cohorts.items():
            outcome = endpoint["observation_variable"]
            metadata = {
                **base,
                "endpoint": outcome,
                "model_scale": "observation_within_species",
                "cohort": cohort_name,
                "historical_constraint": "species_fixed_effect_via_within_species_demeaning",
                "interpretation_scope": "within_species_environment_association",
            }
            try:
                result, data = fit_within_species(cohort, outcome, args.min_observations)
                coefficients.extend(coefficient_rows(result, OBS_ENV, metadata))
                statuses.append({
                    **metadata,
                    "status": "success",
                    "n_rows": int(len(data)),
                    "n_species": int(data["taxon_name"].nunique()),
                    "r_squared": float(result.rsquared),
                })
            except Exception as error:
                statuses.append({**metadata, "status": "failed", "message": f"{type(error).__name__}: {error}"})

        outcome = endpoint["species_variable"]
        metadata = {
            **base,
            "endpoint": outcome,
            "model_scale": "species_between",
            "cohort": "all_species",
            "historical_constraint": "not_yet_applied_PGLS_required",
            "interpretation_scope": "prephylogenetic_between_species_screen_only",
        }
        try:
            result, data = fit_between_species(species, outcome, args.min_species)
            coefficients.extend(coefficient_rows(result, SPECIES_ENV, metadata))
            statuses.append({
                **metadata,
                "status": "success",
                "n_rows": int(len(data)),
                "n_species": int(len(data)),
                "r_squared": float(result.rsquared),
            })
        except Exception as error:
            statuses.append({**metadata, "status": "failed", "message": f"{type(error).__name__}: {error}"})

    coefficient_table = pd.DataFrame(coefficients)
    if coefficient_table.empty:
        raise RuntimeError("No continuous environment model succeeded")
    coefficient_table["p_fdr_within_tier_scale_cohort"] = coefficient_table.groupby(
        ["analysis_tier", "model_scale", "cohort"]
    )["p_value"].transform(lambda values: multipletests(values, method="fdr_bh")[1])
    coefficient_table["screening_signal_fdr_0_05"] = coefficient_table[
        "p_fdr_within_tier_scale_cohort"
    ].lt(0.05)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    coefficient_table.to_csv(out / "continuous_environment_model_coefficients.csv", index=False, encoding="utf-8-sig")
    status_table = pd.DataFrame(statuses)
    status_table.to_csv(out / "continuous_environment_model_status.csv", index=False, encoding="utf-8-sig")
    report = {
        "n_endpoints": int(len(registry)),
        "n_successful_models": int(status_table["status"].eq("success").sum()),
        "n_failed_models": int(status_table["status"].eq("failed").sum()),
        "n_coefficient_rows": int(len(coefficient_table)),
        "n_main_screening_signals": int((coefficient_table["analysis_tier"].eq("main") & coefficient_table["screening_signal_fdr_0_05"]).sum()),
        "n_auxiliary_screening_signals": int((coefficient_table["analysis_tier"].eq("auxiliary") & coefficient_table["screening_signal_fdr_0_05"]).sum()),
        "semantic_status": (
            "Pre-phylogenetic screening. Within-species models remove species-level historical "
            "differences by demeaning. Between-species coefficients are not final ecological "
            "claims and require later PGLS/tree sensitivity. Auxiliary endpoints remain exploratory."
        ),
    }
    (out / "continuous_environment_model_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
