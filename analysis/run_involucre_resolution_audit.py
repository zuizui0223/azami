#!/usr/bin/env python3
"""Reanalyse auxiliary involucre proxies with image-quality covariates."""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


ENDPOINTS = (
    "involucre_projection_roughness",
    "involucre_spread_fraction",
    "spine_relative_length_max_proxy",
)
CLIMATE = (
    "env_chelsa_bio01_native",
    "env_chelsa_bio04_native",
    "env_chelsa_bio12_native",
    "env_chelsa_bio15_native",
)
QUALITY = ("log_min_dimension", "log1p_sharpness")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--head", required=True)
    parser.add_argument("--environment", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-observations", type=int, default=100)
    return parser.parse_args()


def zscore(values: pd.Series) -> pd.Series:
    standard_deviation = float(values.std(ddof=0))
    if not math.isfinite(standard_deviation) or standard_deviation <= 0:
        raise ValueError("Variable has no usable variance")
    return (values - values.mean()) / standard_deviation


def prepare_observations(head: pd.DataFrame, environment: pd.DataFrame) -> pd.DataFrame:
    usable = head.loc[head["involucre_status"].eq("usable")].copy()
    for column in ["original_min_dimension", "original_sharpness", *ENDPOINTS]:
        usable[column] = pd.to_numeric(usable[column], errors="coerce")
    rows = []
    for obs_id, group in usable.groupby("obs_id", sort=True):
        row = {
            "obs_id": str(obs_id),
            "taxon_name": group["taxon_name"].iloc[0],
            "n_usable_heads": int(len(group)),
            "min_dimension_median": float(group["original_min_dimension"].median()),
            "sharpness_median": float(group["original_sharpness"].median()),
        }
        for endpoint in ENDPOINTS:
            row[endpoint] = float(group[endpoint].median())
        rows.append(row)
    observations = pd.DataFrame(rows)
    observations["obs_id"] = observations["obs_id"].astype(str)
    environment = environment.copy()
    environment["obs_id"] = environment["obs_id"].astype(str)
    payload = environment[[
        "obs_id", "taxon_name", "coordinate_precision_tier", *CLIMATE
    ]].drop_duplicates("obs_id")
    joined = observations.merge(
        payload.drop(columns="taxon_name"), on="obs_id", how="left", validate="one_to_one"
    )
    joined["log_min_dimension"] = np.log(joined["min_dimension_median"])
    joined["log1p_sharpness"] = np.log1p(joined["sharpness_median"])
    joined["resolution_stratum"] = pd.cut(
        joined["min_dimension_median"],
        bins=[150.0 - 1e-9, 200.0, 300.0, np.inf],
        labels=["150_199", "200_299", "300_plus"],
        right=False,
    ).astype("string")
    return joined


def fit_within(frame: pd.DataFrame, endpoint: str, minimum: int) -> tuple[object, pd.DataFrame]:
    columns = ["taxon_name", endpoint, *CLIMATE, *QUALITY]
    data = frame[columns].copy()
    for column in [endpoint, *CLIMATE, *QUALITY]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna()
    counts = data.groupby("taxon_name").size()
    data = data.loc[data["taxon_name"].isin(counts.loc[counts.ge(2)].index)].copy()
    if len(data) < minimum or data["taxon_name"].nunique() < 10:
        raise ValueError("Insufficient repeated observations")
    response = data[endpoint] - data.groupby("taxon_name")[endpoint].transform("mean")
    design = pd.DataFrame(index=data.index)
    for predictor in [*CLIMATE, *QUALITY]:
        centred = data[predictor] - data.groupby("taxon_name")[predictor].transform("mean")
        design[predictor] = zscore(centred)
    response = zscore(response)
    design = sm.add_constant(design)
    result = sm.OLS(response, design).fit(
        cov_type="cluster", cov_kwds={"groups": data["taxon_name"]}
    )
    return result, data


def coefficient_rows(result: object, endpoint: str, cohort: str, data: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for predictor in [*CLIMATE, *QUALITY]:
        rows.append({
            "endpoint": endpoint,
            "cohort": cohort,
            "predictor": predictor,
            "beta_standardized": float(result.params[predictor]),
            "se_cluster_taxon": float(result.bse[predictor]),
            "ci_low": float(result.conf_int().loc[predictor, 0]),
            "ci_high": float(result.conf_int().loc[predictor, 1]),
            "p_value": float(result.pvalues[predictor]),
            "n_observations": int(len(data)),
            "n_taxa": int(data["taxon_name"].nunique()),
        })
    return rows


def run_models(observations: pd.DataFrame, minimum: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    strict = observations.loc[
        observations["coordinate_precision_tier"].isin(["high_le_1km", "moderate_1_to_10km"])
    ].copy()
    rows: list[dict[str, object]] = []
    statuses: list[dict[str, object]] = []
    cohorts = {"all_resolution_adjusted": strict}
    for label in ("150_199", "200_299", "300_plus"):
        cohorts[f"resolution_{label}"] = strict.loc[strict["resolution_stratum"].eq(label)].copy()
    for cohort_name, cohort in cohorts.items():
        cohort_minimum = minimum if cohort_name == "all_resolution_adjusted" else 30
        for endpoint in ENDPOINTS:
            try:
                fit, data = fit_within(cohort, endpoint, cohort_minimum)
                rows.extend(coefficient_rows(fit, endpoint, cohort_name, data))
                statuses.append({
                    "endpoint": endpoint,
                    "cohort": cohort_name,
                    "status": "success",
                    "n_observations": int(len(data)),
                    "n_taxa": int(data["taxon_name"].nunique()),
                })
            except Exception as error:
                statuses.append({
                    "endpoint": endpoint,
                    "cohort": cohort_name,
                    "status": "failed",
                    "message": f"{type(error).__name__}: {error}",
                })
    coefficients = pd.DataFrame(rows)
    coefficients["q_fdr_bh_climate_12"] = np.nan
    full_climate = coefficients.loc[
        coefficients["cohort"].eq("all_resolution_adjusted")
        & coefficients["predictor"].isin(CLIMATE)
    ]
    if len(full_climate) == 12:
        reject, q_values, _, _ = multipletests(
            full_climate["p_value"].to_numpy(float), alpha=0.05, method="fdr_bh"
        )
        coefficients.loc[full_climate.index, "q_fdr_bh_climate_12"] = q_values
        coefficients.loc[full_climate.index, "fdr_significant_climate_0_05"] = reject
    else:
        coefficients["fdr_significant_climate_0_05"] = False
    return coefficients, pd.DataFrame(statuses)


def build_audit(coefficients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for endpoint in ENDPOINTS:
        full = coefficients.loc[
            coefficients["endpoint"].eq(endpoint)
            & coefficients["cohort"].eq("all_resolution_adjusted")
            & coefficients["predictor"].eq("env_chelsa_bio04_native")
        ].iloc[0]
        strata = coefficients.loc[
            coefficients["endpoint"].eq(endpoint)
            & coefficients["cohort"].str.startswith("resolution_")
            & coefficients["predictor"].eq("env_chelsa_bio04_native")
        ]
        positive_strata = bool(len(strata) > 0 and (strata["beta_standardized"] > 0).all())
        retained = bool(
            float(full["beta_standardized"]) > 0
            and bool(full.get("fdr_significant_climate_0_05", False))
            and positive_strata
        )
        rows.append({
            "endpoint": endpoint,
            "adjusted_bio4_beta": float(full["beta_standardized"]),
            "adjusted_bio4_ci_low": float(full["ci_low"]),
            "adjusted_bio4_ci_high": float(full["ci_high"]),
            "adjusted_bio4_q": float(full["q_fdr_bh_climate_12"]),
            "successful_resolution_strata": int(len(strata)),
            "resolution_stratum_beta_min": float(strata["beta_standardized"].min()) if len(strata) else np.nan,
            "resolution_stratum_beta_max": float(strata["beta_standardized"].max()) if len(strata) else np.nan,
            "positive_in_all_successful_strata": positive_strata,
            "strict_resolution_control_retained": retained,
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    head = pd.read_csv(args.head, low_memory=False)
    environment = pd.read_csv(args.environment, low_memory=False)
    observations = prepare_observations(head, environment)
    coefficients, statuses = run_models(observations, args.min_observations)
    audit = build_audit(coefficients)
    observations.to_csv(output / "involucre_resolution_observations.csv", index=False)
    coefficients.to_csv(output / "involucre_resolution_adjusted_models.csv", index=False)
    statuses.to_csv(output / "involucre_resolution_model_status.csv", index=False)
    audit.to_csv(output / "involucre_resolution_headline_audit.csv", index=False)
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete_resolution_and_sharpness_bias_control",
        "n_usable_heads": int(head["involucre_status"].eq("usable").sum()),
        "n_observations": int(len(observations)),
        "n_le_10km_observations": int(observations["coordinate_precision_tier"].isin([
            "high_le_1km", "moderate_1_to_10km"
        ]).sum()),
        "quality_covariates": list(QUALITY),
        "resolution_strata": ["150_199", "200_299", "300_plus"],
        "strict_retention_rule": (
            "positive full-data BIO4 coefficient with BH q<0.05 across 12 climate rows and positive BIO4 "
            "coefficient in every resolution stratum that has sufficient data"
        ),
        "n_endpoints_retained": int(audit["strict_resolution_control_retained"].sum()),
    }
    (output / "involucre_resolution_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
