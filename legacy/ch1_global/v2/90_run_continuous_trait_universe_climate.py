#!/usr/bin/env python3
"""Run registry-driven within-taxon climate models for continuous traits.

Linear endpoints use taxon-demeaned standardized slopes with taxon-clustered
standard errors.  Circular endpoints are tested once per sine/cosine pair with a
within-taxon permutation of the predictor, avoiding the false claim that hue is
two independent biological traits.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--traits-long", required=True)
    parser.add_argument("--environment", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--predictors",
        nargs="+",
        default=["env_chelsa_bio01_native", "env_chelsa_bio04_native", "env_chelsa_bio12_native", "env_chelsa_bio15_native"],
    )
    parser.add_argument("--minimum-observations", type=int, default=100)
    parser.add_argument("--minimum-taxa", type=int, default=10)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260826)
    return parser.parse_args()


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _demean_standardize(frame: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="coerce")
    centred = values - values.groupby(frame["taxon_name"]).transform("mean")
    standard_deviation = float(centred.std(ddof=0))
    if not np.isfinite(standard_deviation) or standard_deviation <= 0:
        raise ValueError(f"No within-taxon variation in {column}")
    return centred / standard_deviation


def _fit_component(data: pd.DataFrame, response: str, predictor: str) -> dict[str, Any]:
    work = data[["taxon_name", response, predictor]].copy()
    for column in (response, predictor):
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    counts = work.groupby("taxon_name").size()
    work = work[work["taxon_name"].isin(counts[counts >= 2].index)].copy()
    y = _demean_standardize(work, response)
    x = _demean_standardize(work, predictor)
    result = sm.OLS(y.to_numpy(), x.to_numpy()[:, None]).fit(
        cov_type="cluster", cov_kwds={"groups": work["taxon_name"].to_numpy()}
    )
    return {
        "beta_std_within": float(result.params[0]),
        "standard_error_clustered": float(result.bse[0]),
        "confidence_low_95": float(result.conf_int()[0, 0]),
        "confidence_high_95": float(result.conf_int()[0, 1]),
        "p_value": float(result.pvalues[0]),
        "n_observations": int(len(work)),
        "n_taxa": int(work["taxon_name"].nunique()),
    }


def _permuted_within_taxon(values: np.ndarray, taxa: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    result = values.copy()
    for taxon in np.unique(taxa):
        index = np.flatnonzero(taxa == taxon)
        result[index] = rng.permutation(result[index])
    return result


def _slope(y: np.ndarray, x: np.ndarray) -> float:
    denominator = float(np.dot(x, x))
    return float(np.dot(x, y) / denominator) if denominator > 0 else float("nan")


def _fit_circular(
    data: pd.DataFrame,
    sine: str,
    cosine: str,
    predictor: str,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    work = data[["taxon_name", sine, cosine, predictor]].copy()
    for column in (sine, cosine, predictor):
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    counts = work.groupby("taxon_name").size()
    work = work[work["taxon_name"].isin(counts[counts >= 2].index)].copy()
    ys = _demean_standardize(work, sine).to_numpy()
    yc = _demean_standardize(work, cosine).to_numpy()
    x = _demean_standardize(work, predictor).to_numpy()
    taxa = work["taxon_name"].astype(str).to_numpy()
    beta_sine = _slope(ys, x)
    beta_cosine = _slope(yc, x)
    magnitude = float(math.hypot(beta_sine, beta_cosine))
    exceed = 0
    for _ in range(permutations):
        xp = _permuted_within_taxon(x, taxa, rng)
        if math.hypot(_slope(ys, xp), _slope(yc, xp)) >= magnitude - 1e-15:
            exceed += 1
    return {
        "beta_sine_std_within": beta_sine,
        "beta_cosine_std_within": beta_cosine,
        "effect_magnitude": magnitude,
        "effect_direction_degrees": float(math.degrees(math.atan2(beta_sine, beta_cosine)) % 360),
        "p_value": float((exceed + 1) / (permutations + 1)),
        "p_value_method": f"within_taxon_predictor_permutation_{permutations}",
        "n_observations": int(len(work)),
        "n_taxa": int(work["taxon_name"].nunique()),
    }


def _join_endpoint(traits: pd.DataFrame, environment: pd.DataFrame, endpoint_id: str, predictors: list[str]) -> pd.DataFrame:
    part = traits[traits["endpoint_id"].eq(endpoint_id)].copy()
    columns = ["obs_id", "taxon_name", "value"]
    part = part[columns].rename(columns={"value": endpoint_id})
    env = environment[["obs_id", *predictors]].copy()
    joined = part.merge(env, on="obs_id", how="inner", validate="one_to_one")
    return joined


def run_models(
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    predictors: list[str],
    minimum_observations: int,
    minimum_taxa: int,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    required = {"obs_id", "taxon_name", "endpoint_id", "analysis_tier", "analysis_eligible", "circular_group", "value"}
    missing = sorted(required.difference(traits.columns))
    if missing:
        raise ValueError(f"Trait long table is missing columns: {missing}")
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Trait long table must be unique by obs_id/endpoint_id")
    if environment["obs_id"].astype(str).duplicated().any():
        raise ValueError("Environment table must be unique by obs_id")
    absent_predictors = sorted(set(predictors).difference(environment.columns))
    if absent_predictors:
        raise ValueError(f"Environment table is missing predictors: {absent_predictors}")

    eligible = traits[_as_bool(traits["analysis_eligible"]) & traits["analysis_tier"].isin(["primary", "candidate"])].copy()
    # Empty CSV fields are parsed as NaN by pandas.  Without normalizing the
    # registry metadata after the on-disk round trip, every ordinary linear
    # endpoint is silently excluded because ``NaN.astype(str)`` becomes
    # ``"nan"`` rather than the contract's empty circular-group marker.
    eligible["circular_group"] = eligible["circular_group"].fillna("").astype(str).str.strip()
    endpoint_meta = eligible.drop_duplicates("endpoint_id").set_index("endpoint_id")
    coverage_rows: list[dict[str, Any]] = []
    result_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(seed)

    linear_endpoints = endpoint_meta[endpoint_meta["circular_group"].astype(str).eq("")].index.tolist()
    for endpoint_id in linear_endpoints:
        data = _join_endpoint(eligible, environment, endpoint_id, predictors)
        for predictor in predictors:
            complete = data.dropna(subset=[endpoint_id, predictor]).copy()
            repeat_counts = complete.groupby("taxon_name").size()
            complete = complete[complete["taxon_name"].isin(repeat_counts[repeat_counts >= 2].index)]
            supported = len(complete) >= minimum_observations and complete["taxon_name"].nunique() >= minimum_taxa
            coverage_rows.append({
                "endpoint_id": endpoint_id,
                "analysis_tier": endpoint_meta.loc[endpoint_id, "analysis_tier"],
                "inferential_unit": "linear_endpoint",
                "predictor": predictor,
                "n_observations": int(len(complete)),
                "n_taxa": int(complete["taxon_name"].nunique()),
                "status": "eligible" if supported else "insufficient_support",
            })
            if not supported:
                continue
            try:
                fit = _fit_component(complete, endpoint_id, predictor)
                result_rows.append({
                    "endpoint_id": endpoint_id,
                    "module": endpoint_meta.loc[endpoint_id, "module"],
                    "analysis_tier": endpoint_meta.loc[endpoint_id, "analysis_tier"],
                    "inferential_unit": "linear_endpoint",
                    "predictor": predictor,
                    "status": "ok",
                    "p_value_method": "taxon_clustered_ols",
                    **fit,
                })
            except Exception as error:
                result_rows.append({
                    "endpoint_id": endpoint_id,
                    "module": endpoint_meta.loc[endpoint_id, "module"],
                    "analysis_tier": endpoint_meta.loc[endpoint_id, "analysis_tier"],
                    "inferential_unit": "linear_endpoint",
                    "predictor": predictor,
                    "status": f"failed:{type(error).__name__}:{error}",
                })

    circular = endpoint_meta[endpoint_meta["circular_group"].astype(str).ne("")]
    for group, part in circular.groupby("circular_group"):
        endpoints = part.index.tolist()
        sine = next(value for value in endpoints if "sin" in value)
        cosine = next(value for value in endpoints if "cos" in value)
        sine_data = _join_endpoint(eligible, environment, sine, predictors)
        cosine_values = eligible[eligible["endpoint_id"].eq(cosine)][["obs_id", "value"]].rename(columns={"value": cosine})
        data = sine_data.merge(cosine_values, on="obs_id", how="inner", validate="one_to_one")
        for predictor in predictors:
            complete = data.dropna(subset=[sine, cosine, predictor]).copy()
            repeat_counts = complete.groupby("taxon_name").size()
            complete = complete[complete["taxon_name"].isin(repeat_counts[repeat_counts >= 2].index)]
            supported = len(complete) >= minimum_observations and complete["taxon_name"].nunique() >= minimum_taxa
            coverage_rows.append({
                "endpoint_id": group,
                "analysis_tier": part.iloc[0]["analysis_tier"],
                "inferential_unit": "circular_joint",
                "predictor": predictor,
                "n_observations": int(len(complete)),
                "n_taxa": int(complete["taxon_name"].nunique()),
                "status": "eligible" if supported else "insufficient_support",
            })
            if not supported:
                continue
            try:
                fit = _fit_circular(complete, sine, cosine, predictor, permutations, rng)
                result_rows.append({
                    "endpoint_id": group,
                    "module": part.iloc[0]["module"],
                    "analysis_tier": part.iloc[0]["analysis_tier"],
                    "inferential_unit": "circular_joint",
                    "predictor": predictor,
                    "status": "ok",
                    **fit,
                })
            except Exception as error:
                result_rows.append({
                    "endpoint_id": group,
                    "module": part.iloc[0]["module"],
                    "analysis_tier": part.iloc[0]["analysis_tier"],
                    "inferential_unit": "circular_joint",
                    "predictor": predictor,
                    "status": f"failed:{type(error).__name__}:{error}",
                })

    coefficients = pd.DataFrame(result_rows)
    coverage = pd.DataFrame(coverage_rows)
    if coefficients.empty:
        coefficients = pd.DataFrame(columns=[
            "endpoint_id", "module", "analysis_tier", "inferential_unit",
            "predictor", "status", "p_value",
        ])
    coefficients["q_fdr_bh_within_tier"] = np.nan
    ok = coefficients["status"].eq("ok") & pd.to_numeric(coefficients["p_value"], errors="coerce").notna()
    for tier, index in coefficients[ok].groupby("analysis_tier").groups.items():
        coefficients.loc[index, "q_fdr_bh_within_tier"] = multipletests(
            coefficients.loc[index, "p_value"].astype(float), method="fdr_bh"
        )[1]
    coefficients["fdr_significant_0_05"] = coefficients["q_fdr_bh_within_tier"].lt(0.05)
    # This registry-wide pass is a discovery screen. Endpoint-specific image
    # quality, seasonal/taxon-composition, native-range and reference-
    # validation gates are separate frozen analyses and cannot be inferred from
    # a small q value in this table.
    coefficients["screening_only"] = True
    coefficients["submission_claim_eligible"] = False
    coefficients["required_gate"] = np.where(
        coefficients["analysis_tier"].eq("candidate"),
        "endpoint_specific_image_quality_resolution_and_botanical_validation",
        "pr69_primary_timing_taxon_composition_native_range_and_external_validation",
    )
    report = {
        "n_trait_rows": int(len(traits)),
        "n_environment_rows": int(len(environment)),
        "n_linear_endpoints_modelled": int(coverage.loc[
            coverage["inferential_unit"].eq("linear_endpoint") & coverage["status"].eq("eligible"),
            "endpoint_id",
        ].nunique()),
        "n_circular_traits_modelled": int(coverage.loc[
            coverage["inferential_unit"].eq("circular_joint") & coverage["status"].eq("eligible"),
            "endpoint_id",
        ].nunique()),
        "n_successful_endpoint_predictor_models": int(coefficients["status"].eq("ok").sum()),
        "n_primary_fdr_signals": int((coefficients["analysis_tier"].eq("primary") & coefficients["fdr_significant_0_05"]).sum()),
        "n_candidate_fdr_signals": int((coefficients["analysis_tier"].eq("candidate") & coefficients["fdr_significant_0_05"]).sum()),
        "multiplicity": "BH within primary and candidate tiers separately; circular hue contributes one joint test per predictor",
        "claim_status": "screening_only_not_submission_claims",
        "submission_claims_created": 0,
        "required_gates": [
            "PR69 timing, dominant-taxon and native-range controls for primary rows",
            "endpoint-specific resolution and image-quality controls for candidate rows",
            "independent continuous botanical reference validation before proxy promotion",
        ],
        "interpretation_boundary": "Cross-sectional within-taxon image-phenotype association; not demonstrated plasticity or adaptation.",
    }
    return coefficients, coverage, report


def main() -> None:
    args = parse_args()
    traits = pd.read_csv(args.traits_long, low_memory=False)
    environment = pd.read_csv(args.environment, low_memory=False)
    coefficients, coverage, report = run_models(
        traits,
        environment,
        args.predictors,
        args.minimum_observations,
        args.minimum_taxa,
        args.permutations,
        args.seed,
    )
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(output / "continuous_trait_universe_climate_coefficients.csv", index=False, encoding="utf-8-sig")
    coverage.to_csv(output / "continuous_trait_universe_coverage.csv", index=False, encoding="utf-8-sig")
    (output / "continuous_trait_universe_climate_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
