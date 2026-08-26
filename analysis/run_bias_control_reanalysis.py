#!/usr/bin/env python3
"""Run outcome-blind seasonal and dominant-taxon sensitivity analyses.

This script does not replace the frozen Chapter 1 primary analysis.  It writes a
new reviewer-bias-control family.  The omission set is ranked only by observation
count, and the seasonal adjustment is fixed as taxon-specific first-harmonic
residualization of both outcome and climate predictor before model fitting.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


TRAITS = (
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
PREDICTORS = (
    "chelsa_bio01",
    "chelsa_bio04",
    "chelsa_bio12",
    "chelsa_bio15",
)
MODEL_DEFINITIONS = (
    "taxon_mean_only",
    "taxon_specific_cyclic_doy",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", required=True)
    parser.add_argument("--claims", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-species-model-n", type=int, default=5)
    parser.add_argument("--top-taxa", type=int, default=10)
    parser.add_argument("--baseline-tolerance", type=float, default=0.00075)
    return parser.parse_args()


def zscore(values: pd.Series) -> pd.Series:
    standard_deviation = float(values.std(ddof=1))
    if not math.isfinite(standard_deviation) or standard_deviation <= 0:
        raise ValueError("Variable has no usable variance")
    return values / standard_deviation


def add_cyclic_date(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    dates = pd.to_datetime(result["observed_on"], errors="coerce")
    if dates.isna().any():
        raise ValueError("Every bias-control row must have a valid observed_on date")
    day = dates.dt.dayofyear.astype(float)
    year_length = np.where(dates.dt.is_leap_year, 366.0, 365.0)
    phase = 2.0 * np.pi * (day - 1.0) / year_length
    result["doy_sin"] = np.sin(phase)
    result["doy_cos"] = np.cos(phase)
    return result


def eligible_complete_cases(
    frame: pd.DataFrame,
    trait: str,
    predictor: str,
    minimum: int,
) -> pd.DataFrame:
    columns = ["taxon_name", "observed_on", "doy_sin", "doy_cos", trait, predictor]
    work = frame[columns].copy()
    for column in (trait, predictor, "doy_sin", "doy_cos"):
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    counts = work.groupby("taxon_name").size()
    keep = counts.loc[counts.ge(minimum)].index
    work = work.loc[work["taxon_name"].isin(keep)].copy()
    variable = work.groupby("taxon_name")[[trait, predictor]].nunique()
    keep = variable.index[(variable[trait] >= 2) & (variable[predictor] >= 2)]
    return work.loc[work["taxon_name"].isin(keep)].copy()


def residualize_taxon_mean(work: pd.DataFrame, column: str) -> pd.Series:
    return work[column] - work.groupby("taxon_name")[column].transform("mean")


def residualize_taxon_cyclic(work: pd.DataFrame, column: str) -> pd.Series:
    output = pd.Series(index=work.index, dtype=float)
    for _, group in work.groupby("taxon_name", sort=False):
        design = np.column_stack([
            np.ones(len(group), dtype=float),
            group["doy_sin"].to_numpy(float),
            group["doy_cos"].to_numpy(float),
        ])
        values = group[column].to_numpy(float)
        coefficients, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
        output.loc[group.index] = values - design @ coefficients
    return output


def fit_model(
    frame: pd.DataFrame,
    trait: str,
    predictor: str,
    minimum: int,
    model_definition: str,
) -> dict[str, object]:
    work = eligible_complete_cases(frame, trait, predictor, minimum)
    if len(work) < 20 or work["taxon_name"].nunique() < 2:
        return {"status": "insufficient"}
    if model_definition == "taxon_mean_only":
        y_residual = residualize_taxon_mean(work, trait)
        x_residual = residualize_taxon_mean(work, predictor)
    elif model_definition == "taxon_specific_cyclic_doy":
        y_residual = residualize_taxon_cyclic(work, trait)
        x_residual = residualize_taxon_cyclic(work, predictor)
    else:
        raise ValueError(f"Unknown model definition: {model_definition}")
    y = zscore(y_residual)
    x = pd.DataFrame({"climate": zscore(x_residual)}, index=work.index)
    fit = sm.OLS(y, x).fit(
        cov_type="cluster", cov_kwds={"groups": work["taxon_name"]}
    )
    return {
        "status": "ok",
        "beta_std_within": float(fit.params["climate"]),
        "se_cluster_taxon": float(fit.bse["climate"]),
        "z": float(fit.tvalues["climate"]),
        "p_value": float(fit.pvalues["climate"]),
        "ci_low": float(fit.conf_int().loc["climate", 0]),
        "ci_high": float(fit.conf_int().loc["climate", 1]),
        "n_observations": int(len(work)),
        "n_taxa": int(work["taxon_name"].nunique()),
    }


def omission_scenarios(frame: pd.DataFrame, top_taxa: int) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    concentration = (
        frame.groupby("taxon_name", as_index=False)
        .size()
        .rename(columns={"size": "n_observations"})
        .sort_values(["n_observations", "taxon_name"], ascending=[False, True], kind="stable")
        .reset_index(drop=True)
    )
    concentration["rank"] = np.arange(1, len(concentration) + 1)
    concentration["fraction"] = concentration["n_observations"] / len(frame)
    selected = concentration.head(top_taxa)["taxon_name"].tolist()
    scenarios: list[dict[str, object]] = [{"scenario": "all_taxa", "omitted_taxa": []}]
    for rank, taxon in enumerate(selected, start=1):
        scenarios.append({
            "scenario": f"omit_rank_{rank:02d}",
            "omitted_taxa": [taxon],
        })
    if len(selected) >= 2:
        scenarios.append({
            "scenario": "omit_top_2_joint",
            "omitted_taxa": selected[:2],
        })
    return concentration, scenarios


def run_models(
    frame: pd.DataFrame,
    scenarios: Iterable[dict[str, object]],
    minimum: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for scenario in scenarios:
        omitted = list(scenario["omitted_taxa"])
        cohort = frame.loc[~frame["taxon_name"].isin(omitted)].copy()
        for model_definition in MODEL_DEFINITIONS:
            for trait in TRAITS:
                for predictor in PREDICTORS:
                    result = fit_model(cohort, trait, predictor, minimum, model_definition)
                    rows.append({
                        "scenario": scenario["scenario"],
                        "omitted_taxa": "|".join(omitted),
                        "model_definition": model_definition,
                        "trait": trait,
                        "predictor": predictor,
                        **result,
                    })
    table = pd.DataFrame(rows)
    table["q_fdr_bh"] = np.nan
    table["fdr_significant_0_05"] = False
    for _, indices in table.loc[table["status"].eq("ok")].groupby(
        ["scenario", "model_definition"]
    ).groups.items():
        p_values = table.loc[indices, "p_value"].to_numpy(float)
        reject, q_values, _, _ = multipletests(p_values, alpha=0.05, method="fdr_bh")
        table.loc[indices, "q_fdr_bh"] = q_values
        table.loc[indices, "fdr_significant_0_05"] = reject
    return table


def expected_supported_rows(claims: dict) -> pd.DataFrame:
    records = []
    frozen = claims["within_species_inference_levels"][
        "exhaustive_spatially_thinned_primary_36_component_models"
    ]
    predictor_names = {
        "BIO1": "chelsa_bio01",
        "BIO4": "chelsa_bio04",
        "BIO12": "chelsa_bio12",
        "BIO15": "chelsa_bio15",
    }
    for row in frozen["supported_non_circular_rows"]:
        endpoint = row["endpoint"]
        trait = {
            "orientation_angle": "orientation_angle_degrees_median",
            "corolla_chroma": "corolla_lab_chroma_median",
            "shape_aspect_ratio": "shape_aspect_ratio_median",
        }[endpoint]
        records.append({
            "endpoint": endpoint,
            "trait": trait,
            "predictor": predictor_names[row["predictor"]],
            "expected_beta": float(row["beta"]),
        })
    return pd.DataFrame(records)


def build_headline_audit(models: pd.DataFrame, expected: pd.DataFrame, tolerance: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    omission_names = [name for name in models["scenario"].unique() if name != "all_taxa"]
    for claim in expected.to_dict("records"):
        selection = models.loc[
            models["trait"].eq(claim["trait"])
            & models["predictor"].eq(claim["predictor"])
        ].copy()
        baseline = selection.loc[
            selection["scenario"].eq("all_taxa")
            & selection["model_definition"].eq("taxon_mean_only")
        ].iloc[0]
        adjusted = selection.loc[
            selection["scenario"].eq("all_taxa")
            & selection["model_definition"].eq("taxon_specific_cyclic_doy")
        ].iloc[0]
        omissions = selection.loc[
            selection["scenario"].isin(omission_names)
            & selection["model_definition"].eq("taxon_specific_cyclic_doy")
            & selection["status"].eq("ok")
        ]
        expected_sign = float(np.sign(claim["expected_beta"]))
        omission_sign_stable = bool(
            len(omissions) == len(omission_names)
            and (np.sign(omissions["beta_std_within"].to_numpy(float)) == expected_sign).all()
        )
        strict_retention = bool(
            abs(float(baseline["beta_std_within"]) - claim["expected_beta"]) <= tolerance
            and np.sign(float(adjusted["beta_std_within"])) == expected_sign
            and bool(adjusted["fdr_significant_0_05"])
            and omission_sign_stable
        )
        rows.append({
            **claim,
            "replicated_baseline_beta": float(baseline["beta_std_within"]),
            "baseline_absolute_difference": abs(float(baseline["beta_std_within"]) - claim["expected_beta"]),
            "season_adjusted_beta": float(adjusted["beta_std_within"]),
            "season_adjusted_ci_low": float(adjusted["ci_low"]),
            "season_adjusted_ci_high": float(adjusted["ci_high"]),
            "season_adjusted_q": float(adjusted["q_fdr_bh"]),
            "omission_beta_min": float(omissions["beta_std_within"].min()),
            "omission_beta_max": float(omissions["beta_std_within"].max()),
            "omission_sign_stable": omission_sign_stable,
            "strict_bias_control_retained": strict_retention,
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.observation, low_memory=False)
    required = {"obs_id", "taxon_name", "observed_on", *TRAITS, *PREDICTORS}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Observation table is missing columns: {sorted(missing)}")
    if frame["obs_id"].astype(str).duplicated().any():
        raise ValueError("Observation table must be unique by obs_id")
    frame = add_cyclic_date(frame)
    concentration, scenarios = omission_scenarios(frame, args.top_taxa)
    models = run_models(frame, scenarios, args.min_species_model_n)
    claims = json.loads(Path(args.claims).read_text(encoding="utf-8"))
    expected = expected_supported_rows(claims)
    audit = build_headline_audit(models, expected, args.baseline_tolerance)

    concentration.to_csv(output / "taxon_concentration.csv", index=False)
    models.to_csv(output / "season_and_taxon_omission_models.csv", index=False)
    audit.to_csv(output / "headline_bias_control_audit.csv", index=False)
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete_bias_control_reanalysis_not_stage_control",
        "n_observations": int(len(frame)),
        "n_taxa": int(frame["taxon_name"].nunique()),
        "observed_on_complete": bool(pd.to_datetime(frame["observed_on"], errors="coerce").notna().all()),
        "top_taxa": concentration.head(args.top_taxa).to_dict("records"),
        "model_definitions": {
            "baseline": "taxon-mean residualization of outcome and one climate predictor",
            "seasonal": "taxon-specific intercept plus sine/cosine day-of-year residualization of outcome and predictor",
            "uncertainty": "taxon-clustered robust standard errors",
            "multiplicity": "BH across 36 rows within each model-definition x omission scenario",
        },
        "strict_retention_rule": (
            "baseline reproduces the frozen beta within tolerance; season-adjusted full-data beta has the same sign "
            "and BH q<0.05; all predeclared top-taxon omission seasonal betas keep the frozen sign"
        ),
        "n_headline_rows_retained": int(audit["strict_bias_control_retained"].sum()),
        "n_headline_rows": int(len(audit)),
        "stage_control_status": "not_closed_day_of_year_is_not_an_image_stage_label",
    }
    (output / "bias_control_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
