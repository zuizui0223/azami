#!/usr/bin/env python3
"""Grade image-quality sensitivity for PR70 continuous candidate traits.

GEB v2 uses the global/category-free continuous-trait screen to describe pattern.
This script is deliberately a sensitivity layer rather than a fail-closed claim
filter.  It fits all four CHELSA predictors together with image resolution and
sharpness, and reports sign behaviour across fixed resolution strata.
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

CLIMATE = ("chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15")
QUALITY = ("log_min_dimension", "log1p_sharpness")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--head", required=True)
    p.add_argument("--environment", required=True)
    p.add_argument("--contract", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-observations", type=int, default=100)
    p.add_argument("--min-taxa", type=int, default=10)
    return p.parse_args()


def zscore(values: pd.Series) -> pd.Series:
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("variable has no usable variance")
    return (values - values.mean()) / sd


def fit_within(frame: pd.DataFrame, endpoint: str, minimum: int, min_taxa: int) -> tuple[Any, pd.DataFrame]:
    columns = ["taxon_name", endpoint, *CLIMATE, *QUALITY]
    data = frame[columns].copy()
    for column in [endpoint, *CLIMATE, *QUALITY]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna()
    counts = data.groupby("taxon_name").size()
    data = data[data["taxon_name"].isin(counts[counts >= 2].index)].copy()
    if len(data) < minimum or data["taxon_name"].nunique() < min_taxa:
        raise ValueError("insufficient repeated observations")
    response = data[endpoint] - data.groupby("taxon_name")[endpoint].transform("mean")
    design = pd.DataFrame(index=data.index)
    for predictor in [*CLIMATE, *QUALITY]:
        centred = data[predictor] - data.groupby("taxon_name")[predictor].transform("mean")
        design[predictor] = zscore(centred)
    response = zscore(response)
    design = sm.add_constant(design)
    result = sm.OLS(response, design).fit(
        cov_type="cluster", cov_kwds={"groups": data["taxon_name"].to_numpy()}
    )
    return result, data


def coefficient_rows(result: Any, endpoint_id: str, module: str, cohort: str, data: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for predictor in [*CLIMATE, *QUALITY]:
        rows.append({
            "endpoint_id": endpoint_id,
            "module": module,
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


def main() -> None:
    args = parse_args()
    head = pd.read_csv(args.head, low_memory=False)
    env = pd.read_csv(args.environment, low_memory=False)
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)

    for frame in (head, env):
        frame["obs_id"] = frame["obs_id"].astype(str)
    if head["obs_id"].duplicated().any():
        raise ValueError("v2 candidate quality layer expects one selected head per observation")
    if env["obs_id"].duplicated().any():
        raise ValueError("environment table must be unique by obs_id")

    candidate = contract[contract["analysis_tier"].eq("candidate")].copy()
    endpoint_meta: list[dict[str, str]] = []
    for row in candidate.to_dict("records"):
        variable = row["measurement_variable"]
        status = f"{variable}_status"
        if variable in head.columns and status in head.columns:
            endpoint_meta.append({
                "endpoint_id": row["endpoint_id"],
                "module": row["module"],
                "variable": variable,
                "status": status,
            })

    env_payload = env[["obs_id", *CLIMATE]].copy()
    rows: list[dict[str, Any]] = []
    statuses: list[dict[str, Any]] = []

    for meta in endpoint_meta:
        usable = head[head[meta["status"]].astype(str).eq("usable")].copy()
        usable[meta["variable"]] = pd.to_numeric(usable[meta["variable"]], errors="coerce")
        usable["min_dimension"] = pd.to_numeric(usable["min_dimension"], errors="coerce")
        usable["sharpness"] = pd.to_numeric(usable["sharpness"], errors="coerce")
        joined = usable[["obs_id", "taxon_name", meta["variable"], "min_dimension", "sharpness"]].merge(
            env_payload, on="obs_id", how="left", validate="one_to_one"
        )
        joined = joined.rename(columns={meta["variable"]: meta["endpoint_id"]})
        joined["log_min_dimension"] = np.log(joined["min_dimension"])
        joined["log1p_sharpness"] = np.log1p(joined["sharpness"])
        joined["resolution_stratum"] = pd.cut(
            joined["min_dimension"],
            bins=[150.0 - 1e-9, 200.0, 300.0, np.inf],
            labels=["150_199", "200_299", "300_plus"],
            right=False,
        ).astype("string")

        cohorts: dict[str, pd.DataFrame] = {"all_quality_adjusted": joined}
        for label in ("150_199", "200_299", "300_plus"):
            cohorts[f"resolution_{label}"] = joined[joined["resolution_stratum"].eq(label)].copy()

        for cohort_name, cohort in cohorts.items():
            minimum = args.min_observations if cohort_name == "all_quality_adjusted" else 30
            try:
                result, data = fit_within(cohort, meta["endpoint_id"], minimum, args.min_taxa)
                rows.extend(coefficient_rows(result, meta["endpoint_id"], meta["module"], cohort_name, data))
                statuses.append({
                    "endpoint_id": meta["endpoint_id"], "module": meta["module"],
                    "cohort": cohort_name, "status": "success",
                    "n_observations": int(len(data)), "n_taxa": int(data["taxon_name"].nunique()),
                })
            except Exception as error:
                statuses.append({
                    "endpoint_id": meta["endpoint_id"], "module": meta["module"],
                    "cohort": cohort_name, "status": "failed",
                    "message": f"{type(error).__name__}: {error}",
                })

    coefficients = pd.DataFrame(rows)
    if coefficients.empty:
        coefficients = pd.DataFrame(columns=[
            "endpoint_id", "module", "cohort", "predictor", "beta_standardized",
            "p_value", "n_observations", "n_taxa",
        ])
    coefficients["q_fdr_bh_candidate_quality"] = np.nan
    full = coefficients[
        coefficients["cohort"].eq("all_quality_adjusted")
        & coefficients["predictor"].isin(CLIMATE)
    ].copy()
    if not full.empty:
        q = multipletests(full["p_value"].astype(float).to_numpy(), method="fdr_bh")[1]
        coefficients.loc[full.index, "q_fdr_bh_candidate_quality"] = q
    coefficients["fdr_significant_quality_0_05"] = coefficients["q_fdr_bh_candidate_quality"].lt(0.05)

    audit_rows: list[dict[str, Any]] = []
    for endpoint_id in sorted({m["endpoint_id"] for m in endpoint_meta}):
        module = next(m["module"] for m in endpoint_meta if m["endpoint_id"] == endpoint_id)
        for predictor in CLIMATE:
            full_row = coefficients[
                coefficients["endpoint_id"].eq(endpoint_id)
                & coefficients["cohort"].eq("all_quality_adjusted")
                & coefficients["predictor"].eq(predictor)
            ]
            if full_row.empty:
                continue
            full_one = full_row.iloc[0]
            strata = coefficients[
                coefficients["endpoint_id"].eq(endpoint_id)
                & coefficients["cohort"].str.startswith("resolution_")
                & coefficients["predictor"].eq(predictor)
            ]
            sign = np.sign(float(full_one["beta_standardized"]))
            same_sign = bool(len(strata) > 0 and (np.sign(strata["beta_standardized"].astype(float)) == sign).all())
            audit_rows.append({
                "endpoint_id": endpoint_id,
                "module": module,
                "predictor": predictor,
                "quality_adjusted_beta": float(full_one["beta_standardized"]),
                "quality_adjusted_ci_low": float(full_one["ci_low"]),
                "quality_adjusted_ci_high": float(full_one["ci_high"]),
                "quality_adjusted_q": float(full_one["q_fdr_bh_candidate_quality"]),
                "successful_resolution_strata": int(len(strata)),
                "same_sign_all_successful_strata": same_sign,
                "stratum_beta_min": float(strata["beta_standardized"].min()) if len(strata) else np.nan,
                "stratum_beta_max": float(strata["beta_standardized"].max()) if len(strata) else np.nan,
                "v2_interpretation": (
                    "quality_robust" if bool(full_one["fdr_significant_quality_0_05"]) and same_sign
                    else "image_sensitive_or_uncertain" if not same_sign and len(strata) > 0
                    else "quality_adjusted_not_fdr_supported"
                ),
            })

    audit = pd.DataFrame(audit_rows)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(out / "candidate_quality_adjusted_coefficients_v2.csv", index=False)
    pd.DataFrame(statuses).to_csv(out / "candidate_quality_model_status_v2.csv", index=False)
    audit.to_csv(out / "candidate_quality_sensitivity_v2.csv", index=False)
    report = {
        "n_registered_candidates": int(len(candidate)),
        "n_candidates_measurable_in_extended_head_table": int(len(endpoint_meta)),
        "n_quality_adjusted_climate_rows": int(len(full)),
        "n_quality_adjusted_fdr_signals": int(full["q_fdr_bh_candidate_quality"].lt(0.05).sum()) if not full.empty else 0,
        "interpretation": (
            "GEB v2 sensitivity layer: quality adjustment and resolution strata grade confidence; "
            "failure of this layer alone does not erase a global continuous image-phenotype pattern."
        ),
    }
    (out / "candidate_quality_sensitivity_report_v2.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
