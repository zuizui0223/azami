#!/usr/bin/env python3
"""Run the locked hemisphere-aware seasonal sensitivity family.

This analysis preserves the earlier raw-day-of-year results and adds two
predeclared checks for taxa sampled in both hemispheres: a half-cycle phase
alignment and separate hemisphere-specific cyclic curves.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import run_bias_control_reanalysis as primary  # noqa: E402


MODEL_DEFINITIONS = (
    "hemisphere_phase_aligned_cyclic_doy",
    "taxon_cyclic_doy_by_hemisphere",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", required=True)
    parser.add_argument("--claims", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-species-model-n", type=int, default=5)
    parser.add_argument("--top-taxa", type=int, default=10)
    return parser.parse_args()


def add_hemisphere_season_terms(frame: pd.DataFrame) -> pd.DataFrame:
    result = primary.add_cyclic_date(frame)
    latitude = pd.to_numeric(result["latitude"], errors="coerce")
    if latitude.isna().any() or (~latitude.between(-90, 90)).any():
        raise ValueError("Every row must have a valid latitude")
    result["hemisphere_south"] = latitude.lt(0).astype(float)
    raw_phase = np.arctan2(result["doy_sin"], result["doy_cos"])
    aligned_phase = raw_phase + np.pi * result["hemisphere_south"]
    result["aligned_doy_sin"] = np.sin(aligned_phase)
    result["aligned_doy_cos"] = np.cos(aligned_phase)
    return result


def eligible_complete_cases(
    frame: pd.DataFrame,
    trait: str,
    predictor: str,
    minimum: int,
) -> pd.DataFrame:
    columns = [
        "taxon_name", "observed_on", "latitude", "hemisphere_south",
        "doy_sin", "doy_cos", "aligned_doy_sin", "aligned_doy_cos",
        trait, predictor,
    ]
    work = frame[columns].copy()
    numeric = columns[2:8] + [trait, predictor]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    counts = work.groupby("taxon_name").size()
    keep = counts.loc[counts.ge(minimum)].index
    work = work.loc[work["taxon_name"].isin(keep)].copy()
    variable = work.groupby("taxon_name")[[trait, predictor]].nunique()
    keep = variable.index[(variable[trait] >= 2) & (variable[predictor] >= 2)]
    return work.loc[work["taxon_name"].isin(keep)].copy()


def residualize_group_design(
    work: pd.DataFrame,
    column: str,
    model_definition: str,
) -> pd.Series:
    output = pd.Series(index=work.index, dtype=float)
    for _, group in work.groupby("taxon_name", sort=False):
        if model_definition == "hemisphere_phase_aligned_cyclic_doy":
            design = np.column_stack([
                np.ones(len(group), dtype=float),
                group["aligned_doy_sin"].to_numpy(float),
                group["aligned_doy_cos"].to_numpy(float),
            ])
        elif model_definition == "taxon_cyclic_doy_by_hemisphere":
            south = group["hemisphere_south"].to_numpy(float)
            sine = group["doy_sin"].to_numpy(float)
            cosine = group["doy_cos"].to_numpy(float)
            if np.unique(south).size == 2:
                design = np.column_stack([
                    np.ones(len(group), dtype=float), sine, cosine,
                    south, south * sine, south * cosine,
                ])
            else:
                design = np.column_stack([
                    np.ones(len(group), dtype=float), sine, cosine,
                ])
        else:
            raise ValueError(f"Unknown model definition: {model_definition}")
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
    y_residual = residualize_group_design(work, trait, model_definition)
    x_residual = residualize_group_design(work, predictor, model_definition)
    y = primary.zscore(y_residual)
    x = pd.DataFrame({"climate": primary.zscore(x_residual)}, index=work.index)
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
            for trait in primary.TRAITS:
                for predictor in primary.PREDICTORS:
                    rows.append({
                        "scenario": scenario["scenario"],
                        "omitted_taxa": "|".join(omitted),
                        "model_definition": model_definition,
                        "trait": trait,
                        "predictor": predictor,
                        **fit_model(cohort, trait, predictor, minimum, model_definition),
                    })
    table = pd.DataFrame(rows)
    table["q_fdr_bh"] = np.nan
    table["fdr_significant_0_05"] = False
    for _, indices in table.loc[table["status"].eq("ok")].groupby(
        ["scenario", "model_definition"]
    ).groups.items():
        reject, q_values, _, _ = multipletests(
            table.loc[indices, "p_value"].to_numpy(float),
            alpha=0.05, method="fdr_bh",
        )
        table.loc[indices, "q_fdr_bh"] = q_values
        table.loc[indices, "fdr_significant_0_05"] = reject
    return table


def build_headline_audit(models: pd.DataFrame, expected: pd.DataFrame) -> pd.DataFrame:
    omission_names = [name for name in models["scenario"].unique() if name != "all_taxa"]
    rows: list[dict[str, object]] = []
    for claim in expected.to_dict("records"):
        record: dict[str, object] = dict(claim)
        expected_sign = float(np.sign(claim["expected_beta"]))
        retained = True
        for model_definition in MODEL_DEFINITIONS:
            selected = models.loc[
                models["trait"].eq(claim["trait"])
                & models["predictor"].eq(claim["predictor"])
                & models["model_definition"].eq(model_definition)
            ]
            full = selected.loc[selected["scenario"].eq("all_taxa")].iloc[0]
            omissions = selected.loc[
                selected["scenario"].isin(omission_names)
                & selected["status"].eq("ok")
            ]
            sign_stable = bool(
                len(omissions) == len(omission_names)
                and (np.sign(omissions["beta_std_within"].to_numpy(float)) == expected_sign).all()
            )
            prefix = model_definition
            record[f"{prefix}_beta"] = float(full["beta_std_within"])
            record[f"{prefix}_ci_low"] = float(full["ci_low"])
            record[f"{prefix}_ci_high"] = float(full["ci_high"])
            record[f"{prefix}_q"] = float(full["q_fdr_bh"])
            record[f"{prefix}_omission_beta_min"] = float(omissions["beta_std_within"].min())
            record[f"{prefix}_omission_beta_max"] = float(omissions["beta_std_within"].max())
            record[f"{prefix}_omission_sign_stable"] = sign_stable
            retained = retained and bool(
                np.sign(float(full["beta_std_within"])) == expected_sign
                and bool(full["fdr_significant_0_05"])
                and sign_stable
            )
        record["hemisphere_season_robust"] = retained
        rows.append(record)
    return pd.DataFrame(rows)


def hemisphere_composition(frame: pd.DataFrame) -> pd.DataFrame:
    result = (
        frame.groupby("taxon_name", as_index=False)
        .agg(
            n_observations=("taxon_name", "size"),
            n_northern=("hemisphere_south", lambda x: int((x == 0).sum())),
            n_southern=("hemisphere_south", lambda x: int((x == 1).sum())),
        )
    )
    result["sampled_both_hemispheres"] = result["n_northern"].gt(0) & result["n_southern"].gt(0)
    return result.sort_values(["n_observations", "taxon_name"], ascending=[False, True])


def main() -> None:
    args = parse_args()
    contract = json.loads(Path(args.contract).read_text(encoding="utf-8"))
    if contract.get("status") != "locked_before_hemisphere_sensitivity_outcome_execution":
        raise ValueError("Hemisphere sensitivity contract is not locked")
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.observation, low_memory=False)
    required = {"obs_id", "taxon_name", "observed_on", "latitude", *primary.TRAITS, *primary.PREDICTORS}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Observation table is missing columns: {sorted(missing)}")
    if frame["obs_id"].astype(str).duplicated().any():
        raise ValueError("Observation table must be unique by obs_id")
    frame = add_hemisphere_season_terms(frame)
    concentration, scenarios = primary.omission_scenarios(frame, args.top_taxa)
    models = run_models(frame, scenarios, args.min_species_model_n)
    claims = json.loads(Path(args.claims).read_text(encoding="utf-8"))
    expected = primary.expected_supported_rows(claims)
    audit = build_headline_audit(models, expected)
    composition = hemisphere_composition(frame)

    models.to_csv(output / "hemisphere_season_models.csv", index=False)
    audit.to_csv(output / "hemisphere_season_headline_audit.csv", index=False)
    composition.to_csv(output / "hemisphere_taxon_composition.csv", index=False)
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete_hemisphere_season_sensitivity_not_stage_control",
        "contract": str(Path(args.contract)),
        "n_observations": int(len(frame)),
        "n_taxa": int(frame["taxon_name"].nunique()),
        "n_southern_observations": int(frame["hemisphere_south"].sum()),
        "n_taxa_sampled_both_hemispheres": int(composition["sampled_both_hemispheres"].sum()),
        "model_definitions": list(MODEL_DEFINITIONS),
        "n_headline_rows_robust": int(audit["hemisphere_season_robust"].sum()),
        "n_headline_rows": int(len(audit)),
        "stage_control_status": "not_closed_day_of_year_is_not_an_image_stage_label",
    }
    (output / "hemisphere_season_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
