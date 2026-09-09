#!/usr/bin/env python3
"""Full 46,276-observation species-level multivariable Chapter 1 lane.

Pipeline:
  frozen 46,276 strict-spatial observations
  -> all 27 registered endpoints (including five historical aggregation repairs)
  -> phenotype-blind VIFstep<=10 retained environmental predictors
  -> one simultaneous multivariable species-level model per trait unit
  -> broad spatial sensitivity using the same environmental model + spherical basis
  -> 52-tree multivariable Pagel-lambda PGLS sensitivity.

This lane deliberately does not split inference into separate within-/among-taxon
models. Observation-level records are the measurement base; the ecological
estimand is a species-level comparison based on taxon medians, which keeps the
spatial and phylogenetic sensitivity stages on the same estimand.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, t as student_t

from analysis.run_geb_v2_full27_environment_atlas import (
    as_bool,
    bh_adjust,
    inferential_units,
    stable_rng,
    standardize,
)
from analysis.run_geb_v2_full27_historical_sensitivity import (
    cholesky_with_jitter,
    load_trees,
    tree_covariance,
)
from analysis.run_geb_v2_full27_spatial_sensitivity import moran_test, spherical_basis


RECOVERED_ENDPOINTS = [
    "visible_floret_fraction",
    "corolla_white_pixel_fraction",
    "corolla_redmagenta_pixel_fraction",
    "corolla_purple_pixel_fraction",
    "corolla_yellow_pixel_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint-contract", required=True, type=Path)
    parser.add_argument("--traits-long", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--recovered-display", required=True, type=Path)
    parser.add_argument("--tree-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--predictors", nargs="+", required=True)
    parser.add_argument("--minimum-trait-observations-per-taxon", type=int, default=5)
    parser.add_argument("--minimum-taxa", type=int, default=30)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--moran-permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260909)
    return parser.parse_args()


def load_and_materialize_traits(
    traits_path: Path,
    recovered_path: Path,
    contract: pd.DataFrame,
    environment_ids: set[str],
) -> tuple[pd.DataFrame, dict[str, int]]:
    columns = [
        "obs_id",
        "taxon_name",
        "endpoint_id",
        "module",
        "analysis_tier",
        "measurement_available",
        "value",
    ]
    traits = pd.read_csv(traits_path, usecols=columns, low_memory=False)
    traits["obs_id"] = traits["obs_id"].astype(str)
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    traits = traits[as_bool(traits["measurement_available"]) & traits["value"].notna()].copy()
    traits = traits[traits["obs_id"].isin(environment_ids)].copy()

    recovered = pd.read_csv(recovered_path, low_memory=False)
    recovered["obs_id"] = recovered["obs_id"].astype(str)
    recovered = recovered[recovered["obs_id"].isin(environment_ids)].copy()

    metadata = contract.set_index("endpoint_id")
    appended: list[pd.DataFrame] = []
    counts: dict[str, int] = {}
    for endpoint in RECOVERED_ENDPOINTS:
        if endpoint not in recovered.columns:
            raise ValueError(f"Recovered table lacks {endpoint}")
        if endpoint not in metadata.index:
            raise ValueError(f"Endpoint contract lacks {endpoint}")
        existing_ids = set(traits.loc[traits["endpoint_id"].eq(endpoint), "obs_id"])
        value = pd.to_numeric(recovered[endpoint], errors="coerce")
        take = recovered.loc[value.notna() & ~recovered["obs_id"].isin(existing_ids), ["obs_id"]].copy()
        take["taxon_name"] = take["obs_id"].map(
            traits.drop_duplicates("obs_id").set_index("obs_id")["taxon_name"]
        )
        missing_taxon = take["taxon_name"].isna()
        if missing_taxon.any():
            # Taxon labels for repaired endpoints come from the environment table later.
            take = take.loc[~missing_taxon].copy()
        row = metadata.loc[endpoint]
        take["endpoint_id"] = endpoint
        take["module"] = row["module"]
        take["analysis_tier"] = row["analysis_tier"]
        take["measurement_available"] = True
        take["value"] = take["obs_id"].map(recovered.set_index("obs_id")[endpoint])
        counts[endpoint] = int(len(take))
        appended.append(take[columns])
    if appended:
        traits = pd.concat([traits, *appended], ignore_index=True)
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Materialized trait table is not unique by obs_id/endpoint_id")
    return traits, counts


def load_environment(path: Path, predictors: list[str]) -> pd.DataFrame:
    columns = ["obs_id", "taxon_name", "latitude", "longitude", *predictors]
    env = pd.read_csv(path, usecols=columns, low_memory=False)
    env["obs_id"] = env["obs_id"].astype(str)
    if len(env) != 46276 or env["obs_id"].nunique() != 46276:
        raise ValueError(f"Expected frozen 46,276 environment rows, found {len(env)}")
    for column in ["latitude", "longitude", *predictors]:
        env[column] = pd.to_numeric(env[column], errors="coerce")
    return env


def taxon_environment(environment: pd.DataFrame, predictors: list[str]) -> pd.DataFrame:
    columns = [*predictors, "latitude", "longitude"]
    result = environment.groupby("taxon_name")[columns].median(numeric_only=True)
    return result.dropna(subset=predictors)


def unit_taxon_data(
    unit: dict[str, Any],
    traits: pd.DataFrame,
    env_taxon: pd.DataFrame,
    minimum_trait_observations: int,
) -> pd.DataFrame:
    members = unit["member_endpoint_ids"]
    part = traits[traits["endpoint_id"].isin(members)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ]
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value").reset_index()
    wide.columns.name = None
    wide = wide.dropna(subset=members)
    counts = wide.groupby("taxon_name").size().rename("n_trait_observations")
    medians = wide.groupby("taxon_name")[members].median().join(counts)
    medians = medians[medians["n_trait_observations"].ge(minimum_trait_observations)]
    return medians.join(env_taxon, how="inner").reset_index()


def standardized_design(data: pd.DataFrame, predictors: list[str], include_space: bool) -> tuple[np.ndarray, list[str]]:
    blocks = [np.ones(len(data))]
    names = ["intercept"]
    for predictor in predictors:
        blocks.append(standardize(data[predictor].to_numpy(float)))
        names.append(predictor)
    if include_space:
        basis = spherical_basis(data["latitude"].to_numpy(float), data["longitude"].to_numpy(float))
        keep = np.std(basis, axis=0, ddof=0) > 1e-10
        basis = basis[:, keep]
        for index in range(basis.shape[1]):
            blocks.append(standardize(basis[:, index]))
            names.append(f"space_basis_{index + 1}")
    return np.column_stack(blocks), names


def standardized_response(data: pd.DataFrame, members: list[str]) -> np.ndarray:
    return np.column_stack([standardize(data[member].to_numpy(float)) for member in members])


def freedman_lane_coefficient(
    responses: np.ndarray,
    design: np.ndarray,
    coefficient_index: int,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float, np.ndarray]:
    y = np.asarray(responses, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    full_inverse = np.linalg.pinv(design)
    full_beta = full_inverse @ y
    observed_beta = np.asarray(full_beta[coefficient_index], dtype=float).reshape(-1)
    observed = float(np.linalg.norm(observed_beta))
    reduced = np.delete(design, coefficient_index, axis=1)
    reduced_beta = np.linalg.pinv(reduced) @ y
    reduced_fitted = reduced @ reduced_beta
    residual = y - reduced_fitted
    weights = full_inverse[coefficient_index]
    constant = np.asarray(weights @ reduced_fitted, dtype=float).reshape(-1)
    exceed = 0
    remaining = permutations
    batch_size = 128
    n = len(y)
    while remaining:
        size = min(batch_size, remaining)
        order = np.broadcast_to(np.arange(n), (size, n)).copy()
        order = rng.permuted(order, axis=1)
        permuted = residual[order]
        score = np.einsum("bnk,n->bk", permuted, weights, optimize=True)
        simulated = np.linalg.norm(constant[None, :] + score, axis=1)
        exceed += int(np.sum(simulated >= observed - 1e-15))
        remaining -= size
    fitted = design @ full_beta
    return observed_beta, float((exceed + 1) / (permutations + 1)), y - fitted


def fit_base_models(
    units: list[dict[str, Any]],
    traits: pd.DataFrame,
    env_taxon: pd.DataFrame,
    predictors: list[str],
    minimum_trait_observations: int,
    minimum_taxa: int,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    rows: list[dict[str, Any]] = []
    cached: dict[str, pd.DataFrame] = {}
    for unit in units:
        data = unit_taxon_data(unit, traits, env_taxon, minimum_trait_observations)
        cached[unit["unit_id"]] = data
        base = {
            "unit_id": unit["unit_id"],
            "member_endpoint_ids": "|".join(unit["member_endpoint_ids"]),
            "inferential_unit": unit["inferential_unit"],
            "module": unit["module"],
            "analysis_tier": unit["analysis_tier"],
            "n_taxa": int(len(data)),
        }
        if len(data) < max(minimum_taxa, len(predictors) + 5):
            for predictor in predictors:
                rows.append({**base, "predictor": predictor, "status": "insufficient_taxa"})
            continue
        design, names = standardized_design(data, predictors, include_space=False)
        responses = standardized_response(data, unit["member_endpoint_ids"])
        for predictor in predictors:
            coefficient_index = names.index(predictor)
            rng = stable_rng(seed, "full46276_multivariable", unit["unit_id"], predictor)
            beta, p_value, _ = freedman_lane_coefficient(
                responses, design, coefficient_index, permutations, rng
            )
            payload: dict[str, Any] = {
                **base,
                "predictor": predictor,
                "status": "ok",
                "p_value": p_value,
                "model_predictor_count": len(predictors),
            }
            if unit["inferential_unit"] == "linear_endpoint":
                payload["beta_std"] = float(beta[0])
            else:
                payload["beta_sine_std"] = float(beta[0])
                payload["beta_cosine_std"] = float(beta[1])
                payload["effect_magnitude"] = float(np.linalg.norm(beta))
                payload["effect_direction_degrees"] = float(math.degrees(math.atan2(beta[0], beta[1])) % 360.0)
            rows.append(payload)
    result = pd.DataFrame(rows)
    ok = result["status"].eq("ok")
    result["q_fdr_bh_global_family"] = np.nan
    result.loc[ok, "q_fdr_bh_global_family"] = bh_adjust(result.loc[ok, "p_value"])
    return result, cached


def fit_spatial_sensitivity(
    base: pd.DataFrame,
    units_by_id: dict[str, dict[str, Any]],
    cached: dict[str, pd.DataFrame],
    predictors: list[str],
    permutations: int,
    moran_permutations: int,
    seed: int,
) -> pd.DataFrame:
    selected = base[
        base["status"].eq("ok") & pd.to_numeric(base["q_fdr_bh_global_family"], errors="coerce").lt(0.05)
    ].copy()
    rows: list[dict[str, Any]] = []
    for _, base_row in selected.iterrows():
        unit = units_by_id[base_row["unit_id"]]
        data = cached[unit["unit_id"]].dropna(subset=["latitude", "longitude", *predictors]).copy()
        design, names = standardized_design(data, predictors, include_space=True)
        responses = standardized_response(data, unit["member_endpoint_ids"])
        coefficient_index = names.index(base_row["predictor"])
        rng = stable_rng(seed, "full46276_multivariable_spatial", unit["unit_id"], base_row["predictor"])
        beta, p_value, residual = freedman_lane_coefficient(
            responses, design, coefficient_index, permutations, rng
        )
        residual_for_moran = residual[:, 0] if residual.shape[1] == 1 else np.linalg.norm(residual, axis=1)
        moran_i, moran_p, moran_n = moran_test(
            residual_for_moran,
            data["latitude"].to_numpy(float),
            data["longitude"].to_numpy(float),
            moran_permutations,
            5000,
            stable_rng(seed, "full46276_multivariable_moran", unit["unit_id"], base_row["predictor"]),
        )
        direction_ok = True
        payload: dict[str, Any] = {
            "unit_id": unit["unit_id"],
            "member_endpoint_ids": "|".join(unit["member_endpoint_ids"]),
            "inferential_unit": unit["inferential_unit"],
            "module": unit["module"],
            "analysis_tier": unit["analysis_tier"],
            "predictor": base_row["predictor"],
            "n_taxa": int(len(data)),
            "base_q_fdr_bh_global_family": float(base_row["q_fdr_bh_global_family"]),
            "spatial_permutation_p_value": p_value,
            "residual_morans_i": moran_i,
            "residual_morans_p_value": moran_p,
            "residual_morans_n": moran_n,
        }
        if unit["inferential_unit"] == "linear_endpoint":
            payload["base_beta_std"] = float(base_row["beta_std"])
            payload["spatial_beta_std"] = float(beta[0])
            direction_ok = bool(np.sign(beta[0]) == np.sign(float(base_row["beta_std"])))
            payload["same_linear_direction_as_base"] = direction_ok
        else:
            payload["spatial_beta_sine_std"] = float(beta[0])
            payload["spatial_beta_cosine_std"] = float(beta[1])
            payload["spatial_effect_magnitude"] = float(np.linalg.norm(beta))
            payload["same_linear_direction_as_base"] = np.nan
        payload["broad_spatial_sensitivity_pass"] = bool(
            p_value < 0.05 and np.isfinite(moran_p) and moran_p >= 0.05 and direction_ok
        )
        rows.append(payload)
    return pd.DataFrame(rows)


def fit_phylo_at_lambda(
    responses: np.ndarray,
    design: np.ndarray,
    base_covariance: np.ndarray,
    lam: float,
    coefficient_index: int,
) -> dict[str, Any]:
    covariance = lam * base_covariance
    np.fill_diagonal(covariance, np.diag(base_covariance))
    covariance = (covariance + covariance.T) / 2.0
    factor, jitter = cholesky_with_jitter(covariance)
    yw = solve_triangular(factor, responses, lower=True)
    xw = solve_triangular(factor, design, lower=True)
    information = xw.T @ xw
    information_inverse = np.linalg.pinv(information)
    beta = information_inverse @ xw.T @ yw
    residual = yw - xw @ beta
    n = len(yw)
    k = yw.shape[1]
    p = design.shape[1]
    if n <= p + 1:
        raise ValueError("Insufficient residual degrees of freedom")
    sigma_ml = residual.T @ residual / n
    sigma_unbiased = residual.T @ residual / (n - p)
    if k == 1:
        sigma_ml = np.atleast_2d(sigma_ml)
        sigma_unbiased = np.atleast_2d(sigma_unbiased)
    sign, logdet_sigma = np.linalg.slogdet(sigma_ml + np.eye(k) * 1e-12)
    if sign <= 0:
        raise ValueError("Non-positive residual covariance")
    logdet_covariance = 2.0 * float(np.sum(np.log(np.diag(factor))))
    log_likelihood = -0.5 * (k * logdet_covariance + n * logdet_sigma + n * k * (1.0 + math.log(2.0 * math.pi)))
    target_beta = np.asarray(beta[coefficient_index], dtype=float).reshape(-1)
    target_cov = float(information_inverse[coefficient_index, coefficient_index]) * sigma_unbiased
    if k == 1:
        se = math.sqrt(max(float(target_cov[0, 0]), 0.0))
        statistic = float(target_beta[0] / se)
        p_value = float(2.0 * student_t.sf(abs(statistic), df=n - p))
    else:
        statistic = float(target_beta.T @ np.linalg.pinv(target_cov) @ target_beta)
        p_value = float(chi2.sf(statistic, df=k))
    return {
        "lambda": float(lam),
        "log_likelihood": float(log_likelihood),
        "beta": target_beta,
        "p_value": p_value,
        "covariance_jitter": float(jitter),
    }


def fit_phylo_multivariable(
    responses: np.ndarray,
    design: np.ndarray,
    covariance: np.ndarray,
    coefficient_index: int,
) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    for lam in np.linspace(0.0, 1.0, 11):
        try:
            candidates.append(fit_phylo_at_lambda(responses, design, covariance, float(lam), coefficient_index))
        except Exception:
            pass
    def objective(lam: float) -> float:
        try:
            return -fit_phylo_at_lambda(responses, design, covariance, lam, coefficient_index)["log_likelihood"]
        except Exception:
            return float("inf")
    optimized = minimize_scalar(objective, bounds=(0.0, 1.0), method="bounded", options={"xatol": 1e-5})
    if optimized.success and np.isfinite(optimized.fun):
        candidates.append(fit_phylo_at_lambda(responses, design, covariance, float(optimized.x), coefficient_index))
    if not candidates:
        raise ValueError("No phylogenetic fit succeeded")
    return max(candidates, key=lambda item: item["log_likelihood"])


def fit_historical_sensitivity(
    spatial: pd.DataFrame,
    units_by_id: dict[str, dict[str, Any]],
    cached: dict[str, pd.DataFrame],
    predictors: list[str],
    tree_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = spatial[spatial["broad_spatial_sensitivity_pass"].fillna(False).astype(bool)].copy()
    trees = load_trees(tree_dir)
    tree_resources: list[tuple[str, int | None, list[str], np.ndarray]] = []
    for scenario, replicate, tree in trees:
        names, covariance = tree_covariance(tree)
        tree_resources.append((scenario, replicate, names, covariance))
    model_rows: list[dict[str, Any]] = []
    for _, pair in selected.iterrows():
        unit = units_by_id[pair["unit_id"]]
        data = cached[unit["unit_id"]].copy()
        data["taxon_key"] = data["taxon_name"].astype(str).str.strip().str.replace(" ", "_", regex=False)
        data = data.set_index("taxon_key", drop=False)
        for scenario, replicate, names, covariance in tree_resources:
            lookup = {name: index for index, name in enumerate(names)}
            common = [name for name in names if name in data.index]
            base = {
                "unit_id": unit["unit_id"],
                "member_endpoint_ids": "|".join(unit["member_endpoint_ids"]),
                "inferential_unit": unit["inferential_unit"],
                "module": unit["module"],
                "analysis_tier": unit["analysis_tier"],
                "predictor": pair["predictor"],
                "scenario": scenario,
                "replicate": replicate,
                "n_taxa": len(common),
            }
            if len(common) < max(30, len(predictors) + 5):
                model_rows.append({**base, "status": "insufficient_taxa"})
                continue
            ordered = data.loc[common].copy()
            design, design_names = standardized_design(ordered, predictors, include_space=False)
            responses = standardized_response(ordered, unit["member_endpoint_ids"])
            positions = [lookup[name] for name in common]
            phylo_covariance = covariance[np.ix_(positions, positions)]
            coefficient_index = design_names.index(pair["predictor"])
            try:
                fit = fit_phylo_multivariable(responses, design, phylo_covariance, coefficient_index)
                beta = fit.pop("beta")
                payload: dict[str, Any] = {**base, "status": "ok", **fit}
                if unit["inferential_unit"] == "linear_endpoint":
                    payload["pgls_beta_std"] = float(beta[0])
                    payload["same_linear_direction_as_spatial"] = bool(
                        np.sign(beta[0]) == np.sign(float(pair["spatial_beta_std"]))
                    )
                else:
                    payload["pgls_beta_sine_std"] = float(beta[0])
                    payload["pgls_beta_cosine_std"] = float(beta[1])
                    payload["pgls_effect_magnitude"] = float(np.linalg.norm(beta))
                    payload["same_linear_direction_as_spatial"] = np.nan
                model_rows.append(payload)
            except Exception as error:
                model_rows.append({**base, "status": f"failed:{type(error).__name__}:{error}"})
    models = pd.DataFrame(model_rows)
    summary_rows: list[dict[str, Any]] = []
    if not models.empty:
        keys = ["unit_id", "member_endpoint_ids", "inferential_unit", "module", "analysis_tier", "predictor"]
        for group_key, part in models.groupby(keys, sort=True):
            ok = part[part["status"].eq("ok")].copy()
            linear = group_key[2] == "linear_endpoint"
            direction_ok = (
                bool(ok["same_linear_direction_as_spatial"].fillna(False).astype(bool).all())
                if linear and not ok.empty else not linear
            )
            placement_pass = bool(len(ok) == 52 and ok["p_value"].lt(0.05).all() and direction_ok)
            summary_rows.append({
                "unit_id": group_key[0],
                "member_endpoint_ids": group_key[1],
                "inferential_unit": group_key[2],
                "module": group_key[3],
                "analysis_tier": group_key[4],
                "predictor": group_key[5],
                "n_successful_placement_trees": int(len(ok)),
                "n_placement_trees_p_lt_0_05": int(ok["p_value"].lt(0.05).sum()) if not ok.empty else 0,
                "maximum_p_value": float(ok["p_value"].max()) if not ok.empty else np.nan,
                "lambda_min": float(ok["lambda"].min()) if not ok.empty else np.nan,
                "lambda_max": float(ok["lambda"].max()) if not ok.empty else np.nan,
                "historical_placement_sensitivity_pass": placement_pass,
            })
    return models, pd.DataFrame(summary_rows)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    contract = pd.read_csv(args.endpoint_contract, dtype=str, keep_default_na=False)
    if len(contract) != 27:
        raise ValueError(f"Expected 27 endpoint contract rows, found {len(contract)}")
    predictors = list(args.predictors)
    if len(predictors) != 8 or len(set(predictors)) != 8:
        raise ValueError(f"Expected exactly 8 VIFstep-retained predictors, found {predictors}")
    environment = load_environment(args.environment, predictors)
    traits, recovered_counts = load_and_materialize_traits(
        args.traits_long, args.recovered_display, contract, set(environment["obs_id"])
    )
    # Fill taxon labels for any repaired endpoint row that could not inherit one from the original long table.
    taxon_lookup = environment.set_index("obs_id")["taxon_name"]
    missing = traits["taxon_name"].isna()
    if missing.any():
        traits.loc[missing, "taxon_name"] = traits.loc[missing, "obs_id"].map(taxon_lookup)
    units = inferential_units(contract)
    if len(units) != 26:
        raise ValueError(f"Expected 26 inferential units from 27 endpoints, found {len(units)}")
    units_by_id = {unit["unit_id"]: unit for unit in units}
    env_taxon = taxon_environment(environment, predictors)

    base, cached = fit_base_models(
        units, traits, env_taxon, predictors,
        args.minimum_trait_observations_per_taxon, args.minimum_taxa,
        args.permutations, args.seed,
    )
    base.to_csv(args.out_dir / "full46276_multivariable_base.csv", index=False, encoding="utf-8-sig")

    spatial = fit_spatial_sensitivity(
        base, units_by_id, cached, predictors,
        args.permutations, args.moran_permutations, args.seed,
    )
    spatial.to_csv(args.out_dir / "full46276_multivariable_spatial.csv", index=False, encoding="utf-8-sig")

    phylo_models, phylo_summary = fit_historical_sensitivity(
        spatial, units_by_id, cached, predictors, args.tree_dir
    )
    phylo_models.to_csv(args.out_dir / "full46276_multivariable_phylo_models.csv", index=False, encoding="utf-8-sig")
    phylo_summary.to_csv(args.out_dir / "full46276_multivariable_phylo_summary.csv", index=False, encoding="utf-8-sig")

    base_sig = base[base["status"].eq("ok") & pd.to_numeric(base["q_fdr_bh_global_family"], errors="coerce").lt(0.05)]
    spatial_pass = spatial[spatial.get("broad_spatial_sensitivity_pass", pd.Series(dtype=bool)).fillna(False).astype(bool)] if not spatial.empty else spatial
    phylo_pass = phylo_summary[phylo_summary.get("historical_placement_sensitivity_pass", pd.Series(dtype=bool)).fillna(False).astype(bool)] if not phylo_summary.empty else phylo_summary
    report = {
        "analysis_id": "ch1_full46276_vifstep10_multivariable_species_level_v1",
        "status": "complete",
        "estimand": "species-level taxon-median trait/environment association using one simultaneous 8-predictor model per trait unit",
        "source_observations": 46276,
        "environment_taxa": int(environment["taxon_name"].nunique()),
        "registered_endpoints": 27,
        "inferential_units": 26,
        "recovered_endpoint_rows": recovered_counts,
        "predictors": predictors,
        "minimum_trait_observations_per_taxon": args.minimum_trait_observations_per_taxon,
        "permutations": args.permutations,
        "base_fdr_signals": int(len(base_sig)),
        "spatial_passes": int(len(spatial_pass)),
        "historical_pairs_entered": int(len(phylo_summary)),
        "historical_passes": int(len(phylo_pass)),
        "historical_pass_pairs": phylo_pass[[
            "unit_id", "predictor", "module", "analysis_tier",
            "n_successful_placement_trees", "n_placement_trees_p_lt_0_05",
            "maximum_p_value", "lambda_min", "lambda_max",
        ]].to_dict("records") if not phylo_pass.empty else [],
        "claim_boundary": [
            "The 46,276 observations are the frozen spatially thinned measurement base; coefficients are species-level taxon-median comparisons.",
            "All retained environmental predictors enter the same model simultaneously; coefficients are mutually adjusted conditional associations, not causal effects.",
            "The spatial stage adds the existing second-order spherical basis and residual Moran diagnostic; it is not an SPDE process model.",
            "The phylogenetic stage repeats the same 8-predictor species-level model as Pagel-lambda PGLS across 52 audited placement trees; the trees are grafted and do not resolve reticulation.",
            "Closed colour-composition endpoints remain descriptive rather than independent biological discoveries.",
        ],
    }
    (args.out_dir / "full46276_multivariable_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
