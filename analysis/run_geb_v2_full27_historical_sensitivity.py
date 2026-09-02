#!/usr/bin/env python3
"""Historical-placement sensitivity for spatially retained GEB v2 pairs.

The audited GBOTB/V.PhyloMaker placement trees are reused as uncertainty
hypotheses.  This script does not claim that those grafted trees resolve the
species history of Cirsium.  Linear endpoints use one-response Pagel-lambda
PGLS; the sine/cosine hue pair uses a joint matrix-normal Pagel-lambda GLS and a
two-degree-of-freedom Wald test.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.linalg import solve_triangular
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, t as student_t


PREDICTORS = [
    "chelsa_bio01",
    "chelsa_bio04",
    "chelsa_bio12",
    "chelsa_bio15",
    "chelsa_rsds_mean",
    "chelsa_vpd_mean",
    "chelsa_sfcwind_mean",
    "chelsa_gsp",
    "chelsa_npp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-contract", required=True, type=Path)
    parser.add_argument("--traits-long", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--spatial-dir", required=True, type=Path)
    parser.add_argument("--tree-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--minimum-taxa", type=int, default=30)
    return parser.parse_args()


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("No finite variation")
    return (values - float(np.mean(values))) / sd


def load_trees(tree_dir: Path) -> list[tuple[str, int | None, Any]]:
    trees: list[tuple[str, int | None, Any]] = [
        ("S1", None, Phylo.read(tree_dir / "gbotb_lcvp_scenario1.tre", "newick")),
        ("S3", None, Phylo.read(tree_dir / "gbotb_lcvp_scenario3.tre", "newick")),
    ]
    random_path = tree_dir / "gbotb_lcvp_scenario2_randomized.trees"
    for index, tree in enumerate(Phylo.parse(random_path, "newick"), start=1):
        trees.append(("S2_random", index, tree))
    if len(trees) != 52:
        raise ValueError(f"Expected 52 placement trees, found {len(trees)}")
    return trees


def tree_covariance(tree: Any) -> tuple[list[str], np.ndarray]:
    terminals = tree.get_terminals()
    names = [terminal.name for terminal in terminals]
    if any(name is None for name in names) or len(set(names)) != len(names):
        raise ValueError("Tree terminal names must be unique and non-empty")
    depths = tree.depths()
    paths: dict[str, list[Any]] = {}
    for terminal in terminals:
        paths[terminal.name] = [tree.root, *tree.get_path(terminal)]
    covariance = np.zeros((len(names), len(names)), dtype=float)
    for i, left in enumerate(names):
        left_path = paths[left]
        covariance[i, i] = float(depths[terminals[i]])
        for j in range(i):
            right_path = paths[names[j]]
            common = tree.root
            for left_node, right_node in zip(left_path, right_path):
                if left_node is not right_node:
                    break
                common = left_node
            value = float(depths.get(common, 0.0))
            covariance[i, j] = value
            covariance[j, i] = value
    positive_diagonal = np.diag(covariance)
    if not np.all(np.isfinite(positive_diagonal)) or np.mean(positive_diagonal) <= 0:
        raise ValueError("Tree covariance has invalid root-to-tip depths")
    covariance /= float(np.mean(positive_diagonal))
    return names, covariance


def cholesky_with_jitter(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    scale = float(np.mean(np.diag(matrix)))
    for multiplier in (0.0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4):
        candidate = matrix.copy()
        jitter = multiplier * scale
        if jitter:
            candidate.flat[:: len(candidate) + 1] += jitter
        try:
            return np.linalg.cholesky(candidate), jitter
        except np.linalg.LinAlgError:
            continue
    raise ValueError("Phylogenetic covariance is not positive definite")


def fit_at_lambda(responses: np.ndarray, design: np.ndarray, base_covariance: np.ndarray, lam: float) -> dict[str, Any]:
    covariance = lam * base_covariance
    np.fill_diagonal(covariance, np.diag(base_covariance))
    covariance = (covariance + covariance.T) / 2.0
    factor, jitter = cholesky_with_jitter(covariance)
    whitened_y = solve_triangular(factor, responses, lower=True)
    whitened_x = solve_triangular(factor, design, lower=True)
    information = whitened_x.T @ whitened_x
    information_inverse = np.linalg.pinv(information)
    beta = information_inverse @ whitened_x.T @ whitened_y
    residual = whitened_y - whitened_x @ beta
    n, responses_n = responses.shape
    parameters_n = design.shape[1]
    if n <= parameters_n + 1:
        raise ValueError("Insufficient residual degrees of freedom")
    sigma_ml = residual.T @ residual / n
    sigma_unbiased = residual.T @ residual / (n - parameters_n)
    sign, logdet_sigma = np.linalg.slogdet(sigma_ml)
    if sign <= 0 or not np.isfinite(logdet_sigma):
        sigma_ml = sigma_ml + np.eye(responses_n) * 1e-10
        sigma_unbiased = sigma_unbiased + np.eye(responses_n) * 1e-10
        sign, logdet_sigma = np.linalg.slogdet(sigma_ml)
    logdet_covariance = 2.0 * float(np.sum(np.log(np.diag(factor))))
    log_likelihood = -0.5 * (
        responses_n * logdet_covariance
        + n * logdet_sigma
        + n * responses_n * (1.0 + math.log(2.0 * math.pi))
    )
    predictor_beta = np.asarray(beta[1], dtype=float).reshape(-1)
    predictor_covariance = float(information_inverse[1, 1]) * sigma_unbiased
    if responses_n == 1:
        standard_error = math.sqrt(max(float(predictor_covariance[0, 0]), 0.0))
        statistic = float(predictor_beta[0] / standard_error)
        p_value = float(2.0 * student_t.sf(abs(statistic), df=n - parameters_n))
        joint_statistic = statistic * statistic
    else:
        joint_statistic = float(
            predictor_beta.T @ np.linalg.pinv(predictor_covariance) @ predictor_beta
        )
        p_value = float(chi2.sf(joint_statistic, df=responses_n))
        standard_error = float("nan")
    return {
        "lambda": float(lam),
        "log_likelihood": float(log_likelihood),
        "beta": predictor_beta,
        "standard_error": standard_error,
        "joint_statistic": joint_statistic,
        "p_value": p_value,
        "covariance_jitter": float(jitter),
        "n_taxa": int(n),
    }


def fit_pagel(responses: np.ndarray, predictor: np.ndarray, covariance: np.ndarray) -> dict[str, Any]:
    y = np.asarray(responses, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    y = np.column_stack([standardize(y[:, column]) for column in range(y.shape[1])])
    x = standardize(np.asarray(predictor, dtype=float))
    design = np.column_stack([np.ones(len(x)), x])

    def objective(lam: float) -> float:
        try:
            return -fit_at_lambda(y, design, covariance, lam)["log_likelihood"]
        except Exception:
            return float("inf")

    candidates: list[dict[str, Any]] = []
    for lam in np.linspace(0.0, 1.0, 11):
        try:
            candidates.append(fit_at_lambda(y, design, covariance, float(lam)))
        except Exception:
            pass
    optimized = minimize_scalar(objective, bounds=(0.0, 1.0), method="bounded", options={"xatol": 1e-5})
    if optimized.success and np.isfinite(optimized.fun):
        candidates.append(fit_at_lambda(y, design, covariance, float(optimized.x)))
    if not candidates:
        raise ValueError("No Pagel-lambda fit succeeded")
    return max(candidates, key=lambda fit: fit["log_likelihood"])


def load_observations(traits_path: Path, environment_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    traits = pd.read_csv(
        traits_path,
        usecols=["obs_id", "taxon_name", "endpoint_id", "measurement_available", "value"],
        low_memory=False,
    )
    environment = pd.read_csv(
        environment_path,
        usecols=["obs_id", "taxon_name", *PREDICTORS],
        low_memory=False,
    )
    traits["obs_id"] = traits["obs_id"].astype(str)
    environment["obs_id"] = environment["obs_id"].astype(str)
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    available = traits["measurement_available"].astype(str).str.lower().isin({"true", "1", "yes"})
    return traits[available & traits["value"].notna()].copy(), environment


def pair_taxon_data(row: pd.Series, traits: pd.DataFrame, environment: pd.DataFrame) -> pd.DataFrame:
    members = str(row["member_endpoint_ids"]).split("|")
    part = traits[traits["endpoint_id"].isin(members)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ]
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value").reset_index()
    wide.columns.name = None
    wide = wide.dropna(subset=members)
    counts = wide.groupby("taxon_name").size().rename("n_trait_observations")
    medians = wide.groupby("taxon_name")[members].median().join(counts)
    minimum = int(row["minimum_trait_observations_per_taxon"])
    medians = medians[medians["n_trait_observations"].ge(minimum)]
    env_taxon = environment.groupby("taxon_name")[PREDICTORS].median(numeric_only=True)
    result = medians.join(env_taxon, how="inner").reset_index()
    result["taxon_key"] = result["taxon_name"].str.strip().str.replace(" ", "_", regex=False)
    return result.dropna(subset=[*members, row["predictor"]])


def main() -> int:
    args = parse_args()
    contract = json.loads(args.analysis_contract.read_text(encoding="utf-8"))
    spatial = pd.read_csv(args.spatial_dir / "v2_full27_spatial_among.csv", low_memory=False)
    selected = spatial[
        spatial["status"].eq("ok") & spatial["broad_spatial_sensitivity_pass"].astype(bool)
    ].copy()
    traits, environment = load_observations(args.traits_long, args.environment)
    trees = load_trees(args.tree_dir)
    tree_resources: list[tuple[str, int | None, list[str], np.ndarray]] = []
    for scenario, replicate, tree in trees:
        names, covariance = tree_covariance(tree)
        tree_resources.append((scenario, replicate, names, covariance))

    rows: list[dict[str, Any]] = []
    for _, pair in selected.iterrows():
        members = str(pair["member_endpoint_ids"]).split("|")
        data = pair_taxon_data(pair, traits, environment)
        data = data.set_index("taxon_key", drop=False)
        for scenario, replicate, names, covariance in tree_resources:
            index_lookup = {name: index for index, name in enumerate(names)}
            common = [name for name in names if name in data.index]
            base = {
                "unit_id": pair["unit_id"],
                "member_endpoint_ids": pair["member_endpoint_ids"],
                "inferential_unit": pair["inferential_unit"],
                "predictor": pair["predictor"],
                "environment_block": pair["environment_block"],
                "scenario": scenario,
                "replicate": replicate,
                "n_taxa": len(common),
                "spatial_beta_std": pair.get("spatial_beta_std", np.nan),
            }
            if len(common) < args.minimum_taxa:
                rows.append({**base, "status": "insufficient_taxa"})
                continue
            ordered = data.loc[common]
            positions = [index_lookup[name] for name in common]
            phylogenetic_covariance = covariance[np.ix_(positions, positions)]
            responses = ordered[members].to_numpy(float)
            try:
                fit = fit_pagel(
                    responses,
                    ordered[pair["predictor"]].to_numpy(float),
                    phylogenetic_covariance,
                )
                beta = fit.pop("beta")
                payload: dict[str, Any]
                if len(members) == 1:
                    payload = {
                        "pgls_beta_std": float(beta[0]),
                        "pgls_standard_error": fit.pop("standard_error"),
                        "same_linear_direction_as_spatial": bool(
                            np.sign(beta[0]) == np.sign(float(pair["spatial_beta_std"]))
                        ),
                    }
                else:
                    payload = {
                        "pgls_beta_sine_std": float(beta[0]),
                        "pgls_beta_cosine_std": float(beta[1]),
                        "pgls_effect_magnitude": float(np.linalg.norm(beta)),
                        "pgls_effect_direction_degrees": float(
                            math.degrees(math.atan2(beta[0], beta[1])) % 360
                        ),
                        "same_linear_direction_as_spatial": np.nan,
                    }
                    fit.pop("standard_error")
                rows.append({**base, "status": "ok", **payload, **fit})
            except Exception as error:
                rows.append({**base, "status": f"failed:{type(error).__name__}:{error}"})

    models = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    for keys, part in models.groupby(
        ["unit_id", "member_endpoint_ids", "inferential_unit", "predictor", "environment_block"],
        sort=True,
    ):
        ok = part[part["status"].eq("ok")].copy()
        linear = keys[2] == "linear_endpoint"
        direction_ok = (
            ok["same_linear_direction_as_spatial"].fillna(False).astype(bool).all()
            if linear and not ok.empty
            else not linear
        )
        all_significant = bool(not ok.empty and ok["p_value"].lt(0.05).all())
        placement_pass = bool(len(ok) == 52 and all_significant and direction_ok)
        summary_rows.append(
            {
                "unit_id": keys[0],
                "member_endpoint_ids": keys[1],
                "inferential_unit": keys[2],
                "predictor": keys[3],
                "environment_block": keys[4],
                "n_successful_placement_trees": int(len(ok)),
                "n_placement_trees_p_lt_0_05": int(ok["p_value"].lt(0.05).sum()) if not ok.empty else 0,
                "minimum_p_value": float(ok["p_value"].min()) if not ok.empty else np.nan,
                "maximum_p_value": float(ok["p_value"].max()) if not ok.empty else np.nan,
                "lambda_min": float(ok["lambda"].min()) if not ok.empty else np.nan,
                "lambda_max": float(ok["lambda"].max()) if not ok.empty else np.nan,
                "linear_direction_stable_across_trees": direction_ok if linear else np.nan,
                "historical_placement_sensitivity_pass": placement_pass,
                "candidate_class": (
                    "adaptive_pattern_candidate_under_current_controls"
                    if placement_pass
                    else "spatial_candidate_historical_constraint_not_cleared"
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    models.to_csv(args.out_dir / "v2_full27_historical_placement_models.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(args.out_dir / "v2_full27_historical_placement_summary.csv", index=False, encoding="utf-8-sig")
    report = {
        "analysis_id": "geb_v2_full27_historical_placement_sensitivity_v1",
        "status": "retrospective_selected-row_sensitivity",
        "n_spatial_pairs_entered": int(len(selected)),
        "n_placement_trees": int(len(trees)),
        "n_pairs_passing_all_placement_trees": int(
            summary.get("historical_placement_sensitivity_pass", pd.Series(dtype=bool)).sum()
        ),
        "method": contract["historical_sensitivity"],
        "claim_boundary": "Passing means stability across audited grafting scenarios, not a resolved species-tree correction or proof of adaptation, selection, convergence or mechanism.",
    }
    (args.out_dir / "v2_full27_historical_sensitivity_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
