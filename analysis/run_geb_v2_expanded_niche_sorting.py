#!/usr/bin/env python3
"""Expanded among-taxon environmental sorting for the GEB v2 trait universe."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import run_geb_v2_within_among_alignment as ALIGN


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--predictors", nargs="+", default=ALIGN.DEFAULT_PREDICTORS)
    p.add_argument("--primary-min-trait-observations", type=int, default=5)
    p.add_argument("--sensitivity-min-trait-observations", type=int, default=2)
    p.add_argument("--minimum-common-taxa", type=int, default=20)
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260827)
    return p.parse_args()


def safe_cov(x: np.ndarray) -> np.ndarray:
    c = np.cov(x, rowvar=False)
    if np.ndim(c) == 0:
        c = np.array([[float(c)]])
    c = np.asarray(c, dtype=float)
    return c + np.eye(c.shape[0]) * 1e-8


def bhattacharyya_gaussian(a: np.ndarray, b: np.ndarray) -> float:
    m1, m2 = a.mean(0), b.mean(0)
    c1, c2 = safe_cov(a), safe_cov(b)
    c = (c1 + c2) / 2.0
    diff = (m1 - m2)[:, None]
    term1 = float(np.asarray((diff.T @ np.linalg.pinv(c) @ diff) / 8.0).item())
    _, ld = np.linalg.slogdet(c)
    _, ld1 = np.linalg.slogdet(c1)
    _, ld2 = np.linalg.slogdet(c2)
    distance = term1 + 0.5 * (float(ld) - 0.5 * (float(ld1) + float(ld2)))
    return float(np.exp(-max(distance, 0.0)))


def environmental_pca(taxon_env: pd.DataFrame, predictors: list[str], names: list[str]) -> pd.DataFrame:
    table = taxon_env.set_index("taxon_name").loc[names, predictors].dropna()
    if len(table) != len(names):
        raise ValueError("Common taxa must have complete environmental predictors")
    x = table.to_numpy(float)
    means = x.mean(axis=0)
    sds = x.std(axis=0, ddof=0)
    if np.any(~np.isfinite(sds)) or np.any(sds <= 0):
        raise ValueError("Environmental predictor has no between-taxon variance")
    z = (x - means) / sds
    u, s, _vt = np.linalg.svd(z, full_matrices=False)
    n_pc = min(3, z.shape[1])
    scores = u[:, :n_pc] * s[:n_pc]
    return pd.DataFrame(scores, index=table.index, columns=[f"PC{i+1}" for i in range(n_pc)])


def linear_metrics(values: np.ndarray, env_scores: np.ndarray) -> tuple[float, float, int, int]:
    q1, q3 = np.quantile(values, [0.25, 0.75])
    low = env_scores[values <= q1]
    high = env_scores[values >= q3]
    if min(len(low), len(high)) < 5:
        raise ValueError("Too few taxa in a trait quartile")
    centroid = float(np.linalg.norm(low.mean(0) - high.mean(0)))
    overlap = bhattacharyya_gaussian(low, high)
    return centroid, overlap, int(len(low)), int(len(high))


def multivariate_r2(y: np.ndarray, x: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    y = y - y.mean(axis=0, keepdims=True)
    ys = y.std(axis=0, ddof=0)
    if np.any(~np.isfinite(ys)) or np.any(ys <= 0):
        raise ValueError("Circular component has no between-taxon variance")
    y = y / ys
    beta = np.linalg.pinv(x) @ y
    fitted = x @ beta
    total = float(np.sum(y * y))
    return float(np.sum(fitted * fitted) / total) if total > 0 else float("nan")


def unit_tables(
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    predictors: list[str],
    min_trait_observations: int,
) -> tuple[list[dict[str, Any]], dict[str, pd.DataFrame], pd.DataFrame]:
    taxon_env = ALIGN.taxon_environment(environment, predictors)
    units = ALIGN.inferential_units(traits)
    tables: dict[str, pd.DataFrame] = {}
    for unit in units:
        if unit["inferential_unit"] == "linear_endpoint":
            table = ALIGN.linear_taxon_table(traits, unit["endpoint_id"], taxon_env, min_trait_observations)
        else:
            table = ALIGN.circular_taxon_table(traits, unit["members"], taxon_env, min_trait_observations)
        tables[unit["endpoint_id"]] = table
    return units, tables, taxon_env


def run_scope(
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    predictors: list[str],
    min_trait_observations: int,
    minimum_common_taxa: int,
    permutations: int,
    seed: int,
    scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    units, tables, taxon_env = unit_tables(traits, environment, predictors, min_trait_observations)
    name_sets = [set(table["taxon_name"].astype(str)) for table in tables.values()]
    common = sorted(set.intersection(*name_sets)) if name_sets else []
    if len(common) < minimum_common_taxa:
        raise ValueError(f"Too few common taxa for expanded niche sorting: {len(common)}")
    scores = environmental_pca(taxon_env, predictors, common)
    env_values = scores.to_numpy(float)

    linear_rows: list[dict[str, Any]] = []
    hue_rows: list[dict[str, Any]] = []
    for unit in units:
        table = tables[unit["endpoint_id"]].set_index("taxon_name").loc[common]
        if unit["inferential_unit"] == "linear_endpoint":
            values = table["trait_value"].to_numpy(float)
            observed_centroid, observed_overlap, n_low, n_high = linear_metrics(values, env_values)
            null_centroid = np.empty(permutations, dtype=float)
            null_overlap = np.empty(permutations, dtype=float)
            rng = ALIGN.stable_rng(seed, scope, unit["endpoint_id"], "niche")
            for replicate in range(permutations):
                shuffled = rng.permutation(values)
                cdist, overlap, _, _ = linear_metrics(shuffled, env_values)
                null_centroid[replicate] = cdist
                null_overlap[replicate] = overlap
            p_centroid = float((1 + np.sum(null_centroid >= observed_centroid)) / (permutations + 1))
            p_overlap = float((1 + np.sum(null_overlap <= observed_overlap)) / (permutations + 1))
            linear_rows.append({
                "scope": scope,
                "min_trait_observations_per_taxon": min_trait_observations,
                "endpoint_id": unit["endpoint_id"],
                "module": unit["module"],
                "analysis_tier": unit["analysis_tier"],
                "n_common_taxa": len(common),
                "n_low": n_low,
                "n_high": n_high,
                "observed_centroid_distance": observed_centroid,
                "centroid_permutation_p": p_centroid,
                "observed_bhattacharyya_overlap": observed_overlap,
                "overlap_permutation_p": p_overlap,
                "null_centroid_mean": float(null_centroid.mean()),
                "null_overlap_mean": float(null_overlap.mean()),
            })
        else:
            sine = next(x for x in unit["members"] if "sin" in x)
            cosine = next(x for x in unit["members"] if "cos" in x)
            y = table[[sine, cosine]].to_numpy(float)
            observed = multivariate_r2(y, env_values)
            rng = ALIGN.stable_rng(seed, scope, unit["endpoint_id"], "circular_niche")
            exceed = 0
            for _ in range(permutations):
                if multivariate_r2(y[rng.permutation(len(y))], env_values) >= observed - 1e-15:
                    exceed += 1
            hue_rows.append({
                "scope": scope,
                "min_trait_observations_per_taxon": min_trait_observations,
                "endpoint_id": unit["endpoint_id"],
                "module": unit["module"],
                "analysis_tier": unit["analysis_tier"],
                "n_common_taxa": len(common),
                "joint_environmental_r2": observed,
                "permutation_p": float((exceed + 1) / (permutations + 1)),
                "p_value_method": f"paired_hue_taxon_permutation_{permutations}",
            })

    linear = pd.DataFrame(linear_rows)
    for metric in ["centroid", "overlap"]:
        pcol = f"{metric}_permutation_p"
        qcol = f"{metric}_bh_q_within_tier_scope"
        linear[qcol] = np.nan
        for (_, _tier), index in linear.groupby(["scope", "analysis_tier"]).groups.items():
            linear.loc[index, qcol] = ALIGN.bh_adjust(linear.loc[index, pcol].astype(float))
    linear["centroid_fdr_0_05"] = linear["centroid_bh_q_within_tier_scope"].lt(0.05)
    linear["overlap_fdr_0_05"] = linear["overlap_bh_q_within_tier_scope"].lt(0.05)
    linear["both_metrics_fdr_0_05"] = linear["centroid_fdr_0_05"] & linear["overlap_fdr_0_05"]
    linear["any_metric_fdr_0_05"] = linear["centroid_fdr_0_05"] | linear["overlap_fdr_0_05"]
    hue = pd.DataFrame(hue_rows)
    score_out = scores.reset_index()
    score_out.insert(0, "scope", scope)
    report = {
        "scope": scope,
        "minimum_usable_trait_observations_per_taxon_for_every_endpoint": min_trait_observations,
        "n_common_taxa": len(common),
        "n_linear_endpoints": int(len(linear)),
        "n_circular_joint_units": int(len(hue)),
        "n_linear_both_metric_fdr": int(linear["both_metrics_fdr_0_05"].sum()),
        "n_linear_any_metric_fdr": int(linear["any_metric_fdr_0_05"].sum()),
        "environment_predictors": predictors,
        "permutations": permutations,
    }
    return linear, hue, score_out, report


def main() -> int:
    args = parse_args()
    traits, environment = ALIGN.load_inputs(args.traits_long, args.environment, args.predictors)
    scopes = [(f"common_all_endpoints_min{args.primary_min_trait_observations}", args.primary_min_trait_observations)]
    if args.sensitivity_min_trait_observations != args.primary_min_trait_observations:
        scopes.append((f"common_all_endpoints_min{args.sensitivity_min_trait_observations}", args.sensitivity_min_trait_observations))
    linear_frames = []
    hue_frames = []
    score_frames = []
    reports = []
    for scope, threshold in scopes:
        linear, hue, scores, report = run_scope(
            traits, environment, args.predictors, threshold, args.minimum_common_taxa,
            args.permutations, args.seed, scope,
        )
        linear_frames.append(linear)
        hue_frames.append(hue)
        score_frames.append(scores)
        reports.append(report)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(linear_frames, ignore_index=True).to_csv(args.out_dir / "geb_v2_expanded_niche_sorting_linear.csv", index=False)
    pd.concat(hue_frames, ignore_index=True).to_csv(args.out_dir / "geb_v2_expanded_niche_sorting_hue.csv", index=False)
    pd.concat(score_frames, ignore_index=True).to_csv(args.out_dir / "geb_v2_expanded_niche_environment_pca_scores.csv", index=False)
    output = {
        "scopes": reports,
        "multiplicity": "BH separately by primary/candidate tier and by centroid/overlap metric within each scope; hue is one joint circular test per scope",
        "interpretation_boundary": "Among-taxon environmental sorting only; not within-taxon response, phylogenetic correction, adaptation, or mechanism."
    }
    (args.out_dir / "geb_v2_expanded_niche_sorting_report.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
