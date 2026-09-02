#!/usr/bin/env python3
"""Test whether predeclared process variables add 18D structure beyond core four CHELSA predictors.

The estimand is the same complete-18 capitulum phenotype used by the block
analysis.  Reduced models contain the frozen core four predictors.  Full models
add either all five process-extension predictors (one omnibus test) or one of
four predeclared process blocks.  Significance uses Freedman-Lane permutation of
reduced-model response residual rows.  Within-taxon permutations never move a
residual row between taxa; among-taxon permutations shuffle residual rows across
taxa.

This is an observational redundancy/increment test.  It does not identify a
causal environmental mechanism.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--block-contract", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--thresholds", nargs="+", type=int, default=[5, 2])
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260827)
    return p.parse_args()


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def stable_rng(seed: int, *parts: str) -> np.random.Generator:
    payload = "|".join([str(seed), *parts]).encode()
    return np.random.default_rng(int.from_bytes(hashlib.sha256(payload).digest()[:8], "little"))


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    if len(p) == 0:
        return pd.Series(index=values.index, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * len(p) / np.arange(1, len(p) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(len(p), float)
    out[order] = np.clip(adj, 0, 1)
    return pd.Series(out, index=values.index)


def complete_trait_table(traits: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    required = {
        "obs_id", "taxon_name", "endpoint_id", "analysis_tier",
        "analysis_eligible", "value"
    }
    missing = required.difference(traits.columns)
    if missing:
        raise ValueError(f"Missing trait columns: {sorted(missing)}")
    x = traits[
        traits["analysis_tier"].isin(["primary", "candidate"])
        & as_bool(traits["analysis_eligible"])
    ].copy()
    x["value"] = pd.to_numeric(x["value"], errors="coerce")
    x = x[x["value"].notna()]
    endpoints = sorted(x["endpoint_id"].astype(str).unique())
    if len(endpoints) != 18:
        raise ValueError(f"Expected 18 measured inferential endpoints, found {len(endpoints)}")
    if x.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Trait table must be unique by obs_id/endpoint_id")
    wide = x.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value")
    return wide[endpoints].dropna().reset_index(), endpoints


def weighted_standardize(a: np.ndarray, weights: np.ndarray) -> np.ndarray:
    w = weights / weights.sum()
    mu = (w[:, None] * a).sum(axis=0)
    var = (w[:, None] * (a - mu) ** 2).sum(axis=0)
    sd = np.sqrt(var)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("Zero/invalid weighted standard deviation")
    return (a - mu) / sd


def standardize(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, float)
    sd = a.std(axis=0, ddof=0)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("Zero/invalid standard deviation")
    return (a - a.mean(axis=0)) / sd


def fit_wls(y: np.ndarray, x: np.ndarray, weights: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    sw = np.sqrt(weights)[:, None]
    beta = np.linalg.lstsq(x * sw, y * sw, rcond=None)[0]
    fitted = x @ beta
    residual = y - fitted
    sse = float((weights[:, None] * residual ** 2).sum())
    sst = float((weights[:, None] * y ** 2).sum())
    r2 = 1.0 - sse / sst
    return float(r2), fitted, residual


def fit_ols(y: np.ndarray, x: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    fitted = x @ beta
    residual = y - fitted
    sse = float((residual ** 2).sum())
    sst = float((y ** 2).sum())
    r2 = 1.0 - sse / sst
    return float(r2), fitted, residual


def prepare_within(
    table: pd.DataFrame,
    endpoints: list[str],
    core: list[str],
    extension: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ydf = table[endpoints].astype(float)
    xcore = table[core].astype(float)
    xext = table[extension].astype(float)
    groups = table["taxon_name"].astype(str).to_numpy()
    y = (ydf - ydf.groupby(table["taxon_name"]).transform("mean")).to_numpy(float)
    xc = (xcore - xcore.groupby(table["taxon_name"]).transform("mean")).to_numpy(float)
    xe = (xext - xext.groupby(table["taxon_name"]).transform("mean")).to_numpy(float)
    n = table.groupby("taxon_name").size()
    weights = 1.0 / table["taxon_name"].map(n).to_numpy(float)
    yz = weighted_standardize(y, weights)
    xcz = weighted_standardize(xc, weights)
    xez = weighted_standardize(xe, weights)
    return yz, xcz, xez, weights, groups


def prepare_among(
    table: pd.DataFrame,
    endpoints: list[str],
    core: list[str],
    extension: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    med = table.groupby("taxon_name")[endpoints + core + extension].median().dropna()
    y = standardize(med[endpoints].to_numpy(float))
    xc = standardize(med[core].to_numpy(float))
    xe = standardize(med[extension].to_numpy(float))
    return y, xc, xe, med.index.astype(str).to_numpy()


def nested_effect_wls(y: np.ndarray, xc: np.ndarray, xe: np.ndarray, weights: np.ndarray) -> tuple[float, float, float, float, np.ndarray, np.ndarray]:
    r2_core, fitted_core, residual_core = fit_wls(y, xc, weights)
    full = np.column_stack([xc, xe])
    r2_full, _, _ = fit_wls(y, full, weights)
    delta = max(0.0, r2_full - r2_core)
    partial = delta / max(1e-15, 1.0 - r2_core)
    return r2_core, r2_full, delta, partial, fitted_core, residual_core


def nested_effect_ols(y: np.ndarray, xc: np.ndarray, xe: np.ndarray) -> tuple[float, float, float, float, np.ndarray, np.ndarray]:
    r2_core, fitted_core, residual_core = fit_ols(y, xc)
    full = np.column_stack([xc, xe])
    r2_full, _, _ = fit_ols(y, full)
    delta = max(0.0, r2_full - r2_core)
    partial = delta / max(1e-15, 1.0 - r2_core)
    return r2_core, r2_full, delta, partial, fitted_core, residual_core


def permute_rows_within(values: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = values.copy()
    for group in np.unique(groups):
        idx = np.flatnonzero(groups == group)
        out[idx] = values[rng.permutation(idx)]
    return out


def freedman_lane_within(
    y: np.ndarray,
    xc: np.ndarray,
    xe: np.ndarray,
    weights: np.ndarray,
    groups: np.ndarray,
    observed_delta: float,
    fitted_core: np.ndarray,
    residual_core: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> float:
    exceed = 0
    full = np.column_stack([xc, xe])
    for _ in range(permutations):
        yp = fitted_core + permute_rows_within(residual_core, groups, rng)
        r2c, _, _ = fit_wls(yp, xc, weights)
        r2f, _, _ = fit_wls(yp, full, weights)
        if (r2f - r2c) >= observed_delta - 1e-15:
            exceed += 1
    return float((exceed + 1) / (permutations + 1))


def freedman_lane_among(
    y: np.ndarray,
    xc: np.ndarray,
    xe: np.ndarray,
    observed_delta: float,
    fitted_core: np.ndarray,
    residual_core: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> float:
    exceed = 0
    full = np.column_stack([xc, xe])
    for _ in range(permutations):
        yp = fitted_core + residual_core[rng.permutation(len(residual_core))]
        r2c, _, _ = fit_ols(yp, xc)
        r2f, _, _ = fit_ols(yp, full)
        if (r2f - r2c) >= observed_delta - 1e-15:
            exceed += 1
    return float((exceed + 1) / (permutations + 1))


def test_specs(contract: dict[str, Any]) -> list[dict[str, Any]]:
    incremental = contract["incremental_tests"]
    specs = [{
        "test_id": incremental["omnibus"]["test_id"],
        "test_family": "omnibus",
        "block_id": "all_process_extension",
        "extension_predictors": incremental["omnibus"]["extension_predictors"],
    }]
    for spec in incremental["block_specific"]:
        specs.append({
            "test_id": spec["test_id"],
            "test_family": "block_specific",
            "block_id": spec["block_id"],
            "extension_predictors": spec["extension_predictors"],
        })
    return specs


def main() -> int:
    args = parse_args()
    if args.permutations < 99:
        raise ValueError("At least 99 permutations are required")
    traits = pd.read_csv(args.traits_long, low_memory=False)
    env = pd.read_csv(args.environment, low_memory=False)
    if env["obs_id"].astype(str).duplicated().any():
        raise ValueError("Environment must be unique by obs_id")
    contract = json.loads(args.block_contract.read_text())
    core = contract["core_predictors"]
    complete, endpoints = complete_trait_table(traits)
    merged = complete.merge(env.drop(columns=["taxon_name"], errors="ignore"), on="obs_id", how="left", validate="one_to_one")
    specs = test_specs(contract)
    rows: list[dict[str, Any]] = []

    for threshold in args.thresholds:
        counts = merged.groupby("taxon_name").size()
        retained = counts[counts >= threshold].index
        base = merged[merged["taxon_name"].isin(retained)].copy()
        scope = f"complete18_env_min{threshold}"
        for spec in specs:
            ext = spec["extension_predictors"]
            needed = core + ext
            missing = [x for x in needed if x not in base.columns]
            if missing:
                rows.append({
                    "scope": scope, "scale": "unavailable", **spec,
                    "status": "missing_predictors", "missing_predictors": "|".join(missing),
                })
                continue
            work = base.dropna(subset=needed).copy()
            coverage = len(work) / len(base) if len(base) else 0.0
            if coverage < float(contract.get("minimum_environment_coverage", 0.98)):
                rows.append({
                    "scope": scope, "scale": "unavailable", **spec,
                    "status": "coverage_below_threshold", "environment_complete_fraction": coverage,
                })
                continue

            yw, xcw, xew, weights, groups = prepare_within(work, endpoints, core, ext)
            r2c, r2f, delta, partial, fitted, residual = nested_effect_wls(yw, xcw, xew, weights)
            pw = freedman_lane_within(
                yw, xcw, xew, weights, groups, delta, fitted, residual,
                args.permutations,
                stable_rng(args.seed, scope, spec["test_id"], "within"),
            )
            rows.append({
                "scope": scope, "scale": "within_taxon", **spec,
                "status": "ok", "n_observations_or_taxa": len(work),
                "n_taxa": int(work["taxon_name"].nunique()),
                "n_core_predictors": len(core), "n_extension_predictors": len(ext),
                "r2_core4": r2c, "r2_core4_plus_extension": r2f,
                "delta_r2": delta, "partial_r2": partial,
                "permutation_p": pw,
                "permutation_method": "Freedman-Lane reduced-response residual rows permuted within taxon",
            })

            ya, xca, xea, taxa = prepare_among(work, endpoints, core, ext)
            r2c, r2f, delta, partial, fitted, residual = nested_effect_ols(ya, xca, xea)
            pa = freedman_lane_among(
                ya, xca, xea, delta, fitted, residual,
                args.permutations,
                stable_rng(args.seed, scope, spec["test_id"], "among"),
            )
            rows.append({
                "scope": scope, "scale": "among_taxon", **spec,
                "status": "ok", "n_observations_or_taxa": len(taxa),
                "n_taxa": len(taxa),
                "n_core_predictors": len(core), "n_extension_predictors": len(ext),
                "r2_core4": r2c, "r2_core4_plus_extension": r2f,
                "delta_r2": delta, "partial_r2": partial,
                "permutation_p": pa,
                "permutation_method": "Freedman-Lane reduced-response residual rows permuted across taxa",
            })

    result = pd.DataFrame(rows)
    result["q_bh_block_specific"] = np.nan
    ok_block = result["status"].eq("ok") & result["test_family"].eq("block_specific")
    for (_scope, _scale), idx in result[ok_block].groupby(["scope", "scale"]).groups.items():
        result.loc[idx, "q_bh_block_specific"] = bh_adjust(result.loc[idx, "permutation_p"].astype(float))
    result["supported_0_05"] = False
    omnibus = result["status"].eq("ok") & result["test_family"].eq("omnibus")
    result.loc[omnibus, "supported_0_05"] = result.loc[omnibus, "permutation_p"].lt(0.05)
    block = result["status"].eq("ok") & result["test_family"].eq("block_specific")
    result.loc[block, "supported_0_05"] = result.loc[block, "q_bh_block_specific"].lt(0.05)

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out / "capitulum_environment_incremental_tests.csv", index=False)
    targets = result[result["status"].eq("ok")][[
        "scope", "scale", "test_id", "test_family", "block_id",
        "delta_r2", "partial_r2", "permutation_p", "q_bh_block_specific", "supported_0_05"
    ]].copy()
    targets.insert(0, "target_id", "environment_incremental:" + targets["test_id"].astype(str))
    targets["handoff_status"] = "observational_incremental_environment_target"
    targets.to_csv(out / "capitulum_environment_incremental_eazami_targets.csv", index=False)
    report = {
        "core_predictors": core,
        "n_response_endpoints": len(endpoints),
        "thresholds": args.thresholds,
        "permutations": args.permutations,
        "tests": [s["test_id"] for s in specs],
        "multiplicity": "omnibus unadjusted as one predeclared test per scope/scale; BH across four block-specific tests per scope/scale",
        "interpretation_boundary": "A supported extension contains spatial information beyond the frozen four predictors for the observational 18D estimand; it does not identify a causal process.",
    }
    (out / "capitulum_environment_incremental_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
