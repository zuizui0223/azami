#!/usr/bin/env python3
"""Test predeclared environmental process blocks against the 18D capitulum space.

Within-taxon fits use taxon-centred traits and predictors with inverse taxon
sample-size weights so each taxon contributes equal total weight. Among-taxon
fits use taxon medians. Each environmental block receives one multivariate
permutation test at each scale. Missing process blocks are reported unavailable;
they are never replaced after inspecting trait outcomes.
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
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260827)
    return p.parse_args()


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def stable_rng(seed: int, *parts: str) -> np.random.Generator:
    payload = "|".join([str(seed), *parts]).encode()
    return np.random.default_rng(int.from_bytes(hashlib.sha256(payload).digest()[:8], "little"))


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pvalues.to_numpy(float)
    n = len(p)
    if n == 0:
        return pd.Series(index=pvalues.index, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return pd.Series(out, index=pvalues.index)


def complete_trait_table(traits: pd.DataFrame) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    x = traits[traits["analysis_tier"].isin(["primary", "candidate"]) & as_bool(traits["analysis_eligible"])].copy()
    x["value"] = pd.to_numeric(x["value"], errors="coerce")
    measured = x[x["value"].notna()].copy()
    endpoints = sorted(measured["endpoint_id"].unique())
    if len(endpoints) != 18:
        raise ValueError(f"Expected 18 measured inferential endpoints, found {len(endpoints)}")
    modules = measured.drop_duplicates("endpoint_id").set_index("endpoint_id")["module"].to_dict()
    wide = measured.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value")
    return wide[endpoints].dropna().reset_index(), endpoints, modules


def weighted_standardize(a: np.ndarray, weights: np.ndarray) -> np.ndarray:
    w = weights / weights.sum()
    mu = (w[:, None] * a).sum(0)
    var = (w[:, None] * (a - mu) ** 2).sum(0)
    sd = np.sqrt(var)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("Zero/invalid weighted standard deviation")
    return (a - mu) / sd


def standardize(a: np.ndarray) -> np.ndarray:
    sd = a.std(axis=0, ddof=0)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("Zero/invalid standard deviation")
    return (a - a.mean(axis=0)) / sd


def fit_weighted_multivariate(y: np.ndarray, x: np.ndarray, weights: np.ndarray) -> tuple[float, np.ndarray]:
    yz = weighted_standardize(y, weights)
    xz = weighted_standardize(x, weights)
    sw = np.sqrt(weights)[:, None]
    beta = np.linalg.lstsq(xz * sw, yz * sw, rcond=None)[0]
    fitted = xz @ beta
    r2 = float((weights[:, None] * fitted ** 2).sum() / (weights[:, None] * yz ** 2).sum())
    return r2, beta


def fit_unweighted_multivariate(y: np.ndarray, x: np.ndarray) -> tuple[float, np.ndarray]:
    yz = standardize(y)
    xz = standardize(x)
    beta = np.linalg.lstsq(xz, yz, rcond=None)[0]
    fitted = xz @ beta
    r2 = float((fitted ** 2).sum() / (yz ** 2).sum())
    return r2, beta


def within_arrays(table: pd.DataFrame, endpoints: list[str], predictors: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ydf = table[endpoints].astype(float)
    xdf = table[predictors].astype(float)
    y = (ydf - ydf.groupby(table["taxon_name"]).transform("mean")).to_numpy(float)
    x = (xdf - xdf.groupby(table["taxon_name"]).transform("mean")).to_numpy(float)
    n = table.groupby("taxon_name").size()
    weights = 1.0 / table["taxon_name"].map(n).to_numpy(float)
    groups = table["taxon_name"].astype(str).to_numpy()
    return y, x, weights, groups


def among_arrays(table: pd.DataFrame, endpoints: list[str], predictors: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    med = table.groupby("taxon_name")[endpoints + predictors].median().dropna()
    return med[endpoints].to_numpy(float), med[predictors].to_numpy(float), med.index.astype(str).to_numpy()


def permute_within_rows(x: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = x.copy()
    for group in np.unique(groups):
        idx = np.flatnonzero(groups == group)
        out[idx, :] = x[rng.permutation(idx), :]
    return out


def permutation_p_within(y: np.ndarray, x: np.ndarray, weights: np.ndarray, groups: np.ndarray, observed: float, permutations: int, rng: np.random.Generator) -> float:
    exceed = 0
    for _ in range(permutations):
        xp = permute_within_rows(x, groups, rng)
        r2, _ = fit_weighted_multivariate(y, xp, weights)
        if r2 >= observed - 1e-15:
            exceed += 1
    return float((exceed + 1) / (permutations + 1))


def permutation_p_among(y: np.ndarray, x: np.ndarray, observed: float, permutations: int, rng: np.random.Generator) -> float:
    exceed = 0
    for _ in range(permutations):
        xp = x[rng.permutation(len(x)), :]
        r2, _ = fit_unweighted_multivariate(y, xp)
        if r2 >= observed - 1e-15:
            exceed += 1
    return float((exceed + 1) / (permutations + 1))


def flatten_beta(beta: np.ndarray) -> np.ndarray:
    return beta.reshape(-1)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 else float("nan")


def coefficient_rows(scope: str, block: dict[str, Any], scale: str, predictors: list[str], endpoints: list[str], modules: dict[str, str], beta: np.ndarray) -> list[dict[str, Any]]:
    rows = []
    for pi, predictor in enumerate(predictors):
        for ei, endpoint in enumerate(endpoints):
            rows.append({
                "scope": scope, "block_id": block["block_id"], "tier": block["tier"], "scale": scale,
                "predictor": predictor, "endpoint_id": endpoint, "module": modules[endpoint],
                "beta_standardized": float(beta[pi, ei]),
            })
    return rows


def module_energy_rows(scope: str, block: dict[str, Any], scale: str, endpoints: list[str], modules: dict[str, str], beta: np.ndarray) -> list[dict[str, Any]]:
    endpoint_energy = (beta ** 2).sum(axis=0)
    total = endpoint_energy.sum()
    rows = []
    for module in sorted(set(modules[e] for e in endpoints)):
        idx = [i for i, e in enumerate(endpoints) if modules[e] == module]
        share = float(endpoint_energy[idx].sum() / total) if total > 0 else float("nan")
        expected = len(idx) / len(endpoints)
        rows.append({
            "scope": scope, "block_id": block["block_id"], "tier": block["tier"], "scale": scale,
            "module": module, "n_endpoints": len(idx), "coefficient_energy_share": share,
            "endpoint_count_expected_share": expected,
            "energy_enrichment_vs_endpoint_count": share / expected if expected > 0 else float("nan"),
        })
    return rows


def bootstrap_cosine(table: pd.DataFrame, endpoints: list[str], predictors: list[str], n_boot: int, seed: int, scope: str, block_id: str) -> tuple[float, float, int]:
    if n_boot <= 0:
        return float("nan"), float("nan"), 0
    taxa = table["taxon_name"].drop_duplicates().tolist()
    rng = stable_rng(seed, scope, block_id, "cosine_bootstrap")
    vals = []
    for _ in range(n_boot):
        sampled = rng.choice(taxa, len(taxa), replace=True)
        chunks = []
        for j, taxon in enumerate(sampled):
            g = table[table["taxon_name"].eq(taxon)].copy()
            g["taxon_name"] = f"boot_{j:04d}"
            chunks.append(g)
        b = pd.concat(chunks, ignore_index=True)
        try:
            yw, xw, weights, _ = within_arrays(b, endpoints, predictors)
            _, bw = fit_weighted_multivariate(yw, xw, weights)
            ya, xa, _ = among_arrays(b, endpoints, predictors)
            _, ba = fit_unweighted_multivariate(ya, xa)
            vals.append(cosine(flatten_beta(bw), flatten_beta(ba)))
        except Exception:
            continue
    if not vals:
        return float("nan"), float("nan"), 0
    return float(np.quantile(vals, .025)), float(np.quantile(vals, .975)), len(vals)


def main() -> int:
    args = parse_args()
    traits = pd.read_csv(args.traits_long, low_memory=False)
    env = pd.read_csv(args.environment, low_memory=False)
    contract = json.loads(args.block_contract.read_text())
    blocks = contract["blocks"]
    complete, endpoints, modules = complete_trait_table(traits)
    if env["obs_id"].astype(str).duplicated().any():
        raise ValueError("Environment must be unique by obs_id")
    merged = complete.merge(env.drop(columns=["taxon_name"], errors="ignore"), on="obs_id", how="left", validate="one_to_one")
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    tests, coeffs, energies, geometry, targets, availability = [], [], [], [], [], []
    for threshold in args.thresholds:
        counts = merged.groupby("taxon_name").size()
        keep = counts[counts >= threshold].index
        base = merged[merged["taxon_name"].isin(keep)].copy()
        scope = f"complete18_env_min{threshold}"
        for block in blocks:
            predictors = block["predictors"]
            missing = [x for x in predictors if x not in base.columns]
            finite_ok = not missing and base[predictors].apply(pd.to_numeric, errors="coerce").notna().all(axis=1).mean() >= contract.get("minimum_environment_coverage", 0.98)
            availability.append({"scope": scope, "block_id": block["block_id"], "tier": block["tier"], "available": bool(finite_ok), "missing_predictors": "|".join(missing)})
            if not finite_ok:
                continue
            work = base.dropna(subset=predictors).copy()
            yw, xw, weights, groups = within_arrays(work, endpoints, predictors)
            rw, bw = fit_weighted_multivariate(yw, xw, weights)
            pw = permutation_p_within(yw, xw, weights, groups, rw, args.permutations, stable_rng(args.seed, scope, block["block_id"], "within_perm"))
            ya, xa, taxa = among_arrays(work, endpoints, predictors)
            ra, ba = fit_unweighted_multivariate(ya, xa)
            pa = permutation_p_among(ya, xa, ra, args.permutations, stable_rng(args.seed, scope, block["block_id"], "among_perm"))
            for scale, nobs, r2, p in [("within_taxon", len(work), rw, pw), ("among_taxon", len(taxa), ra, pa)]:
                tests.append({
                    "scope": scope, "block_id": block["block_id"], "tier": block["tier"], "scale": scale,
                    "n_observations_or_taxa": nobs, "n_taxa": work["taxon_name"].nunique(),
                    "n_predictors": len(predictors), "n_response_endpoints": len(endpoints),
                    "multivariate_r2": r2, "permutation_p": p, "construct": block["construct"],
                })
                beta = bw if scale == "within_taxon" else ba
                coeffs += coefficient_rows(scope, block, scale, predictors, endpoints, modules, beta)
                energies += module_energy_rows(scope, block, scale, endpoints, modules, beta)
            cos0 = cosine(flatten_beta(bw), flatten_beta(ba))
            low, high, nboot = bootstrap_cosine(work, endpoints, predictors, args.bootstrap, args.seed, scope, block["block_id"])
            geometry.append({
                "scope": scope, "block_id": block["block_id"], "tier": block["tier"],
                "coefficient_matrix_cosine_within_vs_among": cos0,
                "bootstrap_ci95_low": low, "bootstrap_ci95_high": high,
                "bootstrap_successful_replicates": nboot,
            })
            for scale, val in [("within_taxon", rw), ("among_taxon", ra)]:
                targets.append({"target_id": f"environment_block_r2:{block['block_id']}", "scope": scope, "scale": scale, "value": val, "handoff_status": "observational_environment_block_target"})
            targets.append({"target_id": f"environment_block_cross_scale_cosine:{block['block_id']}", "scope": scope, "scale": "within_vs_among", "value": cos0, "handoff_status": "descriptive_effect_geometry_target"})

    testdf = pd.DataFrame(tests)
    if not testdf.empty:
        testdf["q_bh_across_available_blocks"] = np.nan
        for (_, scale), idx in testdf.groupby(["scope", "scale"]).groups.items():
            testdf.loc[idx, "q_bh_across_available_blocks"] = bh_adjust(testdf.loc[idx, "permutation_p"].astype(float))
        testdf["fdr_supported_0_05"] = testdf["q_bh_across_available_blocks"].lt(.05)
    testdf.to_csv(out / "capitulum_environment_block_tests.csv", index=False)
    pd.DataFrame(coeffs).to_csv(out / "capitulum_environment_block_coefficients.csv", index=False)
    pd.DataFrame(energies).to_csv(out / "capitulum_environment_module_energy.csv", index=False)
    pd.DataFrame(geometry).to_csv(out / "capitulum_environment_cross_scale_geometry.csv", index=False)
    pd.DataFrame(targets).to_csv(out / "capitulum_environment_eazami_targets.csv", index=False)
    pd.DataFrame(availability).to_csv(out / "capitulum_environment_block_availability.csv", index=False)
    report = {
        "environment_contract": contract,
        "endpoint_count": len(endpoints),
        "within_rule": "taxon-centred, endpoint-standardized weighted multivariate OLS with equal total taxon weight",
        "among_rule": "taxon-median standardized multivariate OLS",
        "permutation_rule": "joint row-vector shuffle within taxon for within scale; joint row-vector shuffle among taxa for among scale",
        "multiplicity": "BH across all available predeclared blocks separately within each scope and scale",
        "claim_boundary": "multivariate spatial association geometry, not causal function, selection, adaptation or plasticity",
        "permutations": args.permutations,
        "bootstrap": args.bootstrap,
        "seed": args.seed,
    }
    (out / "capitulum_environment_block_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
