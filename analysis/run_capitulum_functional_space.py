#!/usr/bin/env python3
"""Quantify multilevel organization of the measured GEB-v2 capitulum trait space.

The analysis is pattern-first. It asks whether registered measurement modules
show stronger internal association than between-module association, and whether
within-taxon and among-taxon association geometry is partly shared. Circular hue
is one inferential unit represented by sine/cosine.

Within-taxon correlations use taxon-centred values with inverse taxon sample-size
weights, so each retained taxon contributes equal total weight. Among-taxon
correlations use taxon medians. Functional interpretation is not inferred here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

HUE_COMPONENTS = ["corolla_hue_sin", "corolla_hue_cos"]
HUE_UNIT = "corolla_hue"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--thresholds", nargs="+", type=int, default=[5, 2])
    p.add_argument("--module-permutations", type=int, default=10000)
    p.add_argument("--matrix-permutations", type=int, default=10000)
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260827)
    return p.parse_args()


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def stable_rng(seed: int, *parts: str) -> np.random.Generator:
    payload = "|".join([str(seed), *parts]).encode()
    derived = int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")
    return np.random.default_rng(derived)


def weighted_corr(frame: pd.DataFrame, weights: np.ndarray) -> pd.DataFrame:
    a = frame.to_numpy(float)
    w = np.asarray(weights, float)
    w = w / w.sum()
    mu = (w[:, None] * a).sum(axis=0)
    z = a - mu
    cov = (w[:, None] * z).T @ z
    sd = np.sqrt(np.diag(cov))
    denom = np.outer(sd, sd)
    corr = np.divide(cov, denom, out=np.zeros_like(cov), where=denom > 0)
    np.fill_diagonal(corr, 1.0)
    return pd.DataFrame(corr, index=frame.columns, columns=frame.columns)


def complete_table(traits: pd.DataFrame) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    required = {"obs_id", "taxon_name", "endpoint_id", "module", "analysis_tier", "analysis_eligible", "value"}
    missing = required.difference(traits.columns)
    if missing:
        raise ValueError(f"Missing trait columns: {sorted(missing)}")
    x = traits[traits["analysis_tier"].isin(["primary", "candidate"]) & as_bool(traits["analysis_eligible"])].copy()
    x["value"] = pd.to_numeric(x["value"], errors="coerce")
    measured = x[x["value"].notna()].copy()
    endpoint_ids = sorted(measured["endpoint_id"].unique())
    if len(endpoint_ids) != 18:
        raise ValueError(f"Expected 18 measured inferential endpoints, found {len(endpoint_ids)}")
    if not set(HUE_COMPONENTS).issubset(endpoint_ids):
        raise ValueError("Circular hue components missing")
    module_map = measured.drop_duplicates("endpoint_id").set_index("endpoint_id")["module"].to_dict()
    wide = measured.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value")
    wide = wide[endpoint_ids].dropna().reset_index()
    return wide, endpoint_ids, module_map


def unit_names(endpoint_ids: list[str]) -> list[str]:
    return sorted([x for x in endpoint_ids if x not in HUE_COMPONENTS] + [HUE_UNIT])


def unit_module(unit: str, module_map: dict[str, str]) -> str:
    return module_map[HUE_COMPONENTS[0]] if unit == HUE_UNIT else module_map[unit]


def hue_multiple_r(corr: pd.DataFrame, other: str) -> float:
    r = corr.loc[other, HUE_COMPONENTS].to_numpy(float)
    rx = corr.loc[HUE_COMPONENTS, HUE_COMPONENTS].to_numpy(float)
    r2 = float(r @ np.linalg.pinv(rx) @ r)
    return float(np.sqrt(max(0.0, min(1.0, r2))))


def unit_strength(corr: pd.DataFrame, left: str, right: str) -> float:
    if left == HUE_UNIT:
        return hue_multiple_r(corr, right)
    if right == HUE_UNIT:
        return hue_multiple_r(corr, left)
    return abs(float(corr.loc[left, right]))


def unit_matrix(corr: pd.DataFrame, units: list[str]) -> pd.DataFrame:
    out = pd.DataFrame(np.eye(len(units)), index=units, columns=units, dtype=float)
    for i, left in enumerate(units):
        for right in units[i + 1 :]:
            value = unit_strength(corr, left, right)
            out.loc[left, right] = value
            out.loc[right, left] = value
    return out


def upper_values(matrix: pd.DataFrame) -> np.ndarray:
    a = matrix.to_numpy(float)
    return a[np.triu_indices_from(a, k=1)]


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ar = pd.Series(a).rank(method="average").to_numpy(float)
    br = pd.Series(b).rank(method="average").to_numpy(float)
    return float(np.corrcoef(ar, br)[0, 1])


def module_contrast(matrix: pd.DataFrame, module_by_unit: dict[str, str]) -> tuple[float, float, float]:
    within: list[float] = []
    between: list[float] = []
    names = list(matrix.index)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            target = within if module_by_unit[left] == module_by_unit[right] else between
            target.append(float(matrix.loc[left, right]))
    mw, mb = float(np.mean(within)), float(np.mean(between))
    return mw - mb, mw, mb


def module_permutation_p(matrix: pd.DataFrame, module_by_unit: dict[str, str], permutations: int, rng: np.random.Generator) -> float:
    observed = module_contrast(matrix, module_by_unit)[0]
    names = list(matrix.index)
    labels = np.array([module_by_unit[x] for x in names], dtype=object)
    exceed = 0
    for _ in range(permutations):
        shuffled = rng.permutation(labels)
        mapping = dict(zip(names, shuffled))
        if module_contrast(matrix, mapping)[0] >= observed - 1e-15:
            exceed += 1
    return float((exceed + 1) / (permutations + 1))


def matrix_permutation_p(within: pd.DataFrame, among: pd.DataFrame, permutations: int, rng: np.random.Generator) -> float:
    observed = spearman(upper_values(within), upper_values(among))
    a = within.to_numpy(float)
    b = among.to_numpy(float)
    idx = np.triu_indices_from(a, 1)
    exceed = 0
    for _ in range(permutations):
        perm = rng.permutation(len(among))
        bp = b[np.ix_(perm, perm)]
        if spearman(a[idx], bp[idx]) >= observed - 1e-15:
            exceed += 1
    return float((exceed + 1) / (permutations + 1))


def matrices_for_table(table: pd.DataFrame, endpoint_ids: list[str], units: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    x = table[endpoint_ids].astype(float)
    centred = x - x.groupby(table["taxon_name"]).transform("mean")
    n_by_taxon = table.groupby("taxon_name").size()
    weights = 1.0 / table["taxon_name"].map(n_by_taxon).to_numpy(float)
    within_corr = weighted_corr(centred, weights)
    medians = table.groupby("taxon_name")[endpoint_ids].median()
    among_corr = medians.corr()
    return within_corr, among_corr, unit_matrix(within_corr, units), unit_matrix(among_corr, units)


def eigenspectrum(corr: pd.DataFrame, scope: str, scale: str) -> pd.DataFrame:
    vals = np.linalg.eigvalsh(corr.to_numpy(float))[::-1]
    vals = np.clip(vals, 0, None)
    frac = vals / vals.sum()
    return pd.DataFrame({
        "scope": scope,
        "scale": scale,
        "component": np.arange(1, len(vals) + 1),
        "eigenvalue": vals,
        "variance_fraction": frac,
        "cumulative_variance": np.cumsum(frac),
    })


def long_matrix(matrix: pd.DataFrame, scope: str, scale: str, matrix_type: str) -> pd.DataFrame:
    rows = []
    names = list(matrix.index)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            rows.append({"scope": scope, "scale": scale, "matrix_type": matrix_type, "left": left, "right": right, "value": matrix.loc[left, right]})
    return pd.DataFrame(rows)


def bootstrap_stats(table: pd.DataFrame, endpoint_ids: list[str], units: list[str], module_by_unit: dict[str, str], n_boot: int, seed: int, scope: str) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float], int]:
    if n_boot <= 0:
        nan = (float("nan"), float("nan"))
        return nan, nan, nan, 0
    taxa = table["taxon_name"].drop_duplicates().tolist()
    rng = stable_rng(seed, scope, "bootstrap")
    within_vals, among_vals, similarity_vals = [], [], []
    for _ in range(n_boot):
        sampled = rng.choice(taxa, size=len(taxa), replace=True)
        chunks = []
        for j, taxon in enumerate(sampled):
            g = table[table["taxon_name"].eq(taxon)].copy()
            g["taxon_name"] = f"boot_{j:04d}"
            chunks.append(g)
        b = pd.concat(chunks, ignore_index=True)
        try:
            _, _, wu, au = matrices_for_table(b, endpoint_ids, units)
            within_vals.append(module_contrast(wu, module_by_unit)[0])
            among_vals.append(module_contrast(au, module_by_unit)[0])
            similarity_vals.append(spearman(upper_values(wu), upper_values(au)))
        except Exception:
            continue
    def ci(v: Iterable[float]) -> tuple[float, float]:
        a = np.asarray(list(v), float)
        return float(np.quantile(a, .025)), float(np.quantile(a, .975))
    return ci(within_vals), ci(among_vals), ci(similarity_vals), len(similarity_vals)


def main() -> int:
    args = parse_args()
    traits = pd.read_csv(args.traits_long, low_memory=False)
    complete, endpoint_ids, module_map = complete_table(traits)
    units = unit_names(endpoint_ids)
    module_by_unit = {u: unit_module(u, module_map) for u in units}
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    summaries, integrations, eigens, endpoint_mats, unit_mats, eazami = [], [], [], [], [], []
    for threshold in args.thresholds:
        counts = complete.groupby("taxon_name").size()
        keep = counts[counts >= threshold].index
        table = complete[complete["taxon_name"].isin(keep)].copy()
        scope = f"complete18_min{threshold}"
        wc, ac, wu, au = matrices_for_table(table, endpoint_ids, units)
        wc0, mw0, mb0 = module_contrast(wu, module_by_unit)
        ac0, ma0, mba0 = module_contrast(au, module_by_unit)
        mpw = module_permutation_p(wu, module_by_unit, args.module_permutations, stable_rng(args.seed, scope, "within_module_perm"))
        mpa = module_permutation_p(au, module_by_unit, args.module_permutations, stable_rng(args.seed, scope, "among_module_perm"))
        similarity = spearman(upper_values(wu), upper_values(au))
        sim_p = matrix_permutation_p(wu, au, args.matrix_permutations, stable_rng(args.seed, scope, "matrix_perm"))
        ciw, cia, cis, n_success = bootstrap_stats(table, endpoint_ids, units, module_by_unit, args.bootstrap, args.seed, scope)
        summaries.append({
            "scope": scope,
            "min_complete_observations_per_taxon": threshold,
            "n_complete_observations": len(table),
            "n_taxa": table["taxon_name"].nunique(),
            "n_endpoints": len(endpoint_ids),
            "n_inferential_units": len(units),
            "within_module_contrast": wc0,
            "within_module_contrast_bootstrap_ci95_low": ciw[0],
            "within_module_contrast_bootstrap_ci95_high": ciw[1],
            "among_module_contrast": ac0,
            "among_module_contrast_bootstrap_ci95_low": cia[0],
            "among_module_contrast_bootstrap_ci95_high": cia[1],
            "cross_scale_strength_matrix_spearman": similarity,
            "cross_scale_spearman_bootstrap_ci95_low": cis[0],
            "cross_scale_spearman_bootstrap_ci95_high": cis[1],
            "cross_scale_trait_label_permutation_p": sim_p,
            "bootstrap_successful_replicates": n_success,
        })
        integrations += [
            {"scope": scope, "scale": "within_taxon", "module_contrast": wc0, "mean_within_module_strength": mw0, "mean_between_module_strength": mb0, "module_label_permutation_p": mpw},
            {"scope": scope, "scale": "among_taxon", "module_contrast": ac0, "mean_within_module_strength": ma0, "mean_between_module_strength": mba0, "module_label_permutation_p": mpa},
        ]
        eigens += [eigenspectrum(wc, scope, "within_taxon"), eigenspectrum(ac, scope, "among_taxon")]
        endpoint_mats += [long_matrix(wc, scope, "within_taxon", "signed_endpoint_correlation"), long_matrix(ac, scope, "among_taxon", "signed_endpoint_correlation")]
        unit_mats += [long_matrix(wu, scope, "within_taxon", "inferential_unit_association_strength"), long_matrix(au, scope, "among_taxon", "inferential_unit_association_strength")]
        for target, scale, value, low, high in [
            ("capitulum_within_module_integration_contrast", "within_taxon", wc0, ciw[0], ciw[1]),
            ("capitulum_among_module_integration_contrast", "among_taxon", ac0, cia[0], cia[1]),
            ("capitulum_cross_scale_association_matrix_similarity", "within_vs_among", similarity, cis[0], cis[1]),
        ]:
            eazami.append({"target_id": target, "scope": scope, "scale": scale, "value": value, "ci95_low": low, "ci95_high": high, "handoff_status": "observational_structure_target"})

    pd.DataFrame(summaries).to_csv(out / "capitulum_space_scope_summary.csv", index=False)
    pd.DataFrame(integrations).to_csv(out / "capitulum_space_module_integration.csv", index=False)
    pd.concat(eigens, ignore_index=True).to_csv(out / "capitulum_space_eigenspectra.csv", index=False)
    pd.concat(endpoint_mats, ignore_index=True).to_csv(out / "capitulum_space_endpoint_correlation_matrices.csv", index=False)
    pd.concat(unit_mats, ignore_index=True).to_csv(out / "capitulum_space_unit_strength_matrices.csv", index=False)
    pd.DataFrame(eazami).to_csv(out / "capitulum_space_eazami_targets.csv", index=False)
    report = {
        "n_complete18_all": int(len(complete)),
        "n_complete18_taxa_all": int(complete["taxon_name"].nunique()),
        "endpoint_ids": endpoint_ids,
        "inferential_units": units,
        "within_rule": "taxon-centred endpoint values with inverse taxon sample-size weights; equal total taxon weight",
        "among_rule": "taxon medians",
        "hue_rule": "joint multiple-correlation strength of sine/cosine components",
        "claim_boundary": "phenotypic/measurement-module organization only; not validated functional or genetic modularity",
        "seed": args.seed,
        "module_permutations": args.module_permutations,
        "matrix_permutations": args.matrix_permutations,
        "bootstrap": args.bootstrap,
    }
    (out / "capitulum_space_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
