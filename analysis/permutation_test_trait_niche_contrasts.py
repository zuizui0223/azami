#!/usr/bin/env python3
"""Permutation/null test for Chapter 1 species-level environmental niche contrasts.

The calculation deliberately mirrors the frozen descriptive implementation in
`ch1_global/v2/84_trait_modules_and_niche.py`: species with >=5 observations are
summarized by medians, all nine primary traits and available environment variables
must be complete, environmental PCA is fixed, and low/high trait quartiles are
compared. The null shuffles trait values among the same species, preserving the
environmental availability and trait distribution.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

TRAIT_MODULES = {
    "orientation": ["orientation_angle_degrees_median"],
    "colour": [
        "corolla_lab_lightness_median",
        "corolla_lab_chroma_median",
        "corolla_hue_sin_median",
        "corolla_hue_cos_median",
    ],
    "shape": [
        "shape_aspect_ratio_median",
        "shape_circularity_median",
        "shape_solidity_median",
        "shape_width_cv_median",
    ],
}
TRAITS = [trait for values in TRAIT_MODULES.values() for trait in values]
ENV = [
    "chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
    "topo_elevation", "topo_slope", "topo_roughness",
    "soil_bdod_0_30cm", "soil_cec_0_30cm", "soil_cfvo_0_30cm",
    "soil_clay_0_30cm", "soil_sand_0_30cm", "soil_silt_0_30cm",
    "soil_nitrogen_0_30cm", "soil_phh2o_0_30cm", "soil_soc_0_30cm",
    "soil_ocd_0_30cm",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--direct-backbone-audit", type=Path, default=None)
    p.add_argument("--min-species-n", type=int, default=5)
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260807)
    return p.parse_args()


def safe_cov(x: np.ndarray) -> np.ndarray:
    c = np.cov(x, rowvar=False)
    if np.ndim(c) == 0:
        c = np.array([[float(c)]])
    return c + np.eye(c.shape[0]) * 1e-8


def bhattacharyya_gaussian(a: np.ndarray, b: np.ndarray) -> float:
    m1, m2 = a.mean(0), b.mean(0)
    c1, c2 = safe_cov(a), safe_cov(b)
    c = (c1 + c2) / 2
    diff = (m1 - m2)[:, None]
    inv = np.linalg.pinv(c)
    term1 = float(np.asarray((diff.T @ inv @ diff) / 8).item())
    _, ld = np.linalg.slogdet(c)
    _, ld1 = np.linalg.slogdet(c1)
    _, ld2 = np.linalg.slogdet(c2)
    db = term1 + 0.5 * (float(ld) - 0.5 * (float(ld1) + float(ld2)))
    return float(np.exp(-max(db, 0.0)))


def bh_adjust(values: pd.Series) -> pd.Series:
    p = np.asarray(values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return pd.Series(out, index=values.index)


def species_table(frame: pd.DataFrame, min_species_n: int) -> tuple[pd.DataFrame, list[str]]:
    env = [column for column in ENV if column in frame.columns]
    required = {"taxon_name", *TRAITS, *env}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    counts = frame.groupby("taxon_name").size()
    keep = counts[counts >= min_species_n].index
    sp = frame[frame.taxon_name.isin(keep)].groupby("taxon_name")[TRAITS + env].median(numeric_only=True)
    complete = sp[TRAITS + env].dropna()
    return complete, env


def fixed_environment_pca(sp: pd.DataFrame, env: list[str]) -> np.ndarray:
    z = StandardScaler().fit_transform(sp[env])
    return PCA(n_components=min(3, len(env))).fit_transform(z)


def one_trait_metrics(values: np.ndarray, env_pca: np.ndarray) -> tuple[float, float, int, int]:
    q1, q3 = np.quantile(values, [0.25, 0.75])
    low = env_pca[values <= q1]
    high = env_pca[values >= q3]
    if min(len(low), len(high)) < 5:
        raise ValueError("Too few species in a trait quartile")
    centroid = float(np.linalg.norm(low.mean(0) - high.mean(0)))
    overlap = bhattacharyya_gaussian(low, high)
    return centroid, overlap, len(low), len(high)


def run_scope(
    sp: pd.DataFrame,
    env: list[str],
    scope: str,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    env_pca = fixed_environment_pca(sp, env)
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    null_rows: list[dict[str, object]] = []
    trait_to_module = {t: m for m, traits in TRAIT_MODULES.items() for t in traits}

    for trait_index, trait in enumerate(TRAITS):
        values = sp[trait].to_numpy(float)
        observed_centroid, observed_overlap, n_low, n_high = one_trait_metrics(values, env_pca)
        null_centroid = np.empty(permutations, dtype=float)
        null_overlap = np.empty(permutations, dtype=float)
        for replicate in range(permutations):
            shuffled = rng.permutation(values)
            cdist, overlap, _, _ = one_trait_metrics(shuffled, env_pca)
            null_centroid[replicate] = cdist
            null_overlap[replicate] = overlap
        p_centroid = (1 + np.sum(null_centroid >= observed_centroid)) / (permutations + 1)
        p_overlap = (1 + np.sum(null_overlap <= observed_overlap)) / (permutations + 1)
        rows.append(
            {
                "scope": scope,
                "trait": trait,
                "module": trait_to_module[trait],
                "n_species": len(sp),
                "n_low": n_low,
                "n_high": n_high,
                "observed_centroid_distance": observed_centroid,
                "centroid_permutation_p": p_centroid,
                "observed_bhattacharyya_overlap": observed_overlap,
                "overlap_permutation_p": p_overlap,
                "null_centroid_mean": float(null_centroid.mean()),
                "null_centroid_q025": float(np.quantile(null_centroid, 0.025)),
                "null_centroid_q975": float(np.quantile(null_centroid, 0.975)),
                "null_overlap_mean": float(null_overlap.mean()),
                "null_overlap_q025": float(np.quantile(null_overlap, 0.025)),
                "null_overlap_q975": float(np.quantile(null_overlap, 0.975)),
            }
        )
        # Store compact quantiles rather than all 10,000 values in the submission bundle.
        null_rows.append(
            {
                "scope": scope,
                "trait": trait,
                "centroid_q001": float(np.quantile(null_centroid, 0.001)),
                "centroid_q01": float(np.quantile(null_centroid, 0.01)),
                "centroid_q05": float(np.quantile(null_centroid, 0.05)),
                "centroid_q50": float(np.quantile(null_centroid, 0.50)),
                "centroid_q95": float(np.quantile(null_centroid, 0.95)),
                "centroid_q99": float(np.quantile(null_centroid, 0.99)),
                "overlap_q001": float(np.quantile(null_overlap, 0.001)),
                "overlap_q01": float(np.quantile(null_overlap, 0.01)),
                "overlap_q05": float(np.quantile(null_overlap, 0.05)),
                "overlap_q50": float(np.quantile(null_overlap, 0.50)),
                "overlap_q95": float(np.quantile(null_overlap, 0.95)),
                "overlap_q99": float(np.quantile(null_overlap, 0.99)),
            }
        )

    result = pd.DataFrame(rows)
    result["centroid_bh_q"] = bh_adjust(result["centroid_permutation_p"])
    result["overlap_bh_q"] = bh_adjust(result["overlap_permutation_p"])
    return result, pd.DataFrame(null_rows)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.environment, low_memory=False)
    sp, env = species_table(frame, args.min_species_n)

    all_result, all_null = run_scope(
        sp,
        env,
        "all_complete_species",
        args.permutations,
        args.seed,
    )
    results = [all_result]
    nulls = [all_null]

    direct_names: set[str] = set()
    if args.direct_backbone_audit is not None:
        audit = pd.read_csv(args.direct_backbone_audit)
        required = {"taxon_name", "directly_present_in_GBOTB_extended_LCVP"}
        missing = required.difference(audit.columns)
        if missing:
            raise ValueError(f"Direct-backbone audit lacks columns: {sorted(missing)}")
        direct = audit["directly_present_in_GBOTB_extended_LCVP"]
        if direct.dtype != bool:
            direct = direct.astype(str).str.lower().isin({"true", "1", "yes"})
        direct_names = set(audit.loc[direct, "taxon_name"].astype(str))
        direct_sp = sp.loc[sp.index.intersection(sorted(direct_names))].copy()
        if len(direct_sp) < 20:
            raise ValueError(f"Too few complete direct-backbone taxa: {len(direct_sp)}")
        direct_result, direct_null = run_scope(
            direct_sp,
            env,
            "direct_backbone_complete_species",
            args.permutations,
            args.seed + 100000,
        )
        results.append(direct_result)
        nulls.append(direct_null)

    result = pd.concat(results, ignore_index=True)
    result.to_csv(args.output_dir / "trait_niche_permutation_results.csv", index=False)
    pd.concat(nulls, ignore_index=True).to_csv(
        args.output_dir / "trait_niche_null_quantiles.csv", index=False
    )

    summary = {
        "n_complete_species": int(len(sp)),
        "n_complete_direct_backbone_species": int(
            len(sp.index.intersection(sorted(direct_names))) if direct_names else 0
        ),
        "environment_predictors": env,
        "n_permutations": args.permutations,
        "seed": args.seed,
        "null_rule": "shuffle each trait among the same complete species; keep environmental PCA and availability fixed",
        "multiplicity": "BH separately across the nine primary traits for centroid-distance and overlap tests within each scope",
        "interpretation": "species-level environmental sorting only; not evidence of within-species response or causal adaptation",
    }
    (args.output_dir / "trait_niche_permutation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
