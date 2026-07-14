#!/usr/bin/env python3
"""Submission extensions for Ch1: trait modules and environmental niche contrasts.

Uses the frozen enriched observation table and grouped SPDE-INLA summaries. The
analysis is descriptive/inferential support for the manuscript, not a replacement
for the spatial models.
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
ENV = [
    "chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
    "topo_elevation", "topo_slope", "topo_roughness",
    "soil_bdod_0_30cm", "soil_cec_0_30cm", "soil_cfvo_0_30cm",
    "soil_clay_0_30cm", "soil_sand_0_30cm", "soil_silt_0_30cm",
    "soil_nitrogen_0_30cm", "soil_phh2o_0_30cm", "soil_soc_0_30cm",
    "soil_ocd_0_30cm",
]


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--environment", required=True)
    p.add_argument("--spde-fixed", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-species-n", type=int, default=5)
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
    term1 = float((diff.T @ inv @ diff) / 8)
    _, ld = np.linalg.slogdet(c)
    _, ld1 = np.linalg.slogdet(c1)
    _, ld2 = np.linalg.slogdet(c2)
    db = term1 + 0.5 * (ld - 0.5 * (ld1 + ld2))
    return float(np.exp(-max(db, 0.0)))


def main() -> None:
    a = args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.environment, low_memory=False)
    fx = pd.read_csv(a.spde_fixed, low_memory=False)

    traits = [t for ts in TRAIT_MODULES.values() for t in ts]
    env = [e for e in ENV if e in df.columns]
    needed = {"taxon_name", *traits, *env}
    missing = sorted(needed.difference(df.columns))
    if missing:
        raise SystemExit(f"Missing columns: {missing}")

    # Species-level trait architecture.
    counts = df.groupby("taxon_name").size()
    keep = counts[counts >= a.min_species_n].index
    sp = df[df.taxon_name.isin(keep)].groupby("taxon_name")[traits + env].median(numeric_only=True)
    trait_cc = sp[traits].dropna()
    z = StandardScaler().fit_transform(trait_cc)
    pca = PCA().fit(z)
    load = pd.DataFrame(pca.components_.T, index=traits,
                        columns=[f"PC{i+1}" for i in range(pca.n_components_)])
    load.insert(0, "trait", load.index)
    load.to_csv(out / "species_trait_pca_loadings.csv", index=False)
    pd.DataFrame({
        "component": [f"PC{i+1}" for i in range(pca.n_components_)],
        "variance_explained": pca.explained_variance_ratio_,
        "cumulative_variance": np.cumsum(pca.explained_variance_ratio_),
    }).to_csv(out / "species_trait_pca_variance.csv", index=False)

    # Environmental responsiveness by predeclared trait module.
    fx = fx[fx["term"].ne("Intercept")].copy()
    trait_to_module = {t: m for m, ts in TRAIT_MODULES.items() for t in ts}
    fx["module"] = fx["trait"].map(trait_to_module)
    fx["abs_effect"] = fx["mean"].abs()
    module = fx.groupby(["module", "model_group"], dropna=False).agg(
        n_effects=("mean", "size"),
        median_abs_effect=("abs_effect", "median"),
        mean_abs_effect=("abs_effect", "mean"),
        n_credible_95=("credible_95_excludes_zero", "sum"),
        n_bh_global_005=("q_bh_posterior_tail_global", lambda x: int((x < 0.05).sum())),
    ).reset_index()
    module.to_csv(out / "trait_module_environment_responsiveness.csv", index=False)

    # Niche contrasts: high vs low species quartiles for each trait in environmental PCA space.
    env_cc = sp[env].dropna()
    common = trait_cc.index.intersection(env_cc.index)
    env_z = StandardScaler().fit_transform(env_cc.loc[common])
    env_pca = PCA(n_components=min(3, env_z.shape[1])).fit_transform(env_z)
    niche_rows = []
    for trait in traits:
        vals = sp.loc[common, trait].dropna()
        if len(vals) < 20:
            continue
        q1, q3 = vals.quantile([0.25, 0.75])
        low_names = vals.index[vals <= q1]
        high_names = vals.index[vals >= q3]
        idx = pd.Index(common)
        low = env_pca[idx.get_indexer(low_names)]
        high = env_pca[idx.get_indexer(high_names)]
        if min(len(low), len(high)) < 5:
            continue
        cdist = float(np.linalg.norm(low.mean(0) - high.mean(0)))
        breadth_low = float(np.linalg.det(safe_cov(low)) ** (1 / low.shape[1]))
        breadth_high = float(np.linalg.det(safe_cov(high)) ** (1 / high.shape[1]))
        niche_rows.append({
            "trait": trait,
            "module": trait_to_module[trait],
            "n_low": len(low), "n_high": len(high),
            "low_threshold": float(q1), "high_threshold": float(q3),
            "environmental_centroid_distance": cdist,
            "niche_breadth_low": breadth_low,
            "niche_breadth_high": breadth_high,
            "gaussian_bhattacharyya_overlap": bhattacharyya_gaussian(low, high),
        })
    pd.DataFrame(niche_rows).to_csv(out / "trait_environmental_niche_contrasts.csv", index=False)

    report = {
        "species_total": int(len(sp)),
        "species_complete_all_traits": int(len(trait_cc)),
        "environment_predictors_used": env,
        "trait_modules": TRAIT_MODULES,
        "interpretation": {
            "module": "Compare standardized effect magnitudes and supported effects among orientation, colour, and shape modules.",
            "niche": "High and low species quartiles are compared in the first three environmental PCA axes; overlap is descriptive, not causal.",
        },
    }
    (out / "trait_module_niche_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
