#!/usr/bin/env python3
"""Exploratory v4 functional modules; never modifies the frozen v3 evidence."""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.linalg import lstsq

ENDPOINTS = [
    "orientation_image_vertical_angle", "corolla_lab_lightness", "corolla_lab_chroma",
    "corolla_hue_sin", "corolla_hue_cos", "involucre_length_width_ratio",
    "involucre_apical_taper_ratio", "involucre_basal_taper_ratio",
    "bract_projection_p95", "bract_projection_peak_density"
]
GROUPS = {
    "radiation": ["chelsa_rsds_mean"],
    "thermal": ["chelsa_bio01", "chelsa_bio04"],
    "water_demand": ["chelsa_bio12", "chelsa_bio15", "chelsa_vpd_mean", "chelsa_gsp"],
    "wind": ["chelsa_sfcwind_mean"],
    "productivity": ["chelsa_npp"],
}
VARS = [x for group in GROUPS.values() for x in group]
AXES = {
    "presentation": ["orientation_image_vertical_angle"],
    "colour": ["corolla_lab_lightness", "lab_a_proxy", "lab_b_proxy"],
    "involucre": ["involucre_length_width_ratio", "involucre_apical_taper_ratio",
                  "involucre_basal_taper_ratio", "bract_projection_p95",
                  "bract_projection_peak_density"],
}

def zarr(X):
    X = np.asarray(X, float)
    sd = X.std(axis=0)
    return (X - X.mean(axis=0)) / np.where(sd > 1e-10, sd, 1.)

def fit_r2(X, Y):
    Y = zarr(Y)
    if X.shape[1] == 0:
        return 0.
    X = np.column_stack([np.ones(len(X)), zarr(X)])
    coef = lstsq(X, Y, rcond=None)[0]
    return float(np.clip(1 - np.sum((Y - X @ coef) ** 2) / np.sum(Y ** 2), 0, 1))

def fit_r2_weighted(X, Y, weights):
    """Within-taxon multivariate R2 with equal *total* weight per taxon."""
    if X.shape[1] == 0:
        return 0.
    w = np.asarray(weights, float)
    w = w / w.sum()
    X, Y = np.asarray(X, float), np.asarray(Y, float)
    X = X - np.average(X, weights=w, axis=0)
    Y = Y - np.average(Y, weights=w, axis=0)
    sx = np.sqrt(np.average(X * X, weights=w, axis=0))
    sy = np.sqrt(np.average(Y * Y, weights=w, axis=0))
    X = X / np.where(sx > 1e-10, sx, 1.)
    Y = Y / np.where(sy > 1e-10, sy, 1.)
    sw = np.sqrt(w)[:, None]
    coef = lstsq(X * sw, Y * sw, rcond=None)[0]
    return float(np.clip(1 - np.sum(w[:, None] * (Y - X @ coef) ** 2)
                         / np.sum(w[:, None] * Y ** 2), 0, 1))

def shapley(frame, response, include_spatial=False, weights=None):
    """Descriptive Shapley allocation of *in-sample* multivariate R2; not causal importance."""
    names = list(GROUPS)
    covars = ["latitude", "longitude"] if include_spatial else []
    Y = frame[response].to_numpy(float)
    base = frame[covars].to_numpy(float) if covars else np.empty((len(frame), 0))
    cache = {}
    for k in range(1 << len(names)):
        terms = [p for i, name in enumerate(names) if k & (1 << i) for p in GROUPS[name]]
        arr = np.column_stack([base, frame[terms].to_numpy(float)]) if terms else base
        cache[k] = fit_r2(arr, Y) if weights is None else fit_r2_weighted(arr, Y, weights)
    contrib = {}
    for i, name in enumerate(names):
        score = 0.
        for k in cache:
            if k & (1 << i):
                continue
            n = k.bit_count()
            weight = 1 / (len(names) * math.comb(len(names) - 1, n))
            score += weight * (cache[k | (1 << i)] - cache[k])
        contrib[name] = score
    full = (1 << len(names)) - 1
    return {
        "baseline_r2": cache[0], "full_r2": cache[full],
        "block_shapley_r2": contrib,
        "block_unique_delta_r2": {
            name: cache[full] - cache[full ^ (1 << i)] for i, name in enumerate(names)
        },
    }

def rv(a, b):
    """Ordinary RV. Finite-sample baseline is addressed with matched-n sensitivity."""
    a, b = zarr(a), zarr(b)
    A, B = a.T @ a, b.T @ b
    den = np.sqrt(np.sum(A * A) * np.sum(B * B))
    return float(np.sum((a.T @ b) ** 2) / den) if den > 0 else float("nan")

def construct_cohort(traits, environment, min_taxon=5):
    sub = traits.loc[traits.endpoint_id.isin(ENDPOINTS),
                     ["obs_id", "taxon_name", "endpoint_id", "value"]]
    wide = sub.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id",
                     values="value").reset_index()
    wide = wide.merge(environment[["obs_id", "taxon_name", "latitude", "longitude", *VARS]],
                      on=["obs_id", "taxon_name"], how="inner", validate="one_to_one")
    rad = np.hypot(wide.corolla_hue_sin, wide.corolla_hue_cos)
    # Approximate reconstructed Lab coordinates: medians do not commute with nonlinear
    # transformation. These proxies are NOT calibrated actual flower colours.
    good = rad > 0.1
    wide["lab_a_proxy"] = np.where(good, wide.corolla_lab_chroma *
                                   wide.corolla_hue_cos / rad, np.nan)
    wide["lab_b_proxy"] = np.where(good, wide.corolla_lab_chroma *
                                   wide.corolla_hue_sin / rad, np.nan)
    needed = sorted(set(sum(AXES.values(), [])))
    wide = wide.dropna(subset=needed + VARS + ["latitude", "longitude"]).copy()
    sizes = wide.groupby("taxon_name").size()
    return wide[wide.taxon_name.isin(sizes[sizes >= min_taxon].index)].copy()

def integration(frame, replicates=300, seed=20261020):
    taxa = sorted(frame.taxon_name.unique())
    rng = np.random.default_rng(seed)
    members = sorted(set(sum(AXES.values(), [])))
    taxmed = frame.groupby("taxon_name")[members].median().loc[taxa]
    within = {
        key: frame[cols].to_numpy(float) -
             frame.groupby("taxon_name")[cols].transform("mean").to_numpy(float)
        for key, cols in AXES.items()
    }
    among = {key: taxmed[cols].to_numpy(float) for key, cols in AXES.items()}
    groups = frame.reset_index(drop=True).groupby("taxon_name").indices
    out = []
    for x, y in itertools.combinations(AXES, 2):
        obs_within, obs_among = rv(within[x], within[y]), rv(among[x], among[y])
        diffs = []
        for _ in range(replicates):
            picks = np.array([rng.choice(groups[t]) for t in taxa])
            diffs.append(obs_among - rv(within[x][picks], within[y][picks]))
        out.append({
            "left": x, "right": y, "rv_within_raw": obs_within,
            "rv_among_raw": obs_among,
            "matched_n_diff_median": float(np.median(diffs)),
            "matched_n_diff_ci_low": float(np.quantile(diffs, .025)),
            "matched_n_diff_ci_high": float(np.quantile(diffs, .975)),
            "matched_n_positive_fraction": float(np.mean(np.asarray(diffs) > 0)),
        })
    return out

def fit_univariate_within(frame, drop_taxa=()):
    sub = frame[~frame.taxon_name.isin(drop_taxa)].copy()
    sub["y"] = sub.value - sub.groupby("taxon_name").value.transform("mean")
    sub["x"] = sub.chelsa_rsds_mean - sub.groupby("taxon_name").chelsa_rsds_mean.transform("mean")
    sub = sub[sub.y.notna() & sub.x.notna()]
    x = sub.x.to_numpy() / sub.x.std(ddof=0)
    y = sub.y.to_numpy() / sub.y.std(ddof=0)
    count = sub.groupby("taxon_name").obs_id.transform("size").to_numpy(float)
    return {
        "obs": int(len(sub)), "taxa": int(sub.taxon_name.nunique()),
        "beta_observation_weighted": float(x @ y / (x @ x)),
        "beta_equal_taxon_weight": float((x / count) @ y / ((x / count) @ x)),
    }

def among_chroma_radiation(traits, env, drop_taxa=()):
    """Taxon-median chroma x solar slope, identical inclusion/summary to v3 anchor."""
    c = traits[(traits.endpoint_id == "corolla_lab_chroma") &
               traits.value.notna()].groupby("taxon_name").value.agg(["median", "size"])
    c = c[c["size"] >= 5][["median"]]
    exposures = env.groupby("taxon_name").chelsa_rsds_mean.median()
    common = c.join(exposures).drop(index=list(drop_taxa), errors="ignore").dropna()
    x = zarr(common.chelsa_rsds_mean.to_numpy())
    y = zarr(common["median"].to_numpy())
    return {"n_taxa": int(len(common)), "beta_standardized": float(x @ y / (x @ x))}

def run(traits_path, env_path, outdir):
    outdir.mkdir(exist_ok=True, parents=True)
    traits = pd.read_csv(traits_path, usecols=["obs_id", "taxon_name", "endpoint_id", "value"],
                         dtype={"obs_id": str})
    env = pd.read_csv(env_path, dtype={"obs_id": str})
    if not env.obs_id.is_unique or env.obs_id.nunique() != 46276:
        raise ValueError("Expected frozen 46,276-observation source cohort")
    counts = env.groupby("taxon_name").size().sort_values(ascending=False)
    present = traits[(traits.endpoint_id == "corolla_lab_chroma") & traits.value.notna()][
        ["obs_id", "taxon_name", "value"]]
    thin = present.merge(env[["obs_id", "taxon_name", "chelsa_rsds_mean"]],
                         on=["obs_id", "taxon_name"], validate="one_to_one")
    drop_sets = {
        "all": [], "exclude_vulgare": ["Cirsium vulgare"],
        "exclude_arvense": ["Cirsium arvense"],
        "exclude_top2": list(counts.index[:2]),
        "exclude_top5": list(counts.index[:5]),
    }
    sensitivity = {k: fit_univariate_within(thin, v) for k, v in drop_sets.items()}
    # Unlike raw observation leverage, this tests whether a few widespread taxa
    # control the v3 *between-taxon* colour-radiation headline itself.
    between_sensitivity = {
        k: among_chroma_radiation(traits, env, v)
        for k, v in drop_sets.items()
    }
    cohort = construct_cohort(traits, env)
    tax = cohort.groupby("taxon_name").median(numeric_only=True)
    # Species exposure is the median of ALL environmental observations, as in v3.
    exp = env.groupby("taxon_name")[VARS + ["latitude", "longitude"]].median().loc[tax.index]
    among = tax[list(sum(AXES.values(), []))].join(exp)
    blocks = {
        "among": {name: shapley(among, cols) for name, cols in AXES.items()},
        "among_spatial_latlon": {
            name: shapley(among, cols, include_spatial=True)
            for name, cols in AXES.items()
        },
    }
    # Contrast pooled within-taxon environmental association on the very same
    # 1,777/42 complete cohort. Demean each trait and predictor within taxon.
    # Weight each taxon equally, not each photograph/observation equally.
    centered = cohort.copy()
    allcols = sorted(set(sum(AXES.values(), []) + VARS))
    centered[allcols] = (centered[allcols] -
                         centered.groupby("taxon_name")[allcols].transform("mean"))
    sizes = centered.groupby("taxon_name").obs_id.transform("size").to_numpy(float)
    blocks["within_equal_taxon_weight"] = {
        name: shapley(centered, cols, weights=1 / sizes)
        for name, cols in AXES.items()
    }
    result = {
        "source_env_obs": int(len(env)), "source_taxa": int(len(counts)),
        "top_taxa_counts": counts.head(10).to_dict(),
        "top2_fraction": float(counts.iloc[:2].sum() / len(env)),
        "complete_observations": int(len(cohort)),
        "complete_taxa": int(cohort.taxon_name.nunique()),
        "influence_within_chroma_radiation": sensitivity,
        "influence_among_chroma_radiation": between_sensitivity,
        "environment_blocks": blocks,
        "cross_scale_module_integration": integration(cohort),
        "analysis_scope": "retrospective_exploratory_not_confirmatory",
        "measure_limit": "Lab a/b are approximate recombinations of JPEG median C* and hue, not calibrated reflectance",
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--traits", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.traits, args.environment, args.out_dir), ensure_ascii=False, indent=2))
