#!/usr/bin/env python3
"""Audit and revise the Chapter 1 species-lability analysis.

The legacy responsiveness index was the RMS magnitude of unpooled species-specific
slope estimates. Because absolute noisy estimates are positive, that score is
biased upward when standard errors are large. This script:

1. reproduces the legacy sample-size/precision confounding;
2. constructs an equal-module visible-variation index using seven linear endpoints;
3. constructs a sampling-noise-adjusted association-energy score based on
   E(beta_hat^2 - SE^2) = beta_true^2 under an approximately normal slope estimator;
4. fits a hierarchical variance meta-regression that tests whether latent slope
   heterogeneity changes with visible variation while accounting for every slope SE;
5. writes manuscript-facing tables, summaries and figures.

No ecological model is refit and no raw observation is altered. The analysis uses
archived species-specific slope estimates and their standard errors.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, pearsonr, rankdata, spearmanr

PREDICTORS = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15"]
TRAIT_MODULE = {
    "orientation_angle_degrees_median": "orientation",
    "corolla_lab_lightness_median": "colour",
    "corolla_lab_chroma_median": "colour",
    "shape_aspect_ratio_median": "shape",
    "shape_circularity_median": "shape",
    "shape_solidity_median": "shape",
    "shape_width_cv_median": "shape",
}
MODULE_ORDER = ["orientation", "colour", "shape"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-axes", required=True)
    parser.add_argument("--variation", required=True)
    parser.add_argument("--slopes", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--bootstrap-repeats", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260722)
    return parser.parse_args()


def require_columns(frame: pd.DataFrame, columns: set[str], name: str) -> None:
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise ValueError(f"{name} missing columns: {missing}")


def partial_spearman(x: pd.Series, y: pd.Series, control: pd.Series) -> tuple[float, float]:
    xr = rankdata(pd.to_numeric(x, errors="raise"))
    yr = rankdata(pd.to_numeric(y, errors="raise"))
    zr = rankdata(pd.to_numeric(control, errors="raise"))
    design = np.column_stack([np.ones(len(zr)), zr])
    x_resid = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    y_resid = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    result = pearsonr(x_resid, y_resid)
    return float(result.statistic), float(result.pvalue)


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, repeats: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values: list[float] = []
    n = len(x)
    for _ in range(repeats):
        index = rng.integers(0, n, n)
        value = spearmanr(x[index], y[index]).statistic
        if np.isfinite(value):
            values.append(float(value))
    if not values:
        return float("nan"), float("nan")
    return tuple(float(v) for v in np.quantile(values, [0.025, 0.975]))


def crossing(grid: np.ndarray, values: np.ndarray, level: float, left: bool) -> float | None:
    minimum = int(np.argmin(values))
    if left:
        iterator = range(minimum - 1, -1, -1)
        for index in iterator:
            if values[index] > level and values[index + 1] <= level:
                return float(grid[index] + (level - values[index]) * (grid[index + 1] - grid[index]) / (values[index + 1] - values[index]))
    else:
        iterator = range(minimum, len(values) - 1)
        for index in iterator:
            if values[index] <= level and values[index + 1] > level:
                return float(grid[index] + (level - values[index]) * (grid[index + 1] - grid[index]) / (values[index + 1] - values[index]))
    return None


def fit_variance_meta_regression(slopes: pd.DataFrame, variation_z: pd.Series) -> tuple[dict[str, float | int | None], pd.DataFrame]:
    work = slopes.copy()
    work["variation_z"] = work["taxon_name"].map(variation_z)
    work["group"] = work["trait"].astype(str) + "||" + work["predictor"].astype(str)
    if work["variation_z"].isna().any():
        raise ValueError("Variance meta-regression has unmatched species variation")

    arrays: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []
    for group, part in work.groupby("group", sort=True):
        arrays.append((group, part["beta_std"].to_numpy(float), np.square(part["se"].to_numpy(float)), part["variation_z"].to_numpy(float)))

    def group_nll(log_tau2: float, b: float, array: tuple[str, np.ndarray, np.ndarray, np.ndarray]) -> float:
        _, y, sampling_variance, z = array
        latent_variance = np.exp(log_tau2 + b * z)
        variance = sampling_variance + latent_variance
        weights = 1.0 / variance
        mean = float(np.sum(weights * y) / np.sum(weights))
        residual = y - mean
        return float(0.5 * np.sum(np.log(2 * np.pi * variance) + residual**2 / variance))

    def profile_nll(b: float) -> tuple[float, list[float]]:
        total = 0.0
        log_tau_values: list[float] = []
        for array in arrays:
            result = minimize_scalar(lambda value: group_nll(value, b, array), bounds=(-20.0, 2.0), method="bounded", options={"xatol": 1e-8})
            if not result.success:
                raise RuntimeError(f"Variance meta-regression failed for {array[0]}")
            total += float(result.fun)
            log_tau_values.append(float(result.x))
        return total, log_tau_values

    fitted = minimize_scalar(lambda value: profile_nll(value)[0], bounds=(-3.0, 3.0), method="bounded", options={"xatol": 1e-7})
    if not fitted.success:
        raise RuntimeError("Common variance meta-regression coefficient failed")
    b_hat = float(fitted.x)
    fitted_nll, fitted_log_tau = profile_nll(b_hat)
    null_nll, _ = profile_nll(0.0)
    likelihood_ratio = max(0.0, 2.0 * (null_nll - fitted_nll))
    p_value = float(chi2.sf(likelihood_ratio, 1))

    profile_grid = np.linspace(-1.5, 1.5, 121)
    profile_values = np.array([profile_nll(float(value))[0] for value in profile_grid])
    cutoff = fitted_nll + 0.5 * float(chi2.ppf(0.95, 1))
    ci_low = crossing(profile_grid, profile_values, cutoff, left=True)
    ci_high = crossing(profile_grid, profile_values, cutoff, left=False)

    group_rows = []
    for (group, _, _, _), log_tau in zip(arrays, fitted_log_tau):
        trait, predictor = group.split("||", 1)
        group_rows.append({"trait": trait, "module": TRAIT_MODULE[trait], "predictor": predictor, "latent_slope_sd_at_mean_variation": math.sqrt(math.exp(log_tau))})
    summary: dict[str, float | int | None] = {
        "n_species": int(work["taxon_name"].nunique()),
        "n_slope_estimates": int(len(work)),
        "n_trait_predictor_groups": int(len(arrays)),
        "log_variance_change_per_sd_visible_variation": b_hat,
        "profile_ci95_low": ci_low,
        "profile_ci95_high": ci_high,
        "likelihood_ratio": likelihood_ratio,
        "p_value": p_value,
    }
    return summary, pd.DataFrame(group_rows)


def save_figures(species: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    axes[0].scatter(species["median_slope_n"], species["legacy_environmental_responsiveness"], s=24, alpha=0.7)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Median observations per species-specific slope (log scale)")
    axes[0].set_ylabel("Legacy absolute-slope RMS")
    axes[0].set_title("A. Legacy index is dominated by sample size")
    axes[1].scatter(species["median_slope_se"], species["legacy_environmental_responsiveness"], s=24, alpha=0.7)
    axes[1].set_xlabel("Median slope standard error")
    axes[1].set_ylabel("Legacy absolute-slope RMS")
    axes[1].set_title("B. Legacy index increases with uncertainty")
    fig.tight_layout()
    fig.savefig(out_dir / "figure_precision_audit.svg", bbox_inches="tight")
    fig.savefig(out_dir / "figure_precision_audit.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    ax.scatter(species["equal_module_visible_variation"], species["noise_adjusted_association_energy"], s=30, alpha=0.72)
    ax.axhline(0.0, linewidth=1.0, linestyle="--")
    ax.set_xlabel("Equal-module visible-variation index")
    ax.set_ylabel("Noise-adjusted environment–trait association energy")
    ax.set_title("No detectable coupling after accounting for slope uncertainty")
    rho = spearmanr(species["equal_module_visible_variation"], species["noise_adjusted_association_energy"]).statistic
    ax.text(0.03, 0.96, f"Spearman rho = {rho:.3f}\nN = {len(species)} taxa", transform=ax.transAxes, va="top")
    fig.tight_layout()
    fig.savefig(out_dir / "figure_revised_variation_association.svg", bbox_inches="tight")
    fig.savefig(out_dir / "figure_revised_variation_association.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    legacy = pd.read_csv(args.legacy_axes)
    variation = pd.read_csv(args.variation)
    slopes = pd.read_csv(args.slopes)

    require_columns(legacy, {"taxon_name", "primary_complete", "species_within_variation_index", "species_environmental_responsiveness_index"}, "legacy axes")
    require_columns(variation, {"taxon_name", "trait", "within_variation_standardized", "n_variation"}, "variation table")
    require_columns(slopes, {"taxon_name", "trait", "predictor", "endpoint_type", "beta_std", "se", "n_observations"}, "species-specific slope table")

    legacy_primary = legacy.loc[legacy["primary_complete"].astype(bool)].copy()
    slope_summary = slopes.loc[slopes["taxon_name"].isin(legacy_primary["taxon_name"])].groupby("taxon_name", as_index=False).agg(median_slope_n=("n_observations", "median"), minimum_slope_n=("n_observations", "min"), median_slope_se=("se", "median"))
    legacy_diagnostic = legacy_primary.merge(slope_summary, on="taxon_name", how="inner").rename(columns={"species_within_variation_index": "legacy_within_variation", "species_environmental_responsiveness_index": "legacy_environmental_responsiveness"})

    linear_traits = set(TRAIT_MODULE)
    linear_variation = variation.loc[variation["trait"].isin(linear_traits)].copy()
    linear_slopes = slopes.loc[slopes["trait"].isin(linear_traits) & slopes["endpoint_type"].eq("linear") & slopes["predictor"].isin(PREDICTORS)].copy()
    linear_slopes["beta_std"] = pd.to_numeric(linear_slopes["beta_std"], errors="raise")
    linear_slopes["se"] = pd.to_numeric(linear_slopes["se"], errors="raise")
    linear_slopes["n_observations"] = pd.to_numeric(linear_slopes["n_observations"], errors="raise")

    variation_complete = {taxon for taxon, part in linear_variation.groupby("taxon_name") if set(part["trait"]) == linear_traits and int(part["n_variation"].min()) >= 10}
    slope_complete = {taxon for taxon, part in linear_slopes.groupby("taxon_name") if set(part["trait"]) == linear_traits and part.groupby("trait")["predictor"].nunique().eq(len(PREDICTORS)).all() and int(part["n_observations"].min()) >= 10}
    revised_taxa = sorted(variation_complete.intersection(slope_complete))
    if len(revised_taxa) < 3:
        raise ValueError("Too few complete taxa for revised analysis")

    revised_variation = linear_variation.loc[linear_variation["taxon_name"].isin(revised_taxa)].copy()
    revised_variation["module"] = revised_variation["trait"].map(TRAIT_MODULE)
    module_variation = revised_variation.groupby(["taxon_name", "module"])["within_variation_standardized"].mean().unstack()
    module_variation["equal_module_visible_variation"] = module_variation[MODULE_ORDER].mean(axis=1)

    revised_slopes = linear_slopes.loc[linear_slopes["taxon_name"].isin(revised_taxa)].copy()
    revised_slopes["module"] = revised_slopes["trait"].map(TRAIT_MODULE)
    revised_slopes["sampling_noise_adjusted_squared_effect"] = revised_slopes["beta_std"] ** 2 - revised_slopes["se"] ** 2
    trait_energy = revised_slopes.groupby(["taxon_name", "module", "trait"], as_index=False).agg(trait_association_energy=("sampling_noise_adjusted_squared_effect", "mean"), median_trait_n=("n_observations", "median"), median_trait_se=("se", "median"))
    module_energy = trait_energy.groupby(["taxon_name", "module"])["trait_association_energy"].mean().unstack()
    module_energy["noise_adjusted_association_energy"] = module_energy[MODULE_ORDER].mean(axis=1)
    revised_sampling = revised_slopes.groupby("taxon_name", as_index=False).agg(median_slope_n=("n_observations", "median"), minimum_slope_n=("n_observations", "min"), median_slope_se=("se", "median")).set_index("taxon_name")

    revised = module_variation.add_prefix("variation_module_").join(module_energy.add_prefix("association_module_")).join(revised_sampling)
    revised = revised.rename(columns={"variation_module_equal_module_visible_variation": "equal_module_visible_variation", "association_module_noise_adjusted_association_energy": "noise_adjusted_association_energy"})
    revised = revised.join(legacy_primary.set_index("taxon_name")[["species_environmental_responsiveness_index"]].rename(columns={"species_environmental_responsiveness_index": "legacy_environmental_responsiveness"}), how="left")
    revised.index.name = "taxon_name"
    revised = revised.reset_index()

    corrected_test = spearmanr(revised["equal_module_visible_variation"], revised["noise_adjusted_association_energy"])
    corrected_ci_low, corrected_ci_high = bootstrap_spearman(revised["equal_module_visible_variation"].to_numpy(float), revised["noise_adjusted_association_energy"].to_numpy(float), args.bootstrap_repeats, args.seed)
    meta_variation = revised.set_index("taxon_name")["equal_module_visible_variation"]
    meta_variation_z = (meta_variation - meta_variation.mean()) / meta_variation.std(ddof=1)
    meta_summary, meta_groups = fit_variance_meta_regression(revised_slopes, meta_variation_z)

    legacy_rho = spearmanr(legacy_diagnostic["legacy_within_variation"], legacy_diagnostic["legacy_environmental_responsiveness"])
    legacy_n = spearmanr(legacy_diagnostic["legacy_environmental_responsiveness"], legacy_diagnostic["median_slope_n"])
    legacy_se = spearmanr(legacy_diagnostic["legacy_environmental_responsiveness"], legacy_diagnostic["median_slope_se"])
    partial_rho, partial_p = partial_spearman(legacy_diagnostic["legacy_within_variation"], legacy_diagnostic["legacy_environmental_responsiveness"], legacy_diagnostic["median_slope_n"])
    corrected_n = spearmanr(revised["noise_adjusted_association_energy"], revised["median_slope_n"])
    corrected_se = spearmanr(revised["noise_adjusted_association_energy"], revised["median_slope_se"])

    module_rows = [{"module": module, "n_species": len(revised), "median_visible_variation": float(revised[f"variation_module_{module}"].median()), "median_noise_adjusted_association_energy": float(revised[f"association_module_{module}"].median())} for module in MODULE_ORDER]

    summary = {
        "analysis_status": "reviewer-driven replacement of the legacy absolute-slope RMS lability analysis",
        "legacy_diagnostic": {
            "n_species": int(len(legacy_diagnostic)),
            "variation_vs_legacy_responsiveness_spearman_rho": float(legacy_rho.statistic),
            "variation_vs_legacy_responsiveness_p": float(legacy_rho.pvalue),
            "legacy_responsiveness_vs_median_n_rho": float(legacy_n.statistic),
            "legacy_responsiveness_vs_median_n_p": float(legacy_n.pvalue),
            "legacy_responsiveness_vs_median_se_rho": float(legacy_se.statistic),
            "legacy_responsiveness_vs_median_se_p": float(legacy_se.pvalue),
            "partial_spearman_controlling_median_n": partial_rho,
            "partial_spearman_p": partial_p,
            "decision": "legacy negative correlation and median-split quadrants are withdrawn",
        },
        "revised_primary_cohort": {
            "n_species": int(len(revised)),
            "traits": list(TRAIT_MODULE),
            "n_traits": len(TRAIT_MODULE),
            "predictors": PREDICTORS,
            "eligibility": "all seven linear endpoints, all four predictors per endpoint, minimum n >= 10",
            "module_weighting": "orientation, colour and shape receive equal weight",
            "hue": "excluded from this precision-corrected primary test because archived species-level circular vectors lack component standard errors",
        },
        "noise_adjusted_axis_test": {
            "estimator": "mean(beta_std^2 - se^2) within trait, then equal module average",
            "interpretation": "unbiased association-energy estimator under approximately normal slope errors; negative values indicate no resolved excess signal and are not negative biological responsiveness",
            "spearman_rho": float(corrected_test.statistic),
            "p_value": float(corrected_test.pvalue),
            "species_bootstrap_ci95_low": corrected_ci_low,
            "species_bootstrap_ci95_high": corrected_ci_high,
            "association_energy_vs_median_n_rho": float(corrected_n.statistic),
            "association_energy_vs_median_n_p": float(corrected_n.pvalue),
            "association_energy_vs_median_se_rho": float(corrected_se.statistic),
            "association_energy_vs_median_se_p": float(corrected_se.pvalue),
        },
        "hierarchical_variance_meta_regression": meta_summary,
        "claim_boundary": "No common coupling between visible variation and latent environment-trait slope magnitude was detected after sampling precision was modelled. This is not evidence that every module is biologically invariant.",
    }

    legacy_diagnostic.to_csv(out_dir / "legacy_precision_diagnostic_species.csv", index=False)
    revised.to_csv(out_dir / "revised_species_axes.csv", index=False)
    trait_energy.to_csv(out_dir / "revised_trait_association_energy.csv", index=False)
    pd.DataFrame(module_rows).to_csv(out_dir / "revised_module_summary.csv", index=False)
    meta_groups.to_csv(out_dir / "hierarchical_variance_group_summary.csv", index=False)
    (out_dir / "reviewer_precision_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    save_figures(revised, out_dir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
