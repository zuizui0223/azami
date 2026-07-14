#!/usr/bin/env python3
"""Build two-axis lability summaries and manuscript supplement outputs.

Outputs
-------
- species_lability_axes.csv
- table_s1_all_coefficients.csv
- figure_s1_effect_size_heatmap.png/pdf
- figure_lability_quadrants.png/pdf
- module_environmental_response_summary.csv

The species within-variation index is a robust, trait-standardized median absolute
within-species deviation. The environmental responsiveness index is the root mean
square of standardized within-species slopes across the four predeclared CHELSA
predictors. Circular hue is represented by a joint sine/cosine vector and reported
as an effect magnitude and direction angle; sine and cosine are never interpreted
as separate biological endpoints.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PREDICTORS = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15"]
LINEAR_TRAITS = [
    "orientation_angle_degrees_median",
    "corolla_lab_lightness_median",
    "corolla_lab_chroma_median",
    "shape_aspect_ratio_median",
    "shape_circularity_median",
    "shape_solidity_median",
    "shape_width_cv_median",
]
HUE_SIN = "corolla_hue_sin_median"
HUE_COS = "corolla_hue_cos_median"
MODULE = {
    "orientation_angle_degrees_median": "orientation",
    "corolla_lab_lightness_median": "colour",
    "corolla_lab_chroma_median": "colour",
    "corolla_hue_angle": "colour",
    "shape_aspect_ratio_median": "shape",
    "shape_circularity_median": "shape",
    "shape_solidity_median": "shape",
    "shape_width_cv_median": "shape",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--observations", required=True)
    p.add_argument("--coefficients", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-species-n", type=int, default=5)
    return p.parse_args()


def robust_scale(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    med = x.median()
    mad = (x - med).abs().median()
    scale = 1.4826 * mad
    if not math.isfinite(scale) or scale <= 0:
        scale = x.std(ddof=1)
    if not math.isfinite(scale) or scale <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x - med) / scale


def circular_hue_angle(df: pd.DataFrame) -> pd.Series:
    s = pd.to_numeric(df[HUE_SIN], errors="coerce")
    c = pd.to_numeric(df[HUE_COS], errors="coerce")
    return np.degrees(np.arctan2(s, c)) % 360.0


def circular_distance_deg(x: pd.Series, center: float) -> pd.Series:
    return ((x - center + 180.0) % 360.0 - 180.0).abs()


def species_variation(obs: pd.DataFrame, min_n: int) -> pd.DataFrame:
    obs = obs.copy()
    obs["corolla_hue_angle"] = circular_hue_angle(obs)
    rows = []
    trait_cols = LINEAR_TRAITS + ["corolla_hue_angle"]
    for trait in trait_cols:
        if trait == "corolla_hue_angle":
            values = pd.to_numeric(obs[trait], errors="coerce")
            radians = np.radians(values)
            center = np.degrees(np.arctan2(np.nanmean(np.sin(radians)), np.nanmean(np.cos(radians)))) % 360
            global_dev = circular_distance_deg(values, center)
            denom = float(np.nanmedian(global_dev))
            if not math.isfinite(denom) or denom <= 0:
                denom = float(np.nanstd(global_dev, ddof=1))
            for taxon, g in obs.groupby("taxon_name"):
                x = pd.to_numeric(g[trait], errors="coerce").dropna()
                if len(x) < min_n or not math.isfinite(denom) or denom <= 0:
                    continue
                rad = np.radians(x)
                ctr = np.degrees(np.arctan2(np.mean(np.sin(rad)), np.mean(np.cos(rad)))) % 360
                raw = float(np.median(circular_distance_deg(x, ctr)))
                rows.append({"taxon_name": taxon, "trait": trait, "module": MODULE[trait], "n": len(x), "within_variation_raw": raw, "within_variation_standardized": raw / denom})
        else:
            z = robust_scale(obs[trait])
            tmp = pd.DataFrame({"taxon_name": obs["taxon_name"], "z": z})
            for taxon, g in tmp.groupby("taxon_name"):
                x = g["z"].dropna()
                if len(x) < min_n:
                    continue
                raw = float(np.median(np.abs(x - np.median(x))))
                rows.append({"taxon_name": taxon, "trait": trait, "module": MODULE[trait], "n": len(x), "within_variation_raw": raw, "within_variation_standardized": raw})
    detail = pd.DataFrame(rows)
    if detail.empty:
        return detail
    species = detail.groupby("taxon_name", as_index=False).agg(
        species_within_variation_index=("within_variation_standardized", "mean"),
        n_traits_variation=("trait", "nunique"),
    )
    return detail.merge(species, on="taxon_name", how="left")


def combine_hue_coefficients(coef: pd.DataFrame) -> pd.DataFrame:
    coef = coef.copy()
    ok = coef[coef.get("status", "ok").eq("ok")].copy() if "status" in coef else coef.copy()
    hue = ok[ok["trait"].isin([HUE_SIN, HUE_COS])]
    non = ok[~ok["trait"].isin([HUE_SIN, HUE_COS])].copy()
    non["trait_endpoint"] = non["trait"]
    non["effect_magnitude"] = non["beta_std_within"].abs()
    non["effect_direction_degrees"] = np.where(non["beta_std_within"] >= 0, 0.0, 180.0)
    hue_rows = []
    for pred, g in hue.groupby("predictor"):
        vals = g.set_index("trait")["beta_std_within"]
        if HUE_SIN not in vals or HUE_COS not in vals:
            continue
        bs, bc = float(vals[HUE_SIN]), float(vals[HUE_COS])
        hue_rows.append({
            "trait": "corolla_hue_angle",
            "trait_endpoint": "corolla_hue_angle",
            "predictor": pred,
            "status": "ok",
            "beta_std_within": np.nan,
            "effect_magnitude": math.hypot(bs, bc),
            "effect_direction_degrees": math.degrees(math.atan2(bs, bc)) % 360.0,
            "hue_beta_sin": bs,
            "hue_beta_cos": bc,
            "n_observations": int(g["n_observations"].min()) if "n_observations" in g else np.nan,
            "n_species": int(g["n_species"].min()) if "n_species" in g else np.nan,
        })
    combined = pd.concat([non, pd.DataFrame(hue_rows)], ignore_index=True, sort=False)
    combined["module"] = combined["trait_endpoint"].map(MODULE)
    return combined


def responsiveness(combined: pd.DataFrame) -> pd.DataFrame:
    trait = combined.groupby(["trait_endpoint", "module"], as_index=False).agg(
        environmental_responsiveness_index=("effect_magnitude", lambda x: float(np.sqrt(np.mean(np.square(x.dropna()))))),
        n_predictors=("predictor", "nunique"),
    )
    module = combined.groupby("module", as_index=False).agg(
        module_environmental_responsiveness=("effect_magnitude", lambda x: float(np.sqrt(np.mean(np.square(x.dropna()))))),
        mean_abs_standardized_beta=("effect_magnitude", "mean"),
        n_effects=("effect_magnitude", "count"),
    )
    return trait, module


def plot_heatmap(combined: pd.DataFrame, out: Path) -> None:
    pivot = combined.pivot_table(index="trait_endpoint", columns="predictor", values="effect_magnitude", aggfunc="first")
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    im = ax.imshow(pivot.to_numpy(float), aspect="auto")
    ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=35, ha="right")
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    ax.set_title("Figure S1. Magnitude of standardized within-species effects")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("|standardized beta|; hue = joint vector magnitude")
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            v = pivot.iloc[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center")
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def plot_quadrants(variation: pd.DataFrame, trait_resp: pd.DataFrame, out: Path) -> None:
    vt = variation.groupby(["trait", "module"], as_index=False)["within_variation_standardized"].mean().rename(columns={"trait": "trait_endpoint", "within_variation_standardized": "within_variation_index"})
    p = vt.merge(trait_resp, on=["trait_endpoint", "module"], how="inner")
    xmid = p["within_variation_index"].median()
    ymid = p["environmental_responsiveness_index"].median()
    fig, ax = plt.subplots(figsize=(8, 7))
    for module, g in p.groupby("module"):
        ax.scatter(g["within_variation_index"], g["environmental_responsiveness_index"], label=module, s=70)
    for _, r in p.iterrows():
        ax.annotate(r["trait_endpoint"].replace("_median", ""), (r["within_variation_index"], r["environmental_responsiveness_index"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.axvline(xmid, linestyle="--", linewidth=1)
    ax.axhline(ymid, linestyle="--", linewidth=1)
    ax.set_xlabel("Within-species variation index")
    ax.set_ylabel("Environmental responsiveness index")
    ax.set_title("Trait lability decomposed into two independent axes")
    ax.legend(title="module")
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(args.observations, low_memory=False)
    coef = pd.read_csv(args.coefficients, low_memory=False)
    required_obs = {"taxon_name", *LINEAR_TRAITS, HUE_SIN, HUE_COS}
    required_coef = {"trait", "predictor", "beta_std_within"}
    if missing := sorted(required_obs.difference(obs.columns)):
        raise ValueError(f"Missing observation columns: {missing}")
    if missing := sorted(required_coef.difference(coef.columns)):
        raise ValueError(f"Missing coefficient columns: {missing}")

    variation = species_variation(obs, args.min_species_n)
    variation.to_csv(out / "species_trait_within_variation.csv", index=False, encoding="utf-8-sig")
    species_axes = variation[["taxon_name", "species_within_variation_index", "n_traits_variation"]].drop_duplicates()
    species_axes.to_csv(out / "species_lability_axes.csv", index=False, encoding="utf-8-sig")

    combined = combine_hue_coefficients(coef)
    combined.to_csv(out / "table_s1_all_coefficients.csv", index=False, encoding="utf-8-sig")
    trait_resp, module_resp = responsiveness(combined)
    trait_resp.to_csv(out / "trait_environmental_responsiveness.csv", index=False, encoding="utf-8-sig")
    module_resp.to_csv(out / "module_environmental_response_summary.csv", index=False, encoding="utf-8-sig")

    plot_heatmap(combined, out / "figure_s1_effect_size_heatmap")
    plot_quadrants(variation, trait_resp, out / "figure_lability_quadrants")


if __name__ == "__main__":
    main()
