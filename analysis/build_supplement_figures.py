#!/usr/bin/env python3
"""Build the five submission Supplement figures from committed frozen summary tables.

This script is presentation-only: it never refits a model or changes a frozen
cohort/result. Both PNG and SVG are written for each figure.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tables-dir", type=Path, default=Path("manuscript/supplement/tables"))
    p.add_argument("--out-dir", type=Path, default=Path("manuscript/supplement/figures"))
    return p.parse_args()


def save(fig: plt.Figure, out: Path, stem: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{stem}.png", dpi=240, bbox_inches="tight")
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def figure_s1(tables: Path, out: Path) -> None:
    df = pd.read_csv(tables / "S03_nested_visible_variance.csv")
    labels = df["label"].tolist()
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    left = np.zeros(len(df))
    for column, label in (
        ("among_taxon_fraction", "Among assigned taxon means"),
        ("among_photos_within_taxon_fraction", "Among photographs within taxa"),
        ("among_heads_within_photo_fraction", "Among heads within photograph"),
    ):
        values = df[column].to_numpy(float)
        ax.barh(y, values, left=left, label=label)
        left += values
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of total visible image sums of squares")
    ax.set_title("Figure S1. Nested visible-variance decomposition")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    save(fig, out, "Figure_S1_nested_visible_variance")


def figure_s2(tables: Path, out: Path) -> None:
    df = pd.read_csv(tables / "S04_full_primary_within_taxon_coefficients.csv")
    linear = df.loc[df["endpoint_type"].eq("linear")].copy()
    traits = list(dict.fromkeys(linear["trait_endpoint"]))
    predictors = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15"]
    matrix = np.full((len(traits), len(predictors)), np.nan)
    qmat = np.full_like(matrix, np.nan)
    for i, trait in enumerate(traits):
        for j, predictor in enumerate(predictors):
            row = linear.loc[linear["trait_endpoint"].eq(trait) & linear["predictor"].eq(predictor)]
            if not row.empty:
                matrix[i, j] = float(row.iloc[0]["beta_std_within"])
                qmat[i, j] = float(row.iloc[0]["q_fdr_bh"])
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    image = ax.imshow(matrix, aspect="auto")
    fig.colorbar(image, ax=ax, label="Standardized within-taxon coefficient")
    ax.set_xticks(range(4), ["BIO1", "BIO4", "BIO12", "BIO15"])
    short = [x.replace("_median", "").replace("corolla_lab_", "").replace("shape_", "") for x in traits]
    ax.set_yticks(range(len(short)), short)
    for i in range(len(traits)):
        for j in range(4):
            if np.isfinite(qmat[i, j]) and qmat[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=13)
    ax.set_title("Figure S2. Exhaustive primary climate coefficients (* BH q < 0.05)")
    save(fig, out, "Figure_S2_primary_coefficient_map")


def figure_s3(tables: Path, out: Path) -> None:
    df = pd.read_csv(tables / "S05_spde_model_group_summary.csv")
    traits = list(dict.fromkeys(df["trait"]))
    groups = ["climate", "climate_topography", "climate_soil", "full"]
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    x = np.arange(len(traits))
    width = 0.18
    for offset, group in enumerate(groups):
        sub = df.loc[df["model_group"].eq(group)].set_index("trait")
        values = [float(sub.loc[t, "delta_waic_within_trait"]) for t in traits]
        ax.bar(x + (offset - 1.5) * width, values, width=width, label=group)
    ax.set_xticks(x, [t.replace("_median", "").replace("corolla_lab_", "").replace("shape_", "") for t in traits], rotation=45, ha="right")
    ax.set_ylabel("ΔWAIC within endpoint")
    ax.set_title("Figure S3. Grouped SPDE-INLA model comparison")
    ax.legend(frameon=False, fontsize=8)
    save(fig, out, "Figure_S3_spde_delta_waic")


def figure_s4_s5(tables: Path, out: Path) -> None:
    df = pd.read_csv(tables / "S09_spatial_robustness.csv")
    y = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(7.6, 5.1))
    ax.barh(y, df["residual_morans_I"].astype(float))
    ax.axvline(0, linewidth=0.8)
    ax.set_yticks(y, df["endpoint"])
    ax.invert_yaxis()
    ax.set_xlabel("Residual Moran's I")
    ax.set_title("Figure S4. Residual spatial autocorrelation")
    save(fig, out, "Figure_S4_residual_moran")

    ordered = df.sort_values("minimum_leave_one_region_out_spearman_rho")
    y = np.arange(len(ordered))
    fig, ax = plt.subplots(figsize=(7.6, 5.1))
    ax.barh(y, ordered["minimum_leave_one_region_out_spearman_rho"].astype(float))
    ax.set_yticks(y, ordered["endpoint"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Minimum leave-one-region-out Spearman ρ")
    ax.set_title("Figure S5. Broad-region omission stability")
    save(fig, out, "Figure_S5_leave_one_region_out")


def main() -> None:
    a = args()
    figure_s1(a.tables_dir, a.out_dir)
    figure_s2(a.tables_dir, a.out_dir)
    figure_s3(a.tables_dir, a.out_dir)
    figure_s4_s5(a.tables_dir, a.out_dir)
    print(f"Wrote Supplement figures to {a.out_dir}")


if __name__ == "__main__":
    main()
