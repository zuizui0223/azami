#!/usr/bin/env python3
"""Build supervisor-review figures from frozen Chapter 1 output tables.

This is a presentation layer only. It does not refit models or modify any frozen
scientific result. Each figure is written as both PNG and SVG.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd

FRIENDLY = {
    "orientation_angle_degrees_median": "Orientation",
    "corolla_lab_lightness_median": "Lightness",
    "corolla_lab_chroma_median": "Chroma",
    "corolla_hue_sin_median": "Hue sine",
    "corolla_hue_cos_median": "Hue cosine",
    "shape_aspect_ratio_median": "Aspect ratio",
    "shape_circularity_median": "Circularity",
    "shape_solidity_median": "Solidity",
    "shape_width_cv_median": "Width-profile CV",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-screen", required=True)
    parser.add_argument("--pca-variance", required=True)
    parser.add_argument("--pca-loadings", required=True)
    parser.add_argument("--species-lability", required=True)
    parser.add_argument("--module-summary", required=True)
    parser.add_argument("--niche-contrasts", required=True)
    parser.add_argument("--sample-sensitivity", required=True)
    parser.add_argument("--molecular-coverage", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def save(fig: plt.Figure, out: Path, stem: str) -> None:
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    trait_screen = pd.read_csv(args.trait_screen)
    pca_var = pd.read_csv(args.pca_variance)
    pca_load = pd.read_csv(args.pca_loadings)
    axes = pd.read_csv(args.species_lability)
    modules = pd.read_csv(args.module_summary)
    niche = pd.read_csv(args.niche_contrasts)
    sensitivity = pd.read_csv(args.sample_sensitivity)
    coverage = pd.read_csv(args.molecular_coverage)

    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.axis("off")
    boxes = [
        (0.13, 0.73, "Public observations\n3,725 observations\n216 image-analysis taxa"),
        (0.50, 0.73, "Capitulum detection\n6,626 detected heads\ntight + context crops"),
        (0.87, 0.73, "Continuous traits\norientation, colour, outline\ntrait-specific assessability"),
        (0.18, 0.28, "Species distributions\nwithin/among variance\nPCA and niche contrasts"),
        (0.53, 0.28, "Strict spatial layer\n46,276 thinned observations\n259 input taxa"),
        (0.84, 0.28, "Two-axis lability\n102 complete taxa\nvariation != responsiveness"),
    ]
    for x, y, text in boxes:
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=12,
            bbox=dict(boxstyle="round,pad=0.6", alpha=0.18),
        )
    for start, end in [((0.25, 0.73), (0.38, 0.73)), ((0.62, 0.73), (0.75, 0.73))]:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="->",
                mutation_scale=18,
                linewidth=1.5,
            )
        )
    for start, end in [((0.80, 0.64), (0.24, 0.39)), ((0.86, 0.64), (0.55, 0.39))]:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="->",
                mutation_scale=18,
                linewidth=1.5,
                connectionstyle="arc3,rad=0.08",
            )
        )
    ax.add_patch(
        FancyArrowPatch(
            (0.65, 0.28),
            (0.73, 0.28),
            arrowstyle="->",
            mutation_scale=18,
            linewidth=1.5,
        )
    )
    ax.set_title(
        "Figure 1. Scale-explicit public-image phenomics workflow",
        fontsize=17,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.07,
        "Distinct cohorts answer distinct questions; measurement failure is retained as missingness rather than biological absence.",
        ha="center",
        fontsize=11,
    )
    save(fig, out, "figure1_pipeline")

    partition = trait_screen.copy()
    partition["label"] = partition["trait"].map(FRIENDLY)
    partition = partition.sort_values("within_species_variance_fraction")
    fig, ax = plt.subplots(figsize=(9.2, 6.5))
    y_positions = np.arange(len(partition))
    ax.barh(
        y_positions,
        partition["within_species_variance_fraction"],
        label="Within species",
    )
    ax.barh(
        y_positions,
        partition["between_species_variance_fraction"],
        left=partition["within_species_variance_fraction"],
        label="Among species",
    )
    ax.set_yticks(y_positions, partition["label"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of visible variance")
    ax.set_title("Figure 2. Most visible trait variance occurs within species")
    for index, value in enumerate(partition["within_species_variance_fraction"]):
        ax.text(
            min(value - 0.01, 0.97),
            index,
            f"{value:.3f}",
            ha="right",
            va="center",
            fontsize=9,
        )
    ax.legend(loc="lower right")
    fig.tight_layout()
    save(fig, out, "figure2_variance_partition")

    loadings = pca_load.copy()
    loadings["label"] = loadings["trait"].map(FRIENDLY)
    matrix = loadings[["PC1", "PC2", "PC3"]].to_numpy()
    fig, ax = plt.subplots(figsize=(7.8, 6.6))
    image = ax.imshow(matrix, aspect="auto")
    component_labels = [
        f"{row.component}\n{row.variance_explained * 100:.1f}%"
        for _, row in pca_var.iloc[:3].iterrows()
    ]
    ax.set_xticks(range(3), component_labels)
    ax.set_yticks(range(len(loadings)), loadings["label"])
    ax.set_title("Figure 3. Species-level trait architecture is multidimensional")
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            ax.text(
                column_index,
                row_index,
                f"{matrix[row_index, column_index]:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )
    fig.colorbar(image, ax=ax, label="PCA loading")
    fig.tight_layout()
    save(fig, out, "figure3_pca_loadings")

    complete = axes[axes["primary_complete"].astype(bool)].copy()
    x_values = complete["species_within_variation_index"]
    y_values = complete["species_environmental_responsiveness_index"]
    x_median, y_median = x_values.median(), y_values.median()
    fig, ax = plt.subplots(figsize=(9.3, 6.8))
    ax.scatter(x_values, y_values, s=35, alpha=0.72)
    ax.axvline(x_median, linestyle="--", linewidth=1.2)
    ax.axhline(y_median, linestyle="--", linewidth=1.2)
    ax.set_xlabel("Species within-variation index")
    ax.set_ylabel("Species environmental-responsiveness index")
    ax.set_title(
        "Figure 4. Variation and environmental responsiveness are distinct axes"
    )
    rho = x_values.corr(y_values, method="spearman")
    ax.text(
        0.02,
        0.50,
        f"Spearman rho = {rho:.3f}\nN = {len(complete)} taxa",
        transform=ax.transAxes,
        fontsize=10,
    )
    fig.tight_layout()
    save(fig, out, "figure4_lability_quadrants")

    module = modules.copy()
    module["module"] = pd.Categorical(
        module["module"],
        ["colour", "orientation", "shape"],
        ordered=True,
    )
    module = module.sort_values("module")
    x_positions = np.arange(len(module))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    ax.bar(
        x_positions - width / 2,
        module["median_species_within_variation"],
        width=width,
        label="Within variation",
    )
    ax.bar(
        x_positions + width / 2,
        module["median_species_environmental_responsiveness"],
        width=width,
        label="Environmental responsiveness",
    )
    ax.set_xticks(
        x_positions,
        [value.title() for value in module["module"].astype(str)],
    )
    ax.set_ylabel("Median species-level index")
    ax.set_title(
        "Figure 5. Trait modules occupy different regions of lability space"
    )
    ax.legend()
    fig.tight_layout()
    save(fig, out, "figure5_module_lability")

    contrast = niche.copy()
    contrast["label"] = contrast["trait"].map(FRIENDLY)
    fig, ax = plt.subplots(figsize=(8.8, 6.4))
    ax.scatter(
        contrast["gaussian_bhattacharyya_overlap"],
        contrast["environmental_centroid_distance"],
        s=60,
    )
    for _, row in contrast.iterrows():
        ax.annotate(
            row["label"],
            (
                row["gaussian_bhattacharyya_overlap"],
                row["environmental_centroid_distance"],
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel("Gaussian Bhattacharyya niche overlap")
    ax.set_ylabel("Environmental centroid distance")
    ax.set_title(
        "Figure 6. Environmental sorting is strongest for orientation and colour"
    )
    fig.tight_layout()
    save(fig, out, "figure6_niche_contrasts")

    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    ax.plot(
        sensitivity["min_species_n"],
        sensitivity["variation_spearman"],
        marker="o",
        label="Within-variation ranking",
    )
    ax.plot(
        sensitivity["min_species_n"],
        sensitivity["responsiveness_spearman"],
        marker="o",
        label="Responsiveness ranking",
    )
    ax.set_xticks(sensitivity["min_species_n"])
    ax.set_ylim(0.85, 1.01)
    ax.set_xlabel("Minimum observations per species-trait combination")
    ax.set_ylabel("Spearman correlation with primary ranking")
    ax.set_title(
        "Figure 7. The two-axis result is stable to the minimum sample threshold"
    )
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, out, "figure7_sample_sensitivity")

    n_taxa = 216
    values = {
        "Direct dated-backbone tips": 54,
        "Any nucleotide record": int((coverage["nuccore_all"] > 0).sum()),
        "ITS record": int((coverage["nuccore_its"] > 0).sum()),
        "Plastid record": int((coverage["nuccore_plastid"] > 0).sum()),
        "Complete plastome": int(
            (coverage["nuccore_complete_plastome"] > 0).sum()
        ),
        "SRA record": int((coverage["sra_all"] > 0).sum()),
    }
    labels = list(values)
    counts = list(values.values())
    percentages = [count / n_taxa * 100 for count in counts]
    fig, ax = plt.subplots(figsize=(9.0, 5.7))
    y_positions = np.arange(len(labels))
    ax.barh(y_positions, percentages)
    ax.set_yticks(y_positions, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Coverage among 216 image-analysis taxa (%)")
    ax.set_title(
        "Figure 8. Historical inference is limited by uneven direct and molecular coverage"
    )
    for index, (count, percentage) in enumerate(zip(counts, percentages)):
        ax.text(
            percentage + 1,
            index,
            f"{count}/{n_taxa} ({percentage:.1f}%)",
            va="center",
            fontsize=9,
        )
    fig.tight_layout()
    save(fig, out, "figure8_historical_molecular_coverage")


if __name__ == "__main__":
    main()
