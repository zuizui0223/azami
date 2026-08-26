#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

V2_PATH = Path(__file__).with_name("build_ch1_interpretive_figures_v2.py")
spec = importlib.util.spec_from_file_location("ch1_interpretive_v2", V2_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot load {V2_PATH}")
v2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v2)
core = v2.core


def fixed_plot_figure4(scores, load, variance, out):
    fig, axes = core.plt.subplots(1, 2, figsize=(12.3, 5.9), gridspec_kw={"width_ratios": [1.25, 0.9]})
    ax = axes[0]
    ax.scatter(scores["PC1"], scores["PC2"], s=20, alpha=0.45, linewidths=0)
    ax.axhline(0, linewidth=0.6, linestyle="--")
    ax.axvline(0, linewidth=0.6, linestyle="--")
    ax.set_xlim(-4.2, 3.2)
    ax.set_ylim(-4.3, 3.45)

    scale = 5.3
    labels = {
        "orientation_angle_degrees_median": ("Orientation", -0.45, 1.55),
        "corolla_lab_lightness_median": ("Lightness", -2.00, 2.35),
        "corolla_lab_chroma_median": ("Chroma", 1.55, -2.60),
        "corolla_hue_sin_median": ("Hue sine", -0.75, 2.93),
        "corolla_hue_cos_median": ("Hue cosine", 2.35, 0.48),
        "shape_aspect_ratio_median": ("Aspect ratio", -0.35, -1.05),
        "shape_circularity_median": ("Circularity", -2.55, -2.35),
        "shape_solidity_median": ("Solidity", -2.60, -1.52),
        "shape_width_cv_median": ("Width-profile CV", 2.55, 1.02),
    }
    for _, row in load.iterrows():
        x, y = float(row["PC1"]) * scale, float(row["PC2"]) * scale
        ax.annotate("", xy=(x, y), xytext=(0, 0), arrowprops=dict(arrowstyle="->", linewidth=1.0))
        label, lx, ly = labels[row["trait"]]
        ax.text(lx, ly, label, fontsize=7.5, ha="center", va="center")

    ax.set_xlabel(f"PC1 ({variance[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance[1]*100:.1f}%)")
    ax.set_title("a  Taxon scores and trait directions\ncentral view; 148 taxa complete for all nine endpoints", fontsize=10)

    # Preserve the full score range without allowing one extreme taxon to compress the main panel.
    inset = ax.inset_axes([0.03, 0.05, 0.27, 0.24])
    inset.scatter(scores["PC1"], scores["PC2"], s=7, alpha=0.45, linewidths=0)
    inset.axhline(0, linewidth=0.35); inset.axvline(0, linewidth=0.35)
    inset.set_title("full PC1–PC2 range", fontsize=6.5)
    inset.tick_params(labelsize=5, length=1.5)
    extreme = scores.loc[scores["PC1"].idxmin()]
    inset.annotate(str(extreme["taxon_name"]).replace("Cirsium ", "C. "),
                   (extreme["PC1"], extreme["PC2"]), xytext=(5, 5), textcoords="offset points", fontsize=5.5)
    ax.text(0.02, 0.01, f"Inset retains the PC1 extreme: {extreme['taxon_name']} = {extreme['PC1']:.1f}",
            transform=ax.transAxes, fontsize=6.7, va="bottom")

    ax = axes[1]
    work = load.copy(); work["label"] = work["trait"].map(core.FRIENDLY); work = work.iloc[::-1]
    y = core.np.arange(len(work))
    ax.barh(y, work["PC3"])
    ax.axvline(0, linewidth=0.7)
    ax.set_yticks(y, work["label"])
    ax.set_xlabel(f"PC3 loading ({variance[2]*100:.1f}% variance)")
    ax.set_title("b  Third independent trait dimension", fontsize=10)
    ax.text(0.02, 0.02, f"PC1–PC3 cumulative variance = {variance[:3].sum()*100:.1f}%", transform=ax.transAxes, fontsize=8)
    fig.suptitle("Taxon-level capitulum trait architecture is multidimensional", fontsize=13)
    fig.tight_layout()
    core.save(fig, out, "Figure_4_taxon_trait_architecture")


def fixed_plot_s8(env, niche_summary, out):
    metrics = core.niche_metrics(env)
    support = {}; direct = {}
    for row in niche_summary.itertuples(index=False):
        ep = str(row[0])
        support[ep] = (float(row[1]) < 0.05) or (float(row[2]) < 0.05)
        direct[ep] = (float(row[3]) < 0.05) or (float(row[4]) < 0.05)
    reverse = {
        "orientation_angle_degrees_median": "Orientation", "corolla_lab_lightness_median": "Lightness",
        "corolla_lab_chroma_median": "Chroma", "corolla_hue_sin_median": "Hue sine",
        "corolla_hue_cos_median": "Hue cosine", "shape_aspect_ratio_median": "Aspect ratio",
        "shape_circularity_median": "Circularity", "shape_solidity_median": "Solidity",
        "shape_width_cv_median": "Width-profile CV",
    }
    labels = [reverse[t] for t in metrics["trait"]]
    supported = core.np.array([support.get(x, False) for x in labels])
    direct_supported = core.np.array([direct.get(x, False) for x in labels])
    sizes = core.np.where(supported, 72, 38)

    fig, ax = core.plt.subplots(figsize=(8.4, 6.2))
    ax.scatter(metrics["overlap"], metrics["centroid_distance"], s=sizes, alpha=0.75)
    ring = metrics.loc[direct_supported]
    if len(ring):
        ax.scatter(ring["overlap"], ring["centroid_distance"], s=145, facecolors="none", edgecolors="black", linewidths=1.0)
    for (_, row), label in zip(metrics.iterrows(), labels):
        ax.annotate(label, (row["overlap"], row["centroid_distance"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Gaussian Bhattacharyya overlap (lower = stronger separation)")
    ax.set_ylabel("Environmental centroid distance (higher = stronger separation)")
    ax.set_title("Environmental sorting of low- versus high-trait taxa")
    ax.text(0.01, -0.13,
            "Larger points: at least one all-taxa metric BH q < 0.05. Black ring: direct-backbone subset has at least one supported metric.",
            transform=ax.transAxes, fontsize=7.5)
    fig.tight_layout()
    core.save(fig, out, "Figure_S8_environmental_sorting")
    metrics.to_csv(out / "FigureS8_environmental_sorting_metrics.csv", index=False)


core.plot_figure4 = fixed_plot_figure4
core.plot_s8 = fixed_plot_s8

if __name__ == "__main__":
    core.main()
