#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from PIL import Image

plt.rcParams["svg.hashsalt"] = "azami-figure1-multiscale"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input-figure", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def save(fig, out: Path, stem: str):
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight", metadata={"Date": None})
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight", metadata={"Software": "azami"})
    fig.savefig(out / f"{stem}.pdf", bbox_inches="tight", metadata={"Creator": "azami", "CreationDate": None, "ModDate": None})
    plt.close(fig)


def main():
    a = parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    image = Image.open(a.input_figure).convert("RGB")

    fig = plt.figure(figsize=(13.2, 13.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.6, 1.25], hspace=0.05)

    ax = fig.add_subplot(gs[0])
    ax.imshow(image)
    ax.axis("off")

    ax = fig.add_subplot(gs[1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.01, 0.95, "J   From one measured capitulum to multiscale ecological questions", fontsize=12.5, fontweight="bold", va="top")
    boxes = [
        (0.04, 0.48, "Repeated numerical\nhead phenotypes", "orientation · colour · outline"),
        (0.30, 0.48, "Nested visible variance", "taxon → photograph → head"),
        (0.56, 0.48, "Within-taxon environment", "demeaned climate models + SPDE"),
        (0.82, 0.48, "Among-taxon structure", "niche sorting + PGLS sensitivity"),
    ]
    for x, y, title, subtitle in boxes:
        ax.text(x, y, f"{title}\n{subtitle}", ha="center", va="center", fontsize=9.2,
                bbox=dict(boxstyle="round,pad=0.55", facecolor="white", edgecolor="#666666", linewidth=0.9))
    for x1, x2 in [(0.14, 0.20), (0.40, 0.46), (0.66, 0.72)]:
        ax.add_patch(FancyArrowPatch((x1, 0.48), (x2, 0.48), arrowstyle="->", mutation_scale=14, linewidth=1.0))
    ax.text(0.5, 0.10,
            "Distinct executed cohorts answer distinct questions. Image assessability/failure remains measurement missingness, not biological absence.",
            ha="center", va="center", fontsize=8.5)

    fig.suptitle("Figure 1. From public photographs to continuous multiscale thistle phenomics", fontsize=16, fontweight="bold", y=0.995)
    save(fig, out, "Figure_1_image_to_multiscale_analysis")


if __name__ == "__main__":
    main()
