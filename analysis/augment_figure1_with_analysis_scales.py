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

    # Preserve the validated actual-photo Figure 1 at its natural aspect ratio and
    # append only the analysis-scale strip. The original internal title is retained,
    # so no second super-title is added.
    fig = plt.figure(figsize=(13.2, 11.35))
    gs = fig.add_gridspec(2, 1, height_ratios=[9.2, 2.0], hspace=0.01)

    ax = fig.add_subplot(gs[0])
    ax.imshow(image)
    ax.axis("off")

    ax = fig.add_subplot(gs[1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.01, 0.92, "j   From measured capitula to multiscale ecological questions", fontsize=12.0, fontweight="bold", va="top")
    boxes = [
        (0.10, 0.49, "Repeated numerical\nhead phenotypes", "orientation · colour · outline"),
        (0.36, 0.49, "Nested visible variance", "taxon → photograph → head"),
        (0.62, 0.49, "Within-taxon environment", "demeaned climate models + SPDE"),
        (0.88, 0.49, "Among-taxon structure", "niche sorting + PGLS sensitivity"),
    ]
    for x, y, title, subtitle in boxes:
        ax.text(x, y, f"{title}\n{subtitle}", ha="center", va="center", fontsize=8.9,
                bbox=dict(boxstyle="round,pad=0.48", facecolor="white", edgecolor="#666666", linewidth=0.85))
    for x1, x2 in [(0.19, 0.27), (0.45, 0.53), (0.71, 0.79)]:
        ax.add_patch(FancyArrowPatch((x1, 0.49), (x2, 0.49), arrowstyle="->", mutation_scale=13, linewidth=0.95))
    ax.text(0.5, 0.10,
            "Distinct executed cohorts answer distinct questions; image assessability/failure remains measurement missingness, not biological absence.",
            ha="center", va="center", fontsize=8.2)

    save(fig, out, "Figure_1_image_to_multiscale_analysis")


if __name__ == "__main__":
    main()
