#!/usr/bin/env python3
"""Exploratory Figure 2 from the existing SHA-verified 28-pair table; no new fit."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FROZEN = {
    ("angle", "chroma"): (-0.22512839486857492, -0.028067356918921597),
    ("angle", "aspect"): (-0.24694620634514822, -0.0693870326972044),
    ("chroma", "aspect"): (0.01993972231831634, -0.01824691397495542),
}
LABELS = ["Image angle x chroma", "Image angle x outline aspect",
          "Chroma x outline aspect"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    with args.pairs_csv.open(newline="", encoding="utf-8") as f:
        records = list(csv.DictReader(f))
    if len(records) != 28:
        raise ValueError("Expected original 28 pairs")
    rows = {(r["left"], r["right"]): r for r in records}
    if len(rows) != 28:
        raise ValueError("Duplicate pairs")
    fig, ax = plt.subplots(figsize=(8, 4))
    for i, (pair, expected) in enumerate(FROZEN.items()):
        item = rows[pair]
        if (item["n_photos"], item["n_taxa"]) != ("156683", "174"):
            raise ValueError("Cohort drift")
        estimates = (float(item["adjusted_among_signed_r"]),
                     float(item["adjusted_within_signed_r"]))
        if any(abs(a-b)>1e-9 for a,b in zip(estimates, expected)):
            raise ValueError(f"Original correlation drift: {pair}")
        y = 2-i
        ax.scatter(estimates[0], y+.12, marker="o", s=55, color="tab:blue",
                   label="Among taxa" if i==0 else None)
        ax.scatter(estimates[1], y-.12, marker="s", s=50, color="tab:orange",
                   label="Within taxa" if i==0 else None)
        ax.text(estimates[0]-.015, y+.12, f"{estimates[0]:+.3f}",
                va="center", ha="right", fontsize=9)
        ax.text(estimates[1]+.015, y-.12, f"{estimates[1]:+.3f}",
                va="center", ha="left", fontsize=9)
    ax.axvline(0, linestyle=":", color=".55")
    ax.set_xlim(-.38, .17)
    ax.set_ylim(-.5, 2.65)
    ax.set_yticks((2,1,0), LABELS)
    ax.set_xlabel("Environment-adjusted descriptive photographic correlation (r)")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("Cross-domain trait associations (exploratory)", loc="left")
    fig.text(.01, .012, "Same 156,683 photographs / 174 taxa. No functional or phylogenetic test.\n"
             "Image angle: 0 deg up, 180 deg down relative to the image frame.",
             fontsize=8, color=".4")
    fig.tight_layout(rect=(0,.10,1,1))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(args.output_dir / f"figure2_cross_functional_trait_pairs.{extension}",
                    dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
