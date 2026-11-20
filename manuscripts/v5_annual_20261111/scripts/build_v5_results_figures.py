#!/usr/bin/env python3
"""Regenerate V5 GEB result panels from the adopted *annual five-axis* source CSVs.

Outputs: scalable SVG, vector PDF, high-resolution PNG. This script only
reads committed aggregate results, never individual observation coordinates.
No confidence intervals are invented for the native analogue distance tests.

CLI:
  python manuscripts/v5_annual_20261111/scripts/build_v5_results_figures.py \
    --repo-root . --output manuscripts/v5_annual_20261111/figures/generated
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK = "#20324b"
TEAL = "#0a8e82"
BLUE = "#2475b5"
PURPLE = "#6e40a5"
GREY = "#758295"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

SOURCES = {
    "within": "analysis/v5/results/native_seasonality/v5_within_analogue_seasonal.csv",
    "between": "analysis/v5/results/native_seasonality/v5_cross_analogue_seasonal.csv",
    "slopes": "analysis/v5/results/native_seasonality/v5_signed_pairs_seasonal.csv",
    "head_d": "analysis/v5/results/novel_head_overlap/novel_head_overlap.csv",
    "chroma": "analysis/v5/results/invasion_design/climate_balanced_heads.csv",
}

def get_rows(repo: Path, key: str) -> list[dict[str, str]]:
    with (repo / SOURCES[key]).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))

def ensure(ok: bool, why: str):
    if not ok:
        raise ValueError(why)

def save(fig, out: Path, stem: str):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300, bbox_inches="tight",
                    pad_inches=0.15)
    plt.close(fig)

def build_analogue(repo: Path, out: Path):
    a = get_rows(repo, "within")
    b = get_rows(repo, "between")
    ensure(len(a) == 2 and len(b) == 2, "FIVE_AXIS_GRAIN_COUNTS_CHANGED")
    rows = ([("within", x) for x in a] +
            [("between", x) for x in b])
    val = [float(x["excess"]) for _, x in rows]
    p = [float(x["holm"]) for _, x in rows]
    ensure(np.allclose(val, [.08421227833291306, .23836924606855578,
                             .030937948818972237, .06532956260499098]),
           "ORIGINAL_FIVE_AXIS_EXCESS_CHANGED")
    ensure(np.allclose(p, [.004, .004, .076, .012]),
           "ORIGINAL_HOLM_FAMILY_CHANGED")
    fig, ax = plt.subplots(figsize=(9.5, 3.8))
    y = [3, 2, 1, 0]
    labels = ["Within taxon · 1°", "Within taxon · 5°",
              "Between taxa · 1°", "Between taxa · 5°"]
    colors = [TEAL, TEAL, "#B3BAC5", PURPLE]
    for yy, v, color, pv in zip(y, val, colors, p):
        ax.plot([0, v], [yy, yy], color=color, lw=6, alpha=.55)
        ax.scatter([v], [yy], s=95, color=color, zorder=3)
        ax.text(.267, yy, f"+{v:.3f}  |  Holm P={pv:.3f}",
                va="center", fontsize=10, color=INK)
    ax.set_xlim(-.01, .265)
    ax.set_yticks(y, labels)
    ax.set_ylim(-.65, 3.65)
    ax.grid(axis="x", color="#dce6ed")
    ax.set_axisbelow(True)
    ax.set_xlabel("Excess joint capitulum resemblance (standardised-distance reduction)")
    ax.set_title("Annual five-axis climatic analogues and observed capitulum resemblance",
                 loc="left", fontweight="bold", fontsize=13)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=.26, right=.71, top=.85, bottom=.20)
    fig.text(.035, .02,
             "24 rarefactions; 499 within-taxon whole-vector permutations. "
             "1°: 63 taxa / 7,117 taxon-cells; 5°: 31 taxa / 814 taxon-cells. "
             "No source CI for this summary.", fontsize=8, color=GREY)
    save(fig, out, "figure2_annual_climate_analogue_resemblance")

def build_slopes(repo: Path, out: Path):
    r = get_rows(repo, "slopes")
    matched = {
        (x["axis"], x["trait"], int(x["degrees"])): x
        for x in r
        if x["sample"] == "all_cells" and x["scale"] == "within_pair"
        and (x["axis"], x["trait"]) in {
            ("annual_rain", "angle"),
            ("precipitation_seasonality", "chroma"),
        }
    }
    ensure(len(matched) == 4, "SOURCE_ALL_CELL_WITHIN_PAIR_SLOPES_MISSING")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 3.5))
    specs = [
        ("annual_rain", "angle", "A  BIO12: rainfall → image-referenced angle",
         BLUE),
        ("precipitation_seasonality", "chroma", "B  BIO15: seasonality → photo C*",
         PURPLE),
    ]
    for ax, (axis, trait, title, color) in zip(axes, specs):
        for y, deg in [(1, 1), (0, 5)]:
            row = matched[(axis, trait, deg)]
            beta = float(row["beta"])
            lo, hi = float(row["lo"]), float(row["hi"])
            ax.errorbar(beta, y, xerr=[[beta-lo], [hi-beta]],
                        fmt="o", color=color, lw=2.0, capsize=4, markersize=7)
            ax.text(.405, y, f"{beta:+.3f}", va="center",
                    color=color, fontweight="bold")
        ax.axvline(0, color=GREY, lw=1, ls="--")
        ax.set_yticks([1, 0], ["1°", "5°"])
        ax.set_xlim(-.07, .405)
        ax.set_ylim(-.65, 1.6)
        ax.grid(axis="x", color="#e1e8ee")
        ax.set_axisbelow(True)
        ax.set_title(title, loc="left", fontsize=11.0, fontweight="bold")
        ax.set_xlabel("Standardised slope β")
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=.08, right=.89, bottom=.26, top=.82, wspace=.38)
    fig.text(.04, .06,
             "Pointwise taxon-bootstrap 95% intervals only; exploratory. "
             "Both grains reuse the same photographic source; "
             "neither physical gravity nor pigment concentration was measured.",
             fontsize=8, color=GREY)
    save(fig, out, "figure3_annual_trait_climate_gradients")

def build_native_introduced(repo: Path, out: Path):
    rows = get_rows(repo, "head_d")
    head = {
        x["species"]: x for x in rows
        if float(x["native95_mass"]) == .95 and int(x["n_traits"]) == 3
        and float(x["bandwidth"]) == .5 and x["design"] == "geographic_pairs"
    }
    chroma = {
        x["species"]: x for x in get_rows(repo, "chroma")
        if x["trait"] == "chroma" and int(x["bin_width"]) == 4
    }
    ensure(set(head) == {"Cirsium arvense", "Cirsium vulgare"},
           "ORIGINAL_GEOGRAPHIC_PAIR_COHORT_CHANGED")
    ensure(set(chroma) == set(head), "BALANCED_CHROMA_SOURCE_MISSING")
    fig, ax = plt.subplots(figsize=(9.2, 3.7))
    colors = [BLUE, PURPLE]
    for i, species in enumerate(["Cirsium arvense", "Cirsium vulgare"]):
        x = head[species]
        in_core = float(x["D_native_analogous"])
        outside = float(x["D_native_novel"])
        p = float(x["holm_primary"])
        y = 1 - i
        ax.plot([in_core, outside], [y, y], lw=5, color=colors[i], alpha=.35)
        ax.scatter([in_core], [y], s=95, color=colors[i], zorder=3)
        ax.scatter([outside], [y], s=95, facecolors="white",
                   edgecolors=colors[i], linewidths=2, zorder=4)
        ax.text(.91, y+.12, f"D={in_core:.3f} → {outside:.3f}",
                color=colors[i], fontweight="bold")
        ax.text(.91, y-.19, f"Holm P={p:.3f}", color=GREY)
    ax.set_yticks([1, 0], ["C. arvense", "C. vulgare"])
    ax.set_xlim(.63, .91); ax.set_ylim(-.6, 1.6)
    ax.set_xlabel("Native-reference overlap of photographed three-trait space (D)")
    ax.set_title("Native-core-analogous vs outside-native-core introduced climates",
                 loc="left", fontsize=12, fontweight="bold")
    ax.grid(axis="x", color="#e4e8ed"); ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=.18, right=.79, bottom=.32, top=.80)
    row = chroma["Cirsium vulgare"]
    fig.text(.03, .09,
             f"Introduced minus native C. vulgare chroma: "
             f"{float(row['shift_raw']):+.3f} C*; shared-climate-bin "
             f"frequency balancing: {float(row['shift_balanced']):+.3f} C* "
             "(descriptive, not causal).", fontsize=9.5, color=PURPLE)
    fig.text(.03, .035,
             "Regions: arvense 30 pairs/9 blocks; vulgare 146 pairs/14 blocks. "
             "No familywise-supported outside-core head difference.", fontsize=8, color=GREY)
    save(fig, out, "figure5_native_introduced_head_and_chroma")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path("."))
    p.add_argument("--output", type=Path, required=True)
    opt = p.parse_args()
    repo = opt.repo_root.resolve()
    out = opt.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    build_analogue(repo, out)
    build_slopes(repo, out)
    build_native_introduced(repo, out)
    receipt = {
        "source_hashes": {
            k: hashlib.sha256((repo / v).read_bytes()).hexdigest()
            for k, v in SOURCES.items()
        },
        "statistical_scope": {
            "figure2": "Annual five-axis; original Holm; no invented CIs",
            "figure3": "All eligible cells, within taxon pair, pointwise CIs",
            "figure5": "Region matched D; C* contrasts are descriptive",
        },
        "source_map": (
            "Any observed-location map must use an audit-backed cell table. "
            "The existing v4 plot_sampling_coverage_v4.py reads verified "
            "V4 5-degree cells; these are NOT exact V5 three-trait-cohort points."
        ),
    }
    (out / "figure_source_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

if __name__ == "__main__":
    main()
