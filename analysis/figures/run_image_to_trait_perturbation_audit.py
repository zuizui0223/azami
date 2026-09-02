#!/usr/bin/env python3
"""Render the frozen automated image-to-trait technical-audit summary.

The summary JSON and generated figures are external prepublication artifacts;
both paths must be supplied explicitly.  This script never reconstructs
unavailable per-image/per-perturbation values and does not treat human labels as
truth.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import azami_figstyle as fs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def clean_axis(ax: plt.Axes, *, grid_axis: str = "x") -> None:
    ax.grid(axis=grid_axis, linewidth=0.4, color="#D9D9D9")
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)


def save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fs.savefig(
        fig,
        stem,
        width="double",
        outdir=output_dir,
        metadata={"Creator": "Azami automated technical-audit figure builder"},
    )
    plt.close(fig)


def figure_s1(summary: dict, output_dir: Path) -> None:
    audit = summary["mirror_audit"]
    usable = audit["usable_fraction"]
    labels = ["Colour", "Gross outline", "Orientation"]
    values = np.array([usable["colour"], usable["outline"], usable["orientation"]]) * 100
    p95 = float(audit["orientation_discrepancy_p95_degrees_after_qc"])

    fig = fs.figure(width="double", height=4.65)
    grid = fig.add_gridspec(2, 1, height_ratios=[2.45, 1.25], hspace=0.44)
    fig.subplots_adjust(left=0.20, right=0.97, top=0.97, bottom=0.10)

    ax = fig.add_subplot(grid[0])
    y = np.arange(len(labels))
    ax.barh(y, values, height=0.52, color=fs.C["blue"])
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Post-QC usable measurements (%)")
    clean_axis(ax)
    for yi, value in zip(y, values):
        ax.text(value + 1.1, yi, f"{value:.1f}%", va="center", fontsize=fs.FONT["annot"])
    fs.panel(ax, "a", "Post-QC usability")

    ax2 = fig.add_subplot(grid[1])
    ax2.hlines(0, 0, p95, linewidth=1.6, color=fs.C["orange"])
    ax2.scatter([p95], [0], s=34, zorder=3, color=fs.C["orange"])
    ax2.set_yticks([0], ["Orientation"])
    ax2.set_xlim(0, max(10.0, p95 * 1.8))
    ax2.set_ylim(-0.65, 0.65)
    ax2.set_xlabel("Mirror discrepancy, p95 (°)")
    clean_axis(ax2)
    ax2.text(p95 + 0.25, 0, f"{p95:.2f}°", va="center", fontsize=fs.FONT["annot"])
    fs.panel(ax2, "b", "Orientation repeatability")

    save(fig, output_dir, "Figure_S1_v2_image_to_trait_technical_audit")


def figure_s2(summary: dict, output_dir: Path) -> None:
    audit = summary["perturbation_audit"]
    if int(audit["n_automated_perturbations"]) != 9:
        raise ValueError("Frozen audit no longer reports the expected nine perturbations")

    rhos = [
        float(audit["outline_half_resolution_spearman_rho_approx"]),
        float(audit["orientation_bbox_shift_5pct_spearman_rho"]),
    ]
    error = float(audit["orientation_bbox_shift_5pct_error_p95_degrees"])

    fig = fs.figure(width="double", height=4.45)
    grid = fig.add_gridspec(2, 1, height_ratios=[2.05, 1.25], hspace=0.48)
    fig.subplots_adjust(left=0.24, right=0.97, top=0.97, bottom=0.11)

    ax1 = fig.add_subplot(grid[0])
    labels = ["Outline\nhalf resolution", "Orientation\n5% bbox shift"]
    y = np.arange(2)
    ax1.barh(y, rhos, height=0.50, color=[fs.C["green"], fs.C["orange"]])
    ax1.set_yticks(y, labels)
    ax1.invert_yaxis()
    ax1.set_xlim(0, 1)
    ax1.set_xlabel("Spearman rank agreement with unperturbed measurement")
    clean_axis(ax1)
    for yi, rho in zip(y, rhos):
        prefix = "≈" if yi == 0 else ""
        ax1.text(rho + 0.02, yi, f"{prefix}{rho:.3g}", va="center", fontsize=fs.FONT["annot"])
    fs.panel(ax1, "a", "Representative retained perturbation metrics")

    ax2 = fig.add_subplot(grid[1])
    ax2.hlines(0, 0, error, linewidth=1.6, color=fs.C["orange"])
    ax2.scatter([error], [0], s=34, zorder=3, color=fs.C["orange"])
    ax2.set_yticks([0], ["Orientation\n5% bbox shift"])
    ax2.set_xlim(0, 65)
    ax2.set_ylim(-0.65, 0.65)
    ax2.set_xlabel("Angular discrepancy, p95 (°)")
    clean_axis(ax2)
    ax2.text(error + 1.0, 0, f"{error:.1f}°", va="center", fontsize=fs.FONT["annot"])
    fs.panel(ax2, "b", "Bounding-box displacement sensitivity")

    save(fig, output_dir, "Figure_S2_v2_image_to_trait_perturbation_audit")


def main() -> None:
    args = parse_args()
    fs.use(grid=False)
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    if summary.get("human_ground_truth_used") is not False:
        raise ValueError("This figure route is defined for the fully automated audit only")

    figure_s1(summary, args.output_dir)
    figure_s2(summary, args.output_dir)

    print(json.dumps({
        "status": "PASS",
        "summary": str(args.summary),
        "output_dir": str(args.output_dir),
        "figure_s2_scope": (
            "The frozen summary records nine perturbations but preserves detailed numerical metrics "
            "for representative outline and orientation scenarios only; unavailable scenario-level "
            "values are not reconstructed."
        ),
        "claim_boundary": summary["claim_boundary"],
    }, indent=2))


if __name__ == "__main__":
    main()
