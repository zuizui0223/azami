#!/usr/bin/env python3
"""Rebuild the automated image-to-trait technical-audit Supporting figures.

This script is intentionally limited to the frozen audit summary. It does not
reconstruct unavailable per-image perturbation records or treat human labels as
truth. It visualizes only results fixed in
analysis/ch1/image_to_trait_automated_technical_audit_summary.json.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/ch1/image_to_trait_automated_technical_audit_summary.json"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("manuscript/figures/v2_submission"),
    )
    return p.parse_args()


def save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fixed_time = datetime(2026, 9, 1, tzinfo=timezone.utc)
    fig.savefig(
        output_dir / f"{stem}.pdf",
        bbox_inches="tight",
        facecolor="white",
        metadata={
            "Title": stem,
            "Author": "",
            "Creator": "Azami Chapter 1 automated technical-audit figure builder",
            "CreationDate": fixed_time,
            "ModDate": fixed_time,
        },
    )
    plt.close(fig)


def clean_axis(ax: plt.Axes, *, grid_axis: str = "x") -> None:
    ax.grid(axis=grid_axis, linewidth=0.6, alpha=0.28)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)


def panel_label(fig: plt.Figure, y: float, label: str) -> None:
    # Fixed figure-coordinate x keeps every panel letter on the same far-left edge.
    fig.text(0.015, y, label, ha="left", va="top", fontsize=12, weight="bold")


def figure_s1(summary: dict, output_dir: Path) -> None:
    audit = summary["mirror_audit"]
    usable = audit["usable_fraction"]
    labels = ["Colour", "Gross outline", "Orientation"]
    values = np.array([usable["colour"], usable["outline"], usable["orientation"]]) * 100
    p95 = float(audit["orientation_discrepancy_p95_degrees_after_qc"])

    fig = plt.figure(figsize=(6.6, 5.8))
    grid = fig.add_gridspec(2, 1, height_ratios=[2.6, 1.35], hspace=0.42)
    fig.subplots_adjust(left=0.22, right=0.96, top=0.96, bottom=0.10)

    ax = fig.add_subplot(grid[0])
    y = np.arange(len(labels))
    ax.barh(y, values, height=0.52)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Post-QC usable measurements (%)")
    clean_axis(ax)
    for yi, value in zip(y, values):
        ax.text(value + 1.1, yi, f"{value:.1f}%", va="center", fontsize=9)
    panel_label(fig, 0.965, "(a)")

    ax2 = fig.add_subplot(grid[1])
    ax2.hlines(0, 0, p95, linewidth=2.0)
    ax2.scatter([p95], [0], s=48, zorder=3)
    ax2.set_yticks([0], ["Orientation"])
    ax2.set_xlim(0, max(10.0, p95 * 1.8))
    ax2.set_ylim(-0.65, 0.65)
    ax2.set_xlabel("Mirror discrepancy, p95 (°)")
    clean_axis(ax2)
    ax2.text(p95 + 0.25, 0, f"{p95:.2f}°", va="center", fontsize=10)
    panel_label(fig, 0.405, "(b)")

    save(fig, output_dir, "Figure_S1_v2_image_to_trait_technical_audit")


def figure_s2(summary: dict, output_dir: Path) -> None:
    audit = summary["perturbation_audit"]
    rhos = [
        float(audit["outline_half_resolution_spearman_rho_approx"]),
        float(audit["orientation_bbox_shift_5pct_spearman_rho"]),
    ]
    error = float(audit["orientation_bbox_shift_5pct_error_p95_degrees"])

    fig = plt.figure(figsize=(6.6, 5.4))
    grid = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.35], hspace=0.48)
    fig.subplots_adjust(left=0.25, right=0.96, top=0.96, bottom=0.11)

    ax1 = fig.add_subplot(grid[0])
    labels = ["Outline\nhalf resolution", "Orientation\n5% bbox shift"]
    y = np.arange(2)
    ax1.barh(y, rhos, height=0.50)
    ax1.set_yticks(y, labels)
    ax1.invert_yaxis()
    ax1.set_xlim(0, 1)
    ax1.set_xlabel("Spearman rank agreement with unperturbed measurement")
    clean_axis(ax1)
    for yi, rho in zip(y, rhos):
        prefix = "≈" if yi == 0 else ""
        ax1.text(rho + 0.02, yi, f"{prefix}{rho:.3g}", va="center", fontsize=9)
    panel_label(fig, 0.965, "(a)")

    ax2 = fig.add_subplot(grid[1])
    ax2.hlines(0, 0, error, linewidth=2.0)
    ax2.scatter([error], [0], s=48, zorder=3)
    ax2.set_yticks([0], ["Orientation\n5% bbox shift"])
    ax2.set_xlim(0, 65)
    ax2.set_ylim(-0.65, 0.65)
    ax2.set_xlabel("Angular discrepancy, p95 (°)")
    clean_axis(ax2)
    ax2.text(error + 1.0, 0, f"{error:.1f}°", va="center", fontsize=10)
    panel_label(fig, 0.405, "(b)")

    save(fig, output_dir, "Figure_S2_v2_image_to_trait_perturbation_audit")


def main() -> None:
    args = parse_args()
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    if summary.get("human_ground_truth_used") is not False:
        raise ValueError("This figure route is defined for the fully automated audit only")
    figure_s1(summary, args.output_dir)
    figure_s2(summary, args.output_dir)
    print(json.dumps({
        "status": "PASS",
        "summary": str(args.summary),
        "outputs": [
            str(args.output_dir / "Figure_S1_v2_image_to_trait_technical_audit.png"),
            str(args.output_dir / "Figure_S1_v2_image_to_trait_technical_audit.pdf"),
            str(args.output_dir / "Figure_S2_v2_image_to_trait_perturbation_audit.png"),
            str(args.output_dir / "Figure_S2_v2_image_to_trait_perturbation_audit.pdf"),
        ],
        "claim_boundary": summary["claim_boundary"],
    }, indent=2))


if __name__ == "__main__":
    main()
