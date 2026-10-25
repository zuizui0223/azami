#!/usr/bin/env python3
"""Recreate manuscript summary figures from the locked V4 S1/S2 CSV subsets.

These plots show frozen aggregate statistics only. For the complete 208-test
ledger, 28-edge RV matrices and spatial sampling maps, retrieve and verify
the source artifact described in CLAIM_LEDGER_AND_SUBMISSION_GATES.md.
Nothing here is a fresh regression, independent validation or mechanistic test.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def check_locked_input(assoc: list[dict[str, str]], module: list[dict[str, str]]) -> None:
    """Fail closed if the frozen-display inputs were accidentally replaced."""
    expected = [
        ("S1", "BIO18", 168907, 0.045447, 0.00920635),
        ("S2", "ELEVATION", 26171, -0.173732, 0.00131233),
        ("S2", "RSDS", 27200, 0.134504, 0.01456813),
        ("S2", "RSDS", 27200, -0.142801, 0.01925414),
    ]
    assert len(assoc) == len(expected), "Holm-display associations must remain four"
    for row, (s, e, n, b, p) in zip(assoc, expected):
        assert row["stratum"] == s and row["exposure"] == e
        assert int(row["n"]) == n
        assert np.isclose(float(row["beta"]), b, atol=1e-10)
        assert np.isclose(float(row["holm_p"]), p, atol=1e-10)
        assert float(row["ci_low"]) < float(row["beta"]) < float(row["ci_high"])
    expected_rv = [0.085245, 0.001419, 0.197455, 0.015453]
    assert len(module) == 4, "Four RV medians required"
    assert all(np.isclose(float(x["median_RV"]), v, atol=1e-10) for x, v in zip(module, expected_rv))


def forest(assoc: list[dict[str, str]], dest: Path) -> None:
    labels = [
        "S1 among taxa  |  Head image angle × BIO18",
        "S2 introduced arvense  |  Lightness × elevation",
        "S2 introduced vulgare  |  Chroma × radiation",
        "S2 introduced vulgare  |  Lightness × radiation",
    ]
    n = len(labels)
    y = np.arange(n)[::-1]
    b = np.array([float(r["beta"]) for r in assoc])
    lo = np.array([float(r["ci_low"]) for r in assoc])
    hi = np.array([float(r["ci_high"]) for r in assoc])
    fig, ax = plt.subplots(figsize=(11, 5.1), constrained_layout=True)
    ax.axvline(0, ls="--", lw=1.2, color="0.4")
    ax.errorbar(b, y, xerr=np.vstack([b - lo, hi - b]), fmt="o",
                lw=1.7, capsize=4, markersize=6)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Standardized coefficient β (pointwise 95% CI)")
    ax.set_xlim(-0.29, 0.23)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_title("Four Holm-supported associations in frozen V4 S1–S2")
    for i, row in enumerate(assoc):
        ax.text(0.229, y[i], f"Holm P = {float(row['holm_p']):.4g}  |  n = {row['n']}",
                va="center", ha="right", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})
    fig.text(.012, .01,
             "The four points are the passing display subset of 80 S1 and 128 S2 biological tests. "
             "CIs are pointwise, not Holm-simultaneous. Independent mechanisms are not established.",
             fontsize=8)
    for ext in ("png", "svg"):
        fig.savefig(dest / f"figure2_holm_forest.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def rv_medians(module: list[dict[str, str]], dest: Path) -> None:
    sample = {(r["scale"], r["pair_type"]): float(r["median_RV"]) for r in module}
    xs = np.array([0, 1])
    fig, ax = plt.subplots(figsize=(7.4, 4.5), constrained_layout=True)
    for scale, name in [("within_taxa", "Within taxa"), ("among_taxa", "Among taxa")]:
        values = [sample[(scale, "same_domain")], sample[(scale, "cross_domain")]]
        ax.plot(xs, values, marker="o", markersize=7, linewidth=2, label=name)
        for x, y in zip(xs, values):
            ax.annotate(f"{y:.4f}", (x, y), xytext=(0, 7), textcoords="offset points",
                        ha="center", fontsize=9)
    ax.set_xticks(xs, ["Pairs within domain", "Pairs across domains"])
    ax.set_ylabel("Median RV (descriptive magnitude)")
    ax.set_ylim(0, 0.235)
    ax.set_title("Photographic trait-domain association geometry")
    ax.legend(loc="upper right", frameon=False)
    fig.text(.01, .01,
             "Same 156,683 S1 photographs / 174 taxa; BIO18, BIO05, RSDS, elevation adjusted. "
             "Shared pixels and RV baseline differences remain.",
             fontsize=8)
    for ext in ("png", "svg"):
        fig.savefig(dest / f"figure3_module_RV_medians.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "rendered")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    assoc = rows(ROOT / "table1_holm_associations.csv")
    module = rows(ROOT / "module_summary.csv")
    check_locked_input(assoc, module)
    forest(assoc, args.output_dir)
    rv_medians(module, args.output_dir)
    print(f"Locked-source figures created in {args.output_dir}")


if __name__ == "__main__":
    main()
