#!/usr/bin/env python3
"""Render the current Chapter 1 scale-dependent integration figure.

This is a presentation-only renderer over frozen current-reference outputs. It
must not refit models, resample observations, change construct definitions or
search for additional ecological associations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import numpy as np
import pandas as pd

from analysis.v3.run_construct_scale_upgrade import CORE, MODULES
from legacy.v2.analysis import azami_figstyle as fs

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "reproducibility" / "current_reference"
UPGRADE = REFERENCE / "upgrade"
STEM = "Figure_v3_scale_integration"

SOURCE_PATHS = (
    "upgrade/complete18_construct_integration_within.csv",
    "upgrade/complete18_construct_integration_among.csv",
    "upgrade/complete18_construct_pairwise.csv",
    "upgrade/complete18_taxon_bootstrap.csv",
    "upgrade/construct_scale_contrast_summary.json",
    "upgrade/construct_scale_upgrade_report.json",
)

LABELS = {
    "presentation_angle": "angle",
    "floral_lightness": "lightness",
    "floral_chroma": "chroma",
    "floral_hue": "hue",
    "head_elongation": "elongation",
    "head_compactness": "compactness",
    "involucre_form": "involucre",
    "projection_prominence": "projection\nprominence",
    "projection_pattern": "projection\npattern",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_receipt() -> list[dict]:
    manifest = json.loads((REFERENCE / "manifest.json").read_text(encoding="utf-8"))
    lut = {row["path"]: row for row in manifest["files"]}
    receipt = []
    for rel in SOURCE_PATHS:
        if rel not in lut:
            raise KeyError(f"current-reference manifest does not contain {rel}")
        path = REFERENCE / rel
        actual = sha256(path)
        expected = lut[rel]["sha256"]
        if actual != expected:
            raise ValueError(f"{rel}: SHA-256 mismatch {actual}; expected {expected}")
        receipt.append({"path": rel, "sha256": actual})
    return receipt


def load_inputs() -> dict:
    sources = source_receipt()
    within = pd.read_csv(UPGRADE / "complete18_construct_integration_within.csv", index_col=0)
    among = pd.read_csv(UPGRADE / "complete18_construct_integration_among.csv", index_col=0)
    pairwise = pd.read_csv(UPGRADE / "complete18_construct_pairwise.csv")
    bootstrap = pd.read_csv(UPGRADE / "complete18_taxon_bootstrap.csv")
    summary = json.loads((UPGRADE / "construct_scale_contrast_summary.json").read_text(encoding="utf-8"))
    upgrade = json.loads((UPGRADE / "construct_scale_upgrade_report.json").read_text(encoding="utf-8"))

    expected = list(CORE)
    for name, matrix in (("within", within), ("among", among)):
        if list(matrix.index) != expected or list(matrix.columns) != expected:
            raise ValueError(f"{name} matrix construct order differs from frozen CORE")
        if matrix.shape != (9, 9):
            raise ValueError(f"{name} matrix must be 9 x 9")
        if not np.allclose(np.diag(matrix.to_numpy(float)), 1.0):
            raise ValueError(f"{name} matrix diagonal is not one")

    if len(pairwise) != 36 or len(bootstrap) != 1000:
        raise ValueError("expected 36 construct relations and 1,000 bootstrap replicates")
    if summary["observed"]["relations"] != 36:
        raise ValueError("contrast summary relation count mismatch")
    if upgrade["common_cohort"] != {
        "observations": 1734,
        "taxa": 42,
        "minimum_complete_observations_per_taxon": 5,
    }:
        raise ValueError("common cohort differs from frozen 1,734-observation / 42-taxon cohort")

    upper = np.triu_indices(len(CORE), 1)
    within_upper = within.to_numpy(float)[upper]
    among_upper = among.to_numpy(float)[upper]
    observed_delta = float(np.median(among_upper) - np.median(within_upper))
    if not np.isclose(observed_delta, summary["observed"]["difference_of_medians_among_minus_within"]):
        raise ValueError("matrix-derived observed contrast differs from frozen summary")
    if int(np.sum(among_upper > within_upper)) != summary["observed"]["relations_stronger_among"]:
        raise ValueError("relation-strength count differs from frozen summary")

    bootstrap = bootstrap.copy()
    bootstrap["delta_median_rv"] = bootstrap["median_among_rv"] - bootstrap["median_within_rv"]
    return {
        "within": within,
        "among": among,
        "pairwise": pairwise,
        "bootstrap": bootstrap,
        "summary": summary,
        "upgrade": upgrade,
        "sources": sources,
    }


def _matrix_panel(ax, matrix: pd.DataFrame, title: str, norm: PowerNorm):
    arr = matrix.to_numpy(float).copy()
    np.fill_diagonal(arr, np.nan)
    image = ax.imshow(arr, cmap="viridis", norm=norm, interpolation="nearest")
    ticks = np.arange(len(CORE))
    labels = [LABELS[key] for key in CORE]
    ax.set_xticks(ticks, labels, rotation=55, ha="right", rotation_mode="anchor")
    ax.set_yticks(ticks, labels)
    ax.tick_params(length=0)
    # Module boundaries: presentation | colour | head form | involucre/armature.
    for boundary in (0.5, 3.5, 5.5):
        ax.axhline(boundary, color="white", linewidth=0.8)
        ax.axvline(boundary, color="white", linewidth=0.8)
    ax.set_title(title, loc="left", fontsize=fs.FONT["panel"], fontweight="bold", pad=4)
    return image


def render(out_dir: Path) -> dict:
    data = load_inputs()
    within = data["within"]
    among = data["among"]
    pairwise = data["pairwise"]
    bootstrap = data["bootstrap"]
    summary = data["summary"]
    upgrade = data["upgrade"]

    fs.use(grid=False)
    fig = fs.figure(width="double", height=6.45)
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.93, bottom=0.09, top=0.95, wspace=0.34, hspace=0.38)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    upper = np.triu_indices(len(CORE), 1)
    max_rv = max(
        float(np.max(within.to_numpy(float)[upper])),
        float(np.max(among.to_numpy(float)[upper])),
    )
    norm = PowerNorm(gamma=0.45, vmin=0.0, vmax=max_rv)
    image = _matrix_panel(ax_a, within, "(a) Within taxa", norm)
    _matrix_panel(ax_b, among, "(b) Among taxa", norm)
    cbar = fig.colorbar(image, ax=[ax_a, ax_b], orientation="vertical", fraction=0.035, pad=0.025, aspect=24)
    cbar.set_label("Pairwise RV integration", rotation=270, labelpad=11)

    # Direct relation-by-relation scale contrast.
    x = pairwise["within_taxon_rv"].to_numpy(float)
    y = pairwise["among_taxon_rv"].to_numpy(float)
    positive = y > x
    ax_c.scatter(x[~positive], y[~positive], s=18, alpha=0.85, label="within ≥ among")
    ax_c.scatter(x[positive], y[positive], s=18, alpha=0.85, label="among > within")
    low = min(float(np.min(x)), float(np.min(y))) / 1.8
    high = max(float(np.max(x)), float(np.max(y))) * 1.35
    ax_c.plot([low, high], [low, high], linestyle="--", linewidth=0.8, color=fs.C["rule"], zorder=0)
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlim(low, high)
    ax_c.set_ylim(low, high)
    ax_c.set_xlabel("Within-taxon RV")
    ax_c.set_ylabel("Among-taxon RV")
    ax_c.set_title("(c) Relation-wise scale contrast", loc="left", fontsize=fs.FONT["panel"], fontweight="bold", pad=4)
    alignment = upgrade["common_cohort_matrix_alignment"]
    ax_c.text(
        0.03,
        0.97,
        f"ρ = {alignment['rho']:.3f}; QAP P = {alignment['qap_p_one_sided']:.4f}\n"
        f"{summary['observed']['relations_stronger_among']}/36 relations stronger among taxa",
        transform=ax_c.transAxes,
        ha="left",
        va="top",
        fontsize=fs.FONT["annot"],
    )
    ax_c.legend(loc="lower right", fontsize=fs.FONT["footnote"])

    # Bootstrap uncertainty in the direct among-minus-within strength contrast.
    delta = bootstrap["delta_median_rv"].to_numpy(float)
    ax_d.hist(delta, bins=24, edgecolor="white", linewidth=0.4)
    observed = summary["observed"]["difference_of_medians_among_minus_within"]
    ax_d.axvline(0, linestyle="--", linewidth=0.8, color=fs.C["rule"])
    ax_d.axvline(observed, linewidth=1.2, color=fs.C["black"], label="observed")
    ax_d.set_xlabel("Bootstrap Δ median RV (among − within)")
    ax_d.set_ylabel("Replicates")
    ax_d.set_title("(d) Taxon-bootstrap contrast", loc="left", fontsize=fs.FONT["panel"], fontweight="bold", pad=4)
    boot = summary["taxon_bootstrap"]
    ax_d.text(
        0.97,
        0.97,
        f"median = {boot['difference_of_median_rv_median']:+.3f}\n"
        f"95% interval {boot['difference_of_median_rv_low95']:+.3f} to {boot['difference_of_median_rv_high95']:+.3f}\n"
        f"P(Δ > 0) = {boot['probability_median_among_exceeds_within']:.3f}",
        transform=ax_d.transAxes,
        ha="right",
        va="top",
        fontsize=fs.FONT["annot"],
    )
    ax_d.legend(loc="upper left", fontsize=fs.FONT["footnote"])

    out_dir.mkdir(parents=True, exist_ok=True)
    written = fs.savefig(fig, STEM, width="double", outdir=out_dir)
    plt.close(fig)

    provenance = {
        "schema_version": 1,
        "figure_stem": STEM,
        "figure_role": "current Chapter 1 cross-scale integration main-figure candidate",
        "analysis_id": summary["analysis_id"],
        "common_cohort": upgrade["common_cohort"],
        "construct_order": list(CORE),
        "construct_modules": MODULES,
        "source_files": data["sources"],
        "panels": {
            "a": "within-taxon 9 x 9 construct integration matrix",
            "b": "among-taxon 9 x 9 construct integration matrix on the identical cohort",
            "c": "36 relation-wise within versus among RV values with 1:1 reference",
            "d": "1,000 taxon-bootstrap differences in median RV, among minus within",
        },
        "headline_metrics": {
            "matrix_alignment_rho": upgrade["common_cohort_matrix_alignment"]["rho"],
            "matrix_alignment_qap_p": upgrade["common_cohort_matrix_alignment"]["qap_p_one_sided"],
            "median_within_rv": summary["observed"]["median_within_rv"],
            "median_among_rv": summary["observed"]["median_among_rv"],
            "relations_stronger_among": summary["observed"]["relations_stronger_among"],
            "bootstrap_probability_among_exceeds_within": boot["probability_median_among_exceeds_within"],
        },
        "claim_boundary": (
            "presentation-only visualization of frozen image-defined construct integration; "
            "not evidence that evolution increases integration and not genetic, developmental, "
            "functional or causal modularity"
        ),
        "outputs": [{"path": path.name, "sha256": sha256(path), "size_bytes": path.stat().st_size} for path in written],
    }
    provenance_path = out_dir / f"{STEM}_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance["provenance_path"] = provenance_path.name
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("work/scale-integration-figure"))
    args = parser.parse_args()
    result = render(args.out_dir.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
