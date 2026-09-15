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
ESTIMATOR = REFERENCE / "estimator_validity" / "rv_estimator_validity_summary.json"
ESTIMATOR_SHA256 = "081b37b798de55cc6eb6bb2948b6737a173644fac3ea411505a617857704951e"
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
    estimator_actual = sha256(ESTIMATOR)
    if estimator_actual != ESTIMATOR_SHA256:
        raise ValueError(
            f"estimator-validity summary: SHA-256 mismatch {estimator_actual}; expected {ESTIMATOR_SHA256}"
        )
    receipt.append({
        "path": "estimator_validity/rv_estimator_validity_summary.json",
        "sha256": estimator_actual,
        "run": 34933875866,
        "artifact": 10382387052,
        "artifact_sha256": "345866a7e333f78677ad3797811e2cbe82d3e09ee061f7642df7f4c4d5ec008e",
    })
    return receipt


def load_inputs() -> dict:
    sources = source_receipt()
    within = pd.read_csv(UPGRADE / "complete18_construct_integration_within.csv", index_col=0)
    among = pd.read_csv(UPGRADE / "complete18_construct_integration_among.csv", index_col=0)
    pairwise = pd.read_csv(UPGRADE / "complete18_construct_pairwise.csv")
    bootstrap = pd.read_csv(UPGRADE / "complete18_taxon_bootstrap.csv")
    summary = json.loads((UPGRADE / "construct_scale_contrast_summary.json").read_text(encoding="utf-8"))
    upgrade = json.loads((UPGRADE / "construct_scale_upgrade_report.json").read_text(encoding="utf-8"))
    estimator = json.loads(ESTIMATOR.read_text(encoding="utf-8"))

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
    if estimator["common_cohort"] != {
        "observations": 1734,
        "taxa": 42,
        "constructs": 9,
        "relations": 36,
    }:
        raise ValueError("estimator-validity cohort differs from the frozen common cohort")
    if not estimator["equal_n_within_resampling"]["overall_strength_gate_pass"]:
        raise ValueError("estimator-validity gate does not support stronger-overall-among wording")

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
        "estimator": estimator,
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
    summary = data["summary"]
    upgrade = data["upgrade"]
    estimator = data["estimator"]

    fs.use(grid=False)
    fig = fs.figure(width="double", height=6.45)
    gs = fig.add_gridspec(2, 2, left=0.15, right=0.92, bottom=0.09, top=0.95, wspace=0.38, hspace=0.38)
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
    equal_n = estimator["equal_n_within_resampling"]
    note_box = dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.4)
    ax_c.text(
        0.03,
        0.97,
        f"ρ = {alignment['rho']:.3f}; QAP P = {alignment['qap_p_one_sided']:.4f}\n"
        f"raw: {summary['observed']['relations_stronger_among']}/36 among > within\n"
        f"equal-n median: {equal_n['relations_stronger_among_median']:.0f}/36",
        transform=ax_c.transAxes,
        ha="left",
        va="top",
        fontsize=fs.FONT["annot"],
        bbox=note_box,
    )
    ax_c.legend(loc="lower right", fontsize=fs.FONT["footnote"])

    # Explicit estimator-validity panel. The raw taxon bootstrap is retained as
    # descriptive uncertainty from the frozen estimator; the equal-n row asks
    # whether the strength direction persists when both scales use 42 rows.
    raw_boot = summary["taxon_bootstrap"]
    rows = [
        (
            "Raw taxon\nbootstrap",
            raw_boot["difference_of_median_rv_median"],
            raw_boot["difference_of_median_rv_low95"],
            raw_boot["difference_of_median_rv_high95"],
            raw_boot["probability_median_among_exceeds_within"],
        ),
        (
            "Equal-n\nsensitivity",
            equal_n["among_minus_within_median_rv_median"],
            equal_n["among_minus_within_median_rv_low95"],
            equal_n["among_minus_within_median_rv_high95"],
            equal_n["probability_positive_median_difference"],
        ),
    ]
    for y_pos, (label, centre, low_ci, high_ci, probability) in enumerate(rows[::-1]):
        xerr = np.array([[centre - low_ci], [high_ci - centre]])
        ax_d.errorbar(
            centre,
            y_pos,
            xerr=xerr,
            fmt="o",
            capsize=3,
            linewidth=1.2,
            markersize=5,
        )
        ax_d.text(
            high_ci + 0.004,
            y_pos,
            f"P(Δ>0)={probability:.3f}",
            va="center",
            ha="left",
            fontsize=fs.FONT["annot"],
        )
    ax_d.axvline(0, linestyle="--", linewidth=0.8, color=fs.C["rule"])
    ax_d.set_yticks([0, 1], [rows[1][0], rows[0][0]])
    ax_d.set_xlabel("Δ median RV (among − within)")
    ax_d.set_title("(d) Estimator-validity contrast", loc="left", fontsize=fs.FONT["panel"], fontweight="bold", pad=4)
    ax_d.text(
        0.03,
        0.04,
        "Equal-n: one centred observation per taxon (42 rows)",
        transform=ax_d.transAxes,
        ha="left",
        va="bottom",
        fontsize=fs.FONT["footnote"],
        bbox=note_box,
    )
    ax_d.set_ylim(-0.55, 1.55)

    out_dir.mkdir(parents=True, exist_ok=True)
    written = fs.savefig(fig, STEM, width="double", outdir=out_dir)
    plt.close(fig)

    provenance = {
        "schema_version": 2,
        "figure_stem": STEM,
        "figure_role": "current Chapter 1 cross-scale integration main-figure candidate",
        "analysis_id": summary["analysis_id"],
        "estimator_validity_analysis_id": estimator["analysis_id"],
        "common_cohort": upgrade["common_cohort"],
        "construct_order": list(CORE),
        "construct_modules": MODULES,
        "source_files": data["sources"],
        "display_transforms": {
            "matrix_colour_norm": "PowerNorm(gamma=0.45) on untransformed RV values; shared across panels a-b",
            "relation_scatter_axes": "log-log display only; statistics use untransformed RV values",
            "estimator_validity_axis": "linear point-interval display of frozen raw bootstrap and equal-n summary values",
        },
        "panels": {
            "a": "within-taxon 9 x 9 construct integration matrix",
            "b": "among-taxon 9 x 9 construct integration matrix on the identical cohort",
            "c": "36 raw relation-wise within versus among RV values with 1:1 reference and equal-n median relation-count context",
            "d": "raw taxon-bootstrap and equal-n estimator-validity summaries of among-minus-within median RV",
        },
        "headline_metrics": {
            "matrix_alignment_rho": upgrade["common_cohort_matrix_alignment"]["rho"],
            "matrix_alignment_qap_p": upgrade["common_cohort_matrix_alignment"]["qap_p_one_sided"],
            "raw_median_within_rv": summary["observed"]["median_within_rv"],
            "raw_median_among_rv": summary["observed"]["median_among_rv"],
            "raw_relations_stronger_among": summary["observed"]["relations_stronger_among"],
            "equal_n_difference_median": equal_n["among_minus_within_median_rv_median"],
            "equal_n_difference_low95": equal_n["among_minus_within_median_rv_low95"],
            "equal_n_difference_high95": equal_n["among_minus_within_median_rv_high95"],
            "equal_n_probability_positive": equal_n["probability_positive_median_difference"],
            "equal_n_relations_stronger_among_median": equal_n["relations_stronger_among_median"],
        },
        "claim_boundary": (
            "raw RV values remain frozen descriptive outputs; equal-n sensitivity supports only the qualitative "
            "statement that visible-phenotype integration is stronger overall among taxa. Raw 33/36 is not "
            "treated as an estimator-invariant count. No genetic, developmental, functional or causal modularity "
            "is inferred."
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
