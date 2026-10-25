#!/usr/bin/env python3
"""Summarize the predeclared WCVP taxonomy sensitivity gates.

Pass/fail criteria are declared here before the authority-remapped outcome is
run. They test whether the current GEB headline survives a deterministic change
from source-assigned taxon labels to uniquely resolved accepted WCVP keys.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ANCHORS = (
    ("floral_chroma", "chelsa_rsds_mean", -1),
    ("presentation_angle", "chelsa_bio12", +1),
)


def scalar_anchor(frame: pd.DataFrame, construct: str, predictor: str) -> dict:
    row = frame[(frame.construct_id == construct) & (frame.predictor == predictor)]
    if len(row) != 1:
        raise ValueError(f"expected one row for {construct} x {predictor}; got {len(row)}")
    r = row.iloc[0]
    if str(r["kind"]) != "scalar":
        raise ValueError(f"anchor {construct} x {predictor} is not scalar")
    return {
        "construct": construct,
        "predictor": predictor,
        "n_taxa": int(r["n_taxa"]),
        "beta_std": float(r["beta_std"]),
        "p_value": float(r["p_value"]),
        "q_bh": float(r["q_bh"]),
    }


def summarize(
    prep_report: Path,
    authority_axes: Path,
    authority_upgrade: Path,
    authority_contrast: Path,
    frozen_axes: Path,
    out: Path,
) -> dict:
    prep = json.loads(prep_report.read_text(encoding="utf-8"))
    axes = pd.read_csv(authority_axes)
    frozen = pd.read_csv(frozen_axes)
    upgrade = json.loads(authority_upgrade.read_text(encoding="utf-8"))
    contrast = json.loads(authority_contrast.read_text(encoding="utf-8"))

    authority_anchors = []
    frozen_anchors = []
    anchor_passes = []
    for construct, predictor, expected_sign in ANCHORS:
        current = scalar_anchor(axes, construct, predictor)
        baseline = scalar_anchor(frozen, construct, predictor)
        sign_ok = current["beta_std"] * expected_sign > 0
        fdr_ok = current["q_bh"] < 0.05
        current["expected_sign"] = expected_sign
        current["sign_retained"] = bool(sign_ok)
        current["bh_q_lt_0_05"] = bool(fdr_ok)
        current["gate_pass"] = bool(sign_ok and fdr_ok)
        authority_anchors.append(current)
        frozen_anchors.append(baseline)
        anchor_passes.append(current["gate_pass"])

    alignment = upgrade["common_cohort_matrix_alignment"]
    bootstrap_alignment = upgrade["taxon_bootstrap"]
    geometry_pass = (
        float(alignment["rho"]) > 0
        and float(alignment["qap_p_one_sided"]) < 0.05
        and float(bootstrap_alignment["rho_low95"]) > 0
    )

    observed = contrast["observed"]
    boot = contrast["taxon_bootstrap"]
    strength_pass = (
        float(boot["difference_of_median_rv_low95"]) > 0
        and int(observed["relations_stronger_among"]) > int(observed["relations"]) / 2
    )

    result = {
        "analysis_id": "ch1_v3_wcvp_authority_taxonomy_sensitivity_summary_20260912",
        "gate_declared_before_outcome": True,
        "taxonomy_preparation": prep,
        "authority_common_cohort": upgrade["common_cohort"],
        "authority_geometry": {
            "rho": float(alignment["rho"]),
            "qap_p_one_sided": float(alignment["qap_p_one_sided"]),
            "bootstrap_rho_median": float(bootstrap_alignment["rho_median"]),
            "bootstrap_rho_low95": float(bootstrap_alignment["rho_low95"]),
            "bootstrap_rho_high95": float(bootstrap_alignment["rho_high95"]),
            "gate": "rho > 0 AND QAP P < 0.05 AND bootstrap rho lower 95% bound > 0",
            "gate_pass": bool(geometry_pass),
        },
        "authority_strength": {
            "median_within_rv": float(observed["median_within_rv"]),
            "median_among_rv": float(observed["median_among_rv"]),
            "relations_stronger_among": int(observed["relations_stronger_among"]),
            "relations_total": int(observed["relations"]),
            "bootstrap_difference_median": float(boot["difference_of_median_rv_median"]),
            "bootstrap_difference_low95": float(boot["difference_of_median_rv_low95"]),
            "bootstrap_difference_high95": float(boot["difference_of_median_rv_high95"]),
            "bootstrap_probability_among_exceeds_within": float(boot["probability_median_among_exceeds_within"]),
            "gate": "bootstrap lower 95% bound of among-minus-within median RV > 0 AND majority of observed relations stronger among taxa",
            "gate_pass": bool(strength_pass),
        },
        "frozen_source_assigned_anchors": frozen_anchors,
        "authority_resolved_anchors": authority_anchors,
        "anchor_gate": "each frozen headline anchor retains its predeclared sign AND BH q < 0.05 under authority-resolved grouping",
        "anchor_gate_pass": bool(all(anchor_passes)),
        "headline_taxonomy_robust": bool(geometry_pass and strength_pass and all(anchor_passes)),
        "claim_boundary": (
            "sensitivity to accepted-name resolution and synonym collapse only; a pass does not establish identification accuracy, "
            "species boundaries, genetic lineages or causal evolutionary mechanism"
        ),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prep-report", type=Path, required=True)
    parser.add_argument("--authority-axes", type=Path, required=True)
    parser.add_argument("--authority-upgrade", type=Path, required=True)
    parser.add_argument("--authority-contrast", type=Path, required=True)
    parser.add_argument("--frozen-axes", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(
        args.prep_report,
        args.authority_axes,
        args.authority_upgrade,
        args.authority_contrast,
        args.frozen_axes,
        args.out,
    ), indent=2))


if __name__ == "__main__":
    main()
