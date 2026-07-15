#!/usr/bin/env python3
"""Validate the frozen Chapter 1 manuscript claim registry.

This deliberately uses only the Python standard library so the submission contract
can be checked in a clean runner without scientific Python packages.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--claims", default="manuscript/final_claims.json")
    return parser.parse_args()


def main() -> None:
    claims = json.loads(Path(parse_args().claims).read_text(encoding="utf-8"))

    datasets = claims["datasets"]
    results = claims["primary_results"]
    inference = claims["within_species_inference_levels"]
    history = claims["historical_sensitivity"]
    molecular = claims["molecular_database_coverage"]

    assert datasets["continuous_species_comparison"]["n_taxa"] == 216
    assert datasets["complete_lability_cohort"]["n_taxa"] == 102
    assert sum(results["quadrants"].values()) == 102
    assert -1 <= results["lability_axis_spearman_rho"] <= 1

    strict = inference["strict_le_10km_primary_cohort"]
    expanded = inference["expanded_pooled_coefficient_table"]
    assert strict["n_bh_fdr_signals"] == 0
    assert expanded["n_bh_fdr_signals"] == 4
    assert "different cohorts" in inference["interpretation_rule"].lower()

    assert history["n_input_taxa"] == 216
    assert history["n_direct_backbone_taxa"] == 54
    assert history["n_random_trees"] == 50
    assert history["n_signal_fits"] == 636
    assert history["n_failed_signal_fits"] == 0
    assert history["direct_backbone_supported_non_circular_endpoints"] == 0

    for count in molecular.values():
        assert isinstance(count, int) and 0 <= count <= 216
    assert molecular["complete_plastome"] < molecular["plastid"]

    modules = results["module_medians"]
    assert modules["colour"]["environmental_responsiveness"] == max(
        value["environmental_responsiveness"] for value in modules.values()
    )
    assert modules["shape"]["within_variation"] == max(
        value["within_variation"] for value in modules.values()
    )

    print(
        json.dumps(
            {
                "status": "valid",
                "release": claims["analysis_release"],
                "complete_lability_taxa": datasets["complete_lability_cohort"]["n_taxa"],
                "strict_fdr_signals": strict["n_bh_fdr_signals"],
                "expanded_fdr_signals": expanded["n_bh_fdr_signals"],
                "phylogenetic_fits": history["n_signal_fits"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
