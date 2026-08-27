#!/usr/bin/env python3
"""Validate the submission-facing Azami Chapter 1 v2 claim registry."""
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
    assert claims["release_label"] == "Azami Chapter 1 v2"
    assert claims["freeze_tag"] == "azami-ch1-v2-2026-08-27"
    assert claims["status"] == "submission_hold_pending_independent_validation"
    assert Path(claims["defense_document"]).is_file()

    cohort = claims["cohort"]
    traits = claims["trait_contract"]
    scale = claims["scale_results"]
    sampling = claims["sampling_composition"]
    candidates = claims["adaptive_pattern_candidates_under_current_controls"]

    assert cohort["n_spatial_observations"] == 46276
    assert cohort["n_source_assigned_taxa"] == 259
    assert traits["n_registered_endpoints"] == 27
    assert traits["n_measured_endpoints"] == 22
    assert traits["n_unexecuted_endpoints"] == 5
    assert traits["n_joint_inferential_units"] == 26
    assert traits["raw_fraction_below_taxon_mean_range"] == [
        0.8137321906969928,
        0.9862520151073496,
    ]
    assert scale["primary_cross_scale_classes"] == {
        "both_scales": 7,
        "within_only": 18,
        "among_only": 3,
        "neither": 152,
        "not_comparable": 54,
    }
    assert sampling["n_scenario_rows"] == 674
    assert sampling["n_within_pairs_stable_all_declared_scenarios"] == 16
    assert sampling["n_among_pairs_stable_all_declared_scenarios"] == 10

    assert [row["endpoint"] for row in candidates] == [
        "corolla_lab_chroma",
        "orientation_image_vertical_angle",
    ]
    assert all(row["n_historical_placement_trees_passed"] == 52 for row in candidates)
    assert "not anthocyanin concentration" in candidates[0]["mechanism_boundary"]
    assert "not gravity" in candidates[1]["mechanism_boundary"]
    assert all(row["decisive_next_test"] for row in candidates)

    legacy_keys = {
        "datasets",
        "legacy_precision_audit",
        "revised_primary_results",
        "within_species_inference_levels",
        "auxiliary_involucre",
    }
    assert legacy_keys.isdisjoint(claims), "legacy numerical claim blocks leaked into v2 registry"
    assert claims["submission_blockers"], "submission HOLD must retain open external gates"

    print(json.dumps({
        "status": "valid",
        "release": claims["analysis_release"],
        "registered_endpoints": traits["n_registered_endpoints"],
        "measured_endpoints": traits["n_measured_endpoints"],
        "sampling_scenarios": sampling["n_scenario_rows"],
        "candidate_pairs": [f'{r["endpoint"]}~{r["predictor"]}' for r in candidates],
        "submission_status": claims["status"],
    }, indent=2))


if __name__ == "__main__":
    main()
