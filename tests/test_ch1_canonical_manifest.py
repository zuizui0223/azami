from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_canonical_v2_manifest_paths_exist() -> None:
    payload = json.loads((ROOT / "analysis/ch1/pipeline.json").read_text(encoding="utf-8"))
    assert payload["schema_version"] >= 18
    assert payload["policy"]["active_results"] == "full27_full_environment_only"
    for stage in payload["stages"].values():
        for key in (
            "script", "measurement_script", "validator", "spatial_script",
            "historical_script", "contract", "workflow", "result_ledger",
            "hypothesis_status",
        ):
            value = stage.get(key)
            if value:
                assert (ROOT / value).is_file(), value


def test_submission_manuscript_sources_exist() -> None:
    required = (
        "manuscript/00_title_abstract.md",
        "manuscript/01_introduction.md",
        "manuscript/02_methods.md",
        "manuscript/03_results.md",
        "manuscript/04_discussion.md",
        "manuscript/05_conclusion_and_declarations.md",
        "manuscript/06_references.md",
        "manuscript/07_data_code_availability.md",
        "manuscript/SUBMISSION_MANUSCRIPT.md",
        "manuscript/supplement/SUPPLEMENTARY_INFORMATION.md",
        "manuscript/final_claims.json",
    )
    for relative in required:
        assert (ROOT / relative).is_file(), relative


def test_active_claim_registry_has_no_legacy_numerical_blocks() -> None:
    claims = json.loads((ROOT / "manuscript/final_claims.json").read_text(encoding="utf-8"))
    forbidden = {
        "datasets", "legacy_precision_audit", "revised_primary_results",
        "within_species_inference_levels", "auxiliary_involucre",
    }
    assert forbidden.isdisjoint(claims)
