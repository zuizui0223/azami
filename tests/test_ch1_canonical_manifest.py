from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_canonical_pipeline_paths_exist() -> None:
    pipeline = ROOT / "analysis" / "ch1" / "pipeline.json"
    payload = json.loads(pipeline.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 7
    assert payload["policy"]["move_numbered_scripts"] is False

    path_keys = (
        "script",
        "selection_script",
        "annotation_app",
        "finalizer",
        "evaluation_script",
        "measurement_script",
        "join_script",
        "residual_export_script",
        "region_script",
        "input_builder",
        "audit_script",
        "workflow",
        "evaluation_workflow",
        "audit_workflow",
    )
    for stage in payload["stages"].values():
        for key in path_keys:
            value = stage.get(key)
            if value:
                assert (ROOT / value).is_file(), value


def test_submission_controls_exist() -> None:
    required = (
        "README.md",
        "manuscript/SUBMISSION_MANUSCRIPT.md",
        "manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md",
        "manuscript/FIGURE_TABLE_MAP.md",
        "manuscript/final_claims.json",
        "manuscript/RUNBOOK.md",
        "manuscript/supplement/SUPPLEMENTARY_INFORMATION.md",
        "analysis/ch1/submission_contract.json",
    )
    for relative in required:
        assert (ROOT / relative).is_file(), relative


def test_withdrawn_lability_not_active_release_stage() -> None:
    payload = json.loads((ROOT / "analysis/ch1/pipeline.json").read_text(encoding="utf-8"))
    assert payload["stages"]["legacy_lability_withdrawn"]["status"].startswith("retained for provenance only")
    assert "submission_validation" not in payload["stages"]
    assert "submission_manifest" not in payload["stages"]
    assert "eazami_trait_handoff" not in payload["stages"]
