from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AVAILABILITY = ROOT / "reproducibility" / "material_availability.json"
DURABLE = ROOT / "reproducibility" / "durable_archive_manifest.json"
CHELSA = ROOT / "analysis" / "ch1" / "chelsa_process_environment_sources.json"
RUNBOOK = ROOT / "reproducibility" / "README.md"
VERIFIER = ROOT / "reproducibility" / "verify_local_materials.py"


def test_all_github_materials_marked_present_really_exist() -> None:
    payload = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    for item in payload["github_repository_materials"]:
        assert item["status"] == "present"
        paths = item.get("paths", [item.get("path")])
        for relative in paths:
            assert relative is not None
            assert (ROOT / relative).exists(), relative


def test_durable_numerical_input_sha_matches_drive_manifest() -> None:
    availability = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    durable = json.loads(DURABLE.read_text(encoding="utf-8"))
    by_id = {int(row["artifact_id"]): row for row in durable["direct_files"]}
    for row in availability["durable_owner_archive_materials"]["numerical_inputs"]:
        artifact_id = int(row["artifact_id"])
        assert artifact_id in by_id
        assert row["sha256"] == by_id[artifact_id]["sha256"]
    assert availability["durable_owner_archive_materials"]["status"] == "complete_and_readback_verified"


def test_chelsa_external_boundary_matches_frozen_source_registry() -> None:
    availability = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    registry = json.loads(CHELSA.read_text(encoding="utf-8"))
    external = availability["public_external_source_materials"]
    assert external["status"] == "not_vendored_in_git_or_owner_archive_at_audit_time"
    assert [row["column"] for row in registry["sources"]] == external["layers"]
    assert len(registry["sources"]) == 5
    assert all(row["url"].startswith("https://") for row in registry["sources"])
    conclusion = availability["material_audit_conclusion"]
    assert conclusion["raw_process_raster_reconstruction_fully_offline"] is False
    assert "five frozen CHELSA process rasters" in conclusion["reason"]


def test_windows_runbook_and_local_verifier_are_tracked() -> None:
    assert VERIFIER.is_file()
    text = RUNBOOK.read_text(encoding="utf-8")
    assert r"C:\Users\zuizui\OneDrive - Kyoto University\デスクトップ\azami論文材料" in text
    assert "verify_local_materials.py" in text
    assert "CHELSA process rasters" in text
    assert "PASS 24/24" in text
    assert "PASS 7/7" in text
    assert "PASS 11/11" in text
