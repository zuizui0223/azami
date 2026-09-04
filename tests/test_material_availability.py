from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AVAILABILITY = ROOT / "reproducibility" / "material_availability.json"
DURABLE = ROOT / "reproducibility" / "durable_archive_manifest.json"
PUBLIC_RELEASE = ROOT / "reproducibility" / "public_release_manifest.json"
CHELSA = ROOT / "analysis" / "ch1" / "chelsa_process_environment_sources.json"
RUNBOOK = ROOT / "reproducibility" / "README.md"
VERIFIER = ROOT / "reproducibility" / "verify_local_materials.py"
ZENODO_METADATA = ROOT / "reproducibility" / "zenodo_metadata.json"


def test_all_github_materials_marked_present_really_exist() -> None:
    payload = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    for item in payload["github_repository_materials"]:
        assert item["status"] == "present"
        paths = item.get("paths", [item.get("path")])
        for relative in paths:
            assert relative is not None
            assert (ROOT / relative).exists(), relative


def test_owner_archive_preservation_matches_public_release_hash_contract() -> None:
    availability = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    durable = json.loads(DURABLE.read_text(encoding="utf-8"))
    public = json.loads(PUBLIC_RELEASE.read_text(encoding="utf-8"))
    by_id = {int(row["artifact_id"]): row for row in durable["direct_files"]}
    public_by_name = {row["filename"]: row for row in public["minimum_analysis_inputs"]}
    for row in availability["durable_owner_archive_materials"]["numerical_inputs"]:
        artifact_id = int(row["artifact_id"])
        assert artifact_id in by_id
        assert row["sha256"] == by_id[artifact_id]["sha256"]
        assert row["filename"] in public_by_name
        assert row["sha256"] == public_by_name[row["filename"]]["sha256"]
    assert availability["durable_owner_archive_materials"]["status"] == (
        "complete_and_readback_verified_but_not_public_distribution"
    )


def test_public_release_gate_records_published_doi_pending_anonymous_verification() -> None:
    availability = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    public = json.loads(PUBLIC_RELEASE.read_text(encoding="utf-8"))
    third_party = availability["third_party_reproducibility"]
    assert third_party["status"] == "published_pending_anonymous_redownload_verification"
    assert third_party["full_numerical_reproduction_ready"] is False
    assert third_party["required_analysis_inputs_publicly_downloadable_without_owner_credentials"] is None

    assert public["status"] == "published_pending_anonymous_redownload_verification"
    release = public["public_data_release"]
    assert release["published"] is True
    assert release["doi"] == "10.5281/zenodo.22295791"
    assert release["url"] == "https://zenodo.org/records/22295791"
    assert release["version"] == "v1"
    assert release["anonymous_redownload_verified"] is False
    assert release["publicly_downloadable_without_owner_credentials"] is None
    assert public["code"]["code_ref"] == "584af97b050d15701f26ce1facea212d5b648d4d"
    assert len(public["minimum_analysis_inputs"]) == 4

    bundle = public["staging_bundle"]
    assert bundle["filename"] == "azami_ch1_v2_reproduction_inputs_2026-09-04.zip"
    assert bundle["sha256"] == "50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677"
    assert bundle["size_bytes"] == 56942044
    assert bundle["contains_exact_minimum_analysis_inputs"] is True
    assert bundle["metadata_file"] == "reproducibility/zenodo_metadata.json"
    assert ZENODO_METADATA.is_file()


def test_chelsa_external_boundary_matches_frozen_source_registry() -> None:
    availability = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    public = json.loads(PUBLIC_RELEASE.read_text(encoding="utf-8"))
    registry = json.loads(CHELSA.read_text(encoding="utf-8"))
    external = availability["public_external_source_materials"]
    assert external["status"] == "public_urls_frozen_but_raster_bytes_not_archived_in_public_release"
    assert [row["column"] for row in registry["sources"]] == external["layers"]
    assert [row["column"] for row in registry["sources"]] == public["public_environment_sources"]["layers"]
    assert len(registry["sources"]) == 5
    assert all(row["url"].startswith("https://") for row in registry["sources"])
    assert public["public_environment_sources"]["reconstructed_environment_sha256"] == (
        "e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7"
    )
    conclusion = availability["material_audit_conclusion"]
    assert conclusion["raw_process_raster_reconstruction_fully_offline"] is False
    assert conclusion["third_party_full_numerical_reproduction_ready"] is False


def test_third_party_runbook_has_no_author_local_dependency() -> None:
    assert VERIFIER.is_file()
    assert PUBLIC_RELEASE.is_file()
    text = RUNBOOK.read_text(encoding="utf-8")
    assert "independent reader, reviewer, or researcher" in text
    assert "10.5281/zenodo.22295791" in text
    assert "verify_local_materials.py" in text
    assert "PASS 24/24" in text
    assert "PASS 7/7" in text
    assert "PASS 11/11" in text
    assert "C:\\Users\\zuizui" not in text
    assert "OneDrive - Kyoto University" not in text
    assert "owner-controlled Google Drive" not in text
