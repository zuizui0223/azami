from pathlib import Path
import hashlib
import json

import pytest

from reproducibility.build_current_release_bundle import (
    NATIVE_SHA,
    REFERENCE,
    REFERENCE_MANIFEST,
    input_contract,
    locate_archive,
    release_gaps,
    verified_reference_rows,
)
from reproducibility.release_metadata_contract import (
    SCOPE_ID,
    V2_CONCEPT_DOI,
    V2_RECORD_DOI,
    validate_release_metadata,
)
from reproducibility.run_current_analysis import INPUTS

ROOT = Path(__file__).resolve().parents[1]


def test_release_input_contract_is_exactly_the_current_runner_contract():
    contract = input_contract()
    assert set(contract) == set(INPUTS) == {"continuous", "environment", "spatial", "historical"}
    for role, (artifact_id, archive_sha, members) in INPUTS.items():
        row = contract[role]
        assert row["artifact_id"] == artifact_id
        assert row["archive_sha256"] == archive_sha
        assert set(row["members"]) == set(members)
        for source, (target, member_sha) in members.items():
            assert row["members"][source] == {"bundle_path": target, "sha256": member_sha}
    assert NATIVE_SHA == "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a"


def test_current_reference_release_surface_is_frozen_and_verified():
    rows = verified_reference_rows()
    assert len(rows) == 15
    manifest = json.loads(REFERENCE_MANIFEST.read_text(encoding="utf-8"))
    assert manifest["files"] == rows
    for row in rows:
        path = REFERENCE / row["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]


def test_final_release_fails_closed_on_unfrozen_surfaces():
    assert release_gaps(None, None) == ["final_figure_manifest", "release_metadata"]
    dummy = Path("figure-manifest.json")
    assert release_gaps(dummy, None) == ["release_metadata"]
    assert release_gaps(dummy, Path("release-metadata.json")) == []


def test_release_metadata_template_is_deliberately_not_final():
    template = ROOT / "reproducibility/zenodo_release_metadata.template.json"
    with pytest.raises(ValueError, match="release_approved"):
        validate_release_metadata(template, expected_head="abc123")


def test_release_metadata_contract_accepts_only_approved_head_bound_metadata(tmp_path):
    head = "0123456789abcdef"
    metadata = {
        "schema_version": 1,
        "release_approved": True,
        "scientific_scope": SCOPE_ID,
        "final_code_commit": head,
        "title": "Azami Chapter 1 current numerical and figure release",
        "description": "Frozen inputs, current aggregate references, code, figure provenance and replay receipts for the submitted Chapter 1 analysis.",
        "creators": [{"name": "Example Author"}],
        "archive_strategy": "new_version_existing_concept",
        "preserved_v2_record": {
            "record_doi": V2_RECORD_DOI,
            "concept_doi": V2_CONCEPT_DOI,
            "modify_existing_v2": False,
        },
        "licensing": {
            "software": "MIT",
            "data_and_third_party_strategy": "Preserve source-specific terms and separate software licensing from data and third-party material.",
        },
        "anonymous_redownload_required": True,
        "clean_replay_required": True,
    }
    path = tmp_path / "release.json"
    path.write_text(json.dumps(metadata), encoding="utf-8")
    assert validate_release_metadata(path, expected_head=head) == metadata

    bad = dict(metadata)
    bad["final_code_commit"] = "wrong"
    path.write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="does not match HEAD"):
        validate_release_metadata(path, expected_head=head)

    bad = dict(metadata)
    bad["title"] = "TBD"
    path.write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="title"):
        validate_release_metadata(path, expected_head=head)


def test_archive_discovery_requires_one_exact_artifact_match(tmp_path):
    artifact_id = 9633419268
    try:
        locate_archive(tmp_path, artifact_id)
    except FileNotFoundError as exc:
        assert "found 0" in str(exc)
    else:
        raise AssertionError("missing archive must fail closed")

    one = tmp_path / f"artifact-{artifact_id}-environment.zip"
    one.write_bytes(b"one")
    assert locate_archive(tmp_path, artifact_id) == one

    (tmp_path / f"backup-{artifact_id}.zip").write_bytes(b"two")
    try:
        locate_archive(tmp_path, artifact_id)
    except FileNotFoundError as exc:
        assert "found 2" in str(exc)
    else:
        raise AssertionError("ambiguous archive selection must fail closed")


def test_staging_receipt_tracks_all_recovered_current_only_artifacts():
    receipt = json.loads((ROOT / "reproducibility/CURRENT_RELEASE_STAGING_20260912.json").read_text())
    staged = {row["github_actions_artifact_id"]: row for row in receipt["newly_staged_artifacts"]}
    assert set(staged) == {
        9633419268,
        10130210432,
        10131007603,
        10136229131,
        10135679053,
        10291656193,
    }
    assert staged[10291656193]["drive_file_id"] == "1UnPG1mhjdJKbXFLJDbn6l-TrFtAXgJSd"
    assert staged[10291656193]["archive_sha256"] == "2fed9448c2210af4ded7a4ccc5cbe6f543b64e8f19ba870f5200fa095290e766"
    assert receipt["reference_file_count"] == 15
    assert receipt["scientific_outputs_changed"] is False
    assert receipt["public_release_changed"] is False
