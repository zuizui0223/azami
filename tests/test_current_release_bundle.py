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
from reproducibility.run_current_analysis import INPUTS, commands

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


def test_current_replay_has_eight_stages_and_ends_with_estimator_validity(tmp_path):
    cmds = commands(tmp_path / "inputs", tmp_path / "results")
    assert len(cmds) == 8
    assert cmds[-1][0] == "analysis.v3.run_rv_estimator_validity"
    assert "--equal-n-replicates" in cmds[-1]
    assert cmds[-1][cmds[-1].index("--equal-n-replicates") + 1] == "1000"


def test_current_reference_release_surface_is_frozen_and_verified():
    rows = verified_reference_rows()
    assert len(rows) == 16
    manifest = json.loads(REFERENCE_MANIFEST.read_text(encoding="utf-8"))
    assert manifest["files"] == rows
    for row in rows:
        path = REFERENCE / row["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
    estimator = next(row for row in rows if row["path"] == "estimator_validity/rv_estimator_validity_summary.json")
    assert estimator["artifact"] == 10382387052
    assert estimator["archive_sha256"] == "345866a7e333f78677ad3797811e2cbe82d3e09ee061f7642df7f4c4d5ec008e"


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
        10292140117,
        10382387052,
        10382578412,
    }
    assert staged[10291656193]["drive_file_id"] == "1UnPG1mhjdJKbXFLJDbn6l-TrFtAXgJSd"
    assert staged[10291656193]["archive_sha256"] == "2fed9448c2210af4ded7a4ccc5cbe6f543b64e8f19ba870f5200fa095290e766"
    taxonomy = staged[10292140117]
    assert taxonomy["drive_file_id"] == "1FYSiKtJDBq9m4sWxlMuRBvLirpHj90Ja"
    assert taxonomy["archive_sha256"] == "dfb6eec3001e3a984662d5aba06cda5fa80e144b36ccb4af9fdf45973854edc5"
    assert taxonomy["headline_taxonomy_robust"] is True
    assert taxonomy["contains_exact_native_status"] is True
    assert taxonomy["native_status_sha256_after_permitted_newline_normalization"] == NATIVE_SHA
    estimator = staged[10382387052]
    assert estimator["drive_file_id"] == "1JvjPWBG-EuouuTVgLjm6VYyjk7ui58Sd"
    assert estimator["equal_n_strength_gate_pass"] is True
    figure = staged[10382578412]
    assert figure["drive_file_id"] == "19-xLleKgq2Gj3MjQtidXsb4A6E-yih_4"
    assert figure["png_sha256"] == "a4c103fc23d1c2640ca87601bcf8d05e59826f975d943c72a24cf0c61d31452c"
    assert receipt["reference_file_count"] == 16
    assert receipt["scientific_outputs_changed"] is False
    assert receipt["public_release_changed"] is False
