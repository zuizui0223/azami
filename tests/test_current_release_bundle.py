from pathlib import Path
import hashlib
import json

import pytest

from reproducibility.build_current_release_bundle import (
    NATIVE_SHA,
    TAXONOMY_ARCHIVE_SHA256,
    TAXONOMY_ARTIFACT_ID,
    TAXONOMY_NATIVE_MEMBER,
    REFERENCE,
    REFERENCE_MANIFEST,
    input_contract,
    locate_archive,
    release_gaps,
    is_manuscript_document_path,
    tracked_manuscript_documents,
    verify_manifest_files,
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
    assert TAXONOMY_ARTIFACT_ID == 10292140117
    assert TAXONOMY_ARCHIVE_SHA256 == "dfb6eec3001e3a984662d5aba06cda5fa80e144b36ccb4af9fdf45973854edc5"
    assert TAXONOMY_NATIVE_MEMBER == "input/observation_native_status.csv"


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
        10382954095,
        10382578412,
        10903835882,
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
    replay = staged[10382954095]
    assert replay["drive_file_id"] == "1URtYSHHM-yl3wrTAtTJRkk2FDPiT9PVl"
    assert replay["archive_sha256"] == "2064b63366681aa3bd708fef00d2e7ac50d20a7fad771849b793c0d2b6957a0e"
    assert replay["completed_numerical_stages"] == 8
    assert replay["aggregate_files_compared"] == 16
    assert replay["validation_status"] == "PASS"
    figure = staged[10382578412]
    assert figure["drive_file_id"] == "19-xLleKgq2Gj3MjQtidXsb4A6E-yih_4"
    assert figure["png_sha256"] == "a4c103fc23d1c2640ca87601bcf8d05e59826f975d943c72a24cf0c61d31452c"
    full_figure_surface = staged[10903835882]
    assert full_figure_surface["drive_file_id"] == "1AnrI_UttuQYRzTuhYRnn5Y74RT3Uu5u0"
    assert full_figure_surface["checksum_verification"] == "PASS"
    assert full_figure_surface["main_figure_count"] == 5
    assert full_figure_surface["supporting_figure_count"] == 8
    assert full_figure_surface["document_pagination_validated"] is False
    zenodo = receipt["zenodo_staging_bundle"]
    assert zenodo["status"] == "assembled_and_verified_not_public"
    assert zenodo["source_main_commit"] == "b128d8e3cb4c326c1c038f5d1e7f22ffe37b48b4"
    assert zenodo["github_actions_artifact_id"] == 10919753688
    assert zenodo["github_actions_artifact_sha256"] == "a16053d5dbe50569d08494b594a4e92a6ad011ce98c588f0e06082cddfb1ad3a"
    assert zenodo["inner_bundle_sha256"] == "dee0a37be6b38f7fa72a48053021cb1709031b93ee1bf18043b75892ca55f051"
    assert zenodo["durable_drive_file_id"] == "1H4MXMtmSl26lb8HNrc0t-PrdcT8vxBMH"
    assert zenodo["reference_files"] == 16
    assert zenodo["figure_manifest_files"] == 30
    assert zenodo["release_ready"] is False
    assert zenodo["manuscript_files_included"] is False
    assert zenodo["document_pagination_is_submission_only"] is True
    assert set(zenodo["release_gaps"]) == {"release_metadata"}
    assert zenodo["includes"]["manuscript_docx_pdf"] is False
    policy = receipt["zenodo_release_policy"]
    assert policy["manuscript_upload"] is False
    assert policy["manuscript_files_included"] is False
    assert policy["document_pagination_is_submission_only"] is True
    assert set(policy["current_expected_release_gaps"]) == {"release_metadata"}
    assert policy["code_snapshot_manuscript_guard"] == "fail_closed"
    docqa = receipt["document_figure_qa"]
    assert docqa["supporting_information"]["pages_inspected"] == 13
    assert docqa["supporting_information"]["layout_status"] == "PASS"
    assert docqa["main_manuscript"]["canonical_final_document_found"] is False
    assert docqa["main_manuscript"]["status"] == "OPEN"
    assert docqa["overall_document_pagination_validated"] is False
    assert docqa["zenodo_release_gate"] is False
    assert docqa["scope"] == "journal_submission_only"
    metaprep = receipt["release_metadata_preparation"]
    assert metaprep["status"] == "prepared_not_approved"
    assert metaprep["non_author_fields_prepared"] is True
    assert metaprep["release_approved"] is False
    assert receipt["reference_file_count"] == 16
    assert receipt["scientific_outputs_changed"] is False
    assert receipt["public_release_changed"] is False

def test_figure_manifest_accepts_portable_manifest_relative_files(tmp_path):
    figures = tmp_path / "figures"
    figures.mkdir()
    image = figures / "Figure_1.png"
    image.write_bytes(b"figure-bytes")
    digest = hashlib.sha256(image.read_bytes()).hexdigest()
    manifest = tmp_path / "final_figure_manifest.json"
    manifest.write_text(json.dumps({
        "files": [{"path": "figures/Figure_1.png", "sha256": digest}]
    }), encoding="utf-8")
    rows = verify_manifest_files(manifest)
    assert len(rows) == 1
    assert rows[0]["source_kind"] == "manifest_relative"
    assert rows[0]["source_path"] == image.resolve()
    assert rows[0]["path"] == "figures/Figure_1.png"

def test_document_pagination_is_not_a_zenodo_release_gate(tmp_path):
    manifest = tmp_path / "final_figure_manifest.json"
    manifest.write_text(json.dumps({"document_pagination_validated": False, "files": [
        {"path": "placeholder.png", "sha256": "0" * 64}
    ]}), encoding="utf-8")
    assert release_gaps(manifest, tmp_path / "release-metadata.json") == []


def test_manuscript_documents_are_forbidden_from_zenodo_code_snapshot():
    assert is_manuscript_document_path("submission/Main_manuscript.docx")
    assert is_manuscript_document_path("submission/Supporting_Information.pdf")
    assert is_manuscript_document_path("submission/title_page.doc")
    assert is_manuscript_document_path("submission/cover-letter.pdf")
    assert not is_manuscript_document_path("reproducibility/Figure_1.pdf")
    assert not is_manuscript_document_path("analysis/v3/README.md")
    assert tracked_manuscript_documents() == []

def test_prepared_release_metadata_fills_non_author_fields_but_fails_closed():
    path = ROOT / "reproducibility/zenodo_release_metadata.prepared.json"
    obj = json.loads(path.read_text(encoding="utf-8"))
    assert obj["schema_version"] == 1
    assert obj["release_approved"] is False
    assert obj["scientific_scope"] == SCOPE_ID
    assert obj["title"] == "Azami Chapter 1 current numerical and figure release"
    assert obj["archive_strategy"] == "new_version_existing_concept"
    assert obj["preserved_v2_record"] == {
        "record_doi": V2_RECORD_DOI,
        "concept_doi": V2_CONCEPT_DOI,
        "modify_existing_v2": False,
    }
    assert obj["licensing"]["software"] == "MIT"
    assert "third-party" in obj["licensing"]["data_and_third_party_strategy"].lower()
    assert "TBD" not in path.read_text(encoding="utf-8")
    assert obj["creators"] == []
    with pytest.raises(ValueError, match="release_approved"):
        validate_release_metadata(path, expected_head="not-final")


def test_document_qa_receipt_is_submission_only():
    text = (ROOT / "reproducibility/DOCUMENT_FIGURE_QA_20260927.md").read_text(encoding="utf-8")
    assert "### Layout result: PASS" in text
    assert "13 rendered pages" in text
    assert "Main document QA remains **OPEN**" in text
    assert "journal-submission document QA remains open" in text
    assert "Zenodo release readiness is evaluated independently" in text
    assert "must not be added to the Zenodo bundle" in text

