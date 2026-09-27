from pathlib import Path
import json

import pytest

from reproducibility.build_current_release_bundle import (
    NATIVE_SHA,
    TAXONOMY_ARCHIVE_SHA256,
    TAXONOMY_ARTIFACT_ID,
    TAXONOMY_NATIVE_MEMBER,
    input_contract,
    locate_archive,
    release_gaps,
)
from reproducibility.release_metadata_contract import (
    SCOPE_ID,
    V2_CONCEPT_DOI,
    V2_RECORD_DOI,
    validate_release_metadata,
)
from reproducibility.run_current_analysis import INPUTS, commands

ROOT = Path(__file__).resolve().parents[1]


def test_release_input_contract_matches_current_runner():
    contract = input_contract()
    assert set(contract) == set(INPUTS) == {"continuous", "environment", "spatial", "historical"}
    for role, (artifact_id, archive_sha, members) in INPUTS.items():
        row = contract[role]
        assert row["artifact_id"] == artifact_id
        assert row["archive_sha256"] == archive_sha
        assert set(row["members"]) == set(members)
        for source, (target, member_sha) in members.items():
            assert row["members"][source] == {
                "analysis_path": target,
                "sha256": member_sha,
            }

    assert NATIVE_SHA == "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a"
    assert TAXONOMY_ARTIFACT_ID == 10292140117
    assert TAXONOMY_ARCHIVE_SHA256 == "dfb6eec3001e3a984662d5aba06cda5fa80e144b36ccb4af9fdf45973854edc5"
    assert TAXONOMY_NATIVE_MEMBER == "input/observation_native_status.csv"


def test_current_replay_still_has_eight_stages(tmp_path):
    cmds = commands(tmp_path / "inputs", tmp_path / "results")
    assert len(cmds) == 8
    assert cmds[-1][0] == "analysis.v3.run_rv_estimator_validity"
    assert "--equal-n-replicates" in cmds[-1]
    assert cmds[-1][cmds[-1].index("--equal-n-replicates") + 1] == "1000"


def test_zenodo_release_gate_is_metadata_only():
    assert release_gaps(None) == ["release_metadata"]
    assert release_gaps(Path("release-metadata.json")) == []


def test_archive_discovery_requires_one_exact_artifact_match(tmp_path):
    artifact_id = 9633419268
    with pytest.raises(FileNotFoundError, match="found 0"):
        locate_archive(tmp_path, artifact_id)

    one = tmp_path / f"artifact-{artifact_id}-environment.zip"
    one.write_bytes(b"one")
    assert locate_archive(tmp_path, artifact_id) == one

    (tmp_path / f"backup-{artifact_id}.zip").write_bytes(b"two")
    with pytest.raises(FileNotFoundError, match="found 2"):
        locate_archive(tmp_path, artifact_id)


def test_prepared_metadata_is_data_only_and_inherits_v2_dataset_creator():
    path = ROOT / "reproducibility/zenodo_release_metadata.prepared.json"
    obj = json.loads(path.read_text(encoding="utf-8"))
    assert obj["schema_version"] == 1
    assert obj["release_approved"] is False
    assert obj["scientific_scope"] == SCOPE_ID
    assert obj["title"] == "Azami Chapter 1 v3 numerical analysis input package"
    assert obj["archive_strategy"] == "new_version_existing_concept"
    assert obj["creators"] == [
        {"name": "ZHANG, Ruiqi", "orcid": None, "affiliation": None}
    ]
    assert obj["preserved_v2_record"] == {
        "record_doi": V2_RECORD_DOI,
        "concept_doi": V2_CONCEPT_DOI,
        "modify_existing_v2": False,
    }
    description = obj["description"].lower()
    assert "analysis code" in description
    assert "intentionally excluded" in description
    assert "figures" in description
    assert "manuscript" in description
    assert "github code only" in obj["licensing"]["software"].lower()
    assert "TBD" not in path.read_text(encoding="utf-8")
    with pytest.raises(ValueError, match="release_approved"):
        validate_release_metadata(path, expected_head="not-final")


def test_metadata_validator_accepts_approved_data_only_contract(tmp_path):
    head = "0123456789abcdef"
    metadata = {
        "schema_version": 1,
        "release_approved": True,
        "scientific_scope": SCOPE_ID,
        "final_code_commit": head,
        "title": "Azami Chapter 1 v3 numerical analysis input package",
        "description": "Exact processed analysis inputs paired to GitHub code.",
        "creators": [{"name": "ZHANG, Ruiqi", "orcid": None, "affiliation": None}],
        "archive_strategy": "new_version_existing_concept",
        "preserved_v2_record": {
            "record_doi": V2_RECORD_DOI,
            "concept_doi": V2_CONCEPT_DOI,
            "modify_existing_v2": False,
        },
        "licensing": {
            "software": "MIT (GitHub code only; not included)",
            "data_and_third_party_strategy": "Preserve source-specific terms.",
        },
        "anonymous_redownload_required": True,
        "clean_replay_required": True,
    }
    path = tmp_path / "release.json"
    path.write_text(json.dumps(metadata), encoding="utf-8")
    assert validate_release_metadata(path, expected_head=head) == metadata


def test_staging_receipt_declares_data_only_zenodo_policy():
    receipt = json.loads(
        (ROOT / "reproducibility/CURRENT_RELEASE_STAGING_20260912.json").read_text()
    )
    policy = receipt["zenodo_release_policy"]
    assert policy["archive_role"] == "durable_processed_analysis_inputs"
    assert policy["analysis_code_location"] == "github"
    assert policy["code_included"] is False
    assert policy["manuscript_files_included"] is False
    assert policy["figures_included"] is False
    assert policy["reference_outputs_included"] is False
    assert policy["replay_receipts_included"] is False
    assert policy["analysis_input_count"] == 5
    assert set(policy["current_expected_release_gaps"]) == {"release_metadata"}

    bundle = receipt["zenodo_staging_bundle"]
    assert bundle["status"] == "assembled_and_verified_not_public"
    assert bundle["archive_role"] == "durable_processed_analysis_inputs"
    assert bundle["source_main_commit"] == "810637d9dfdde9c4b810c2472f43f254786ce614"
    assert bundle["workflow_run"] == 36287521925
    assert bundle["github_actions_artifact_id"] == 10920813215
    assert bundle["github_actions_artifact_sha256"] == "0afaed34da6654819bf979e676f50f877a816beac1d281fe92b6735eee5f6faa"
    assert bundle["inner_bundle_sha256"] == "c4ec876206e6f8ef0cd69d126fa31b2b71aa1b1919172aeb371f75bbabcd6d10"
    assert bundle["durable_drive_file_id"] == "1Bqj_5pLAd7DQmC5OLs9x426i8XwHsgZh"
    assert bundle["analysis_input_count"] == 5
    assert bundle["code_included"] is False
    assert bundle["manuscript_files_included"] is False
    assert bundle["figures_included"] is False
    assert bundle["reference_outputs_included"] is False
    assert bundle["replay_receipts_included"] is False
    assert set(bundle["release_gaps"]) == {"release_metadata"}

def test_data_release_approval_is_explicit_and_data_only():
    path = ROOT / "reproducibility/zenodo_data_release_approval.json"
    obj = json.loads(path.read_text(encoding="utf-8"))
    assert obj["schema_version"] == 1
    assert obj["approved"] is True
    assert obj["scientific_scope"] == SCOPE_ID
    assert obj["archive_role"] == "durable_processed_analysis_inputs"
    assert obj["creator"] == "ZHANG, Ruiqi"
    assert obj["archive_strategy"] == "new_version_existing_concept"
    assert obj["preserved_concept_doi"] == V2_CONCEPT_DOI
    boundary = obj["public_payload_boundary"]
    assert boundary["analysis_input_count"] == 5
    for key in (
        "code_included",
        "manuscript_files_included",
        "figures_included",
        "reference_outputs_included",
        "replay_receipts_included",
        "original_photographs_included",
        "detector_training_included",
    ):
        assert boundary[key] is False
    assert "CC BY 4.0" in obj["licensing_approval"]
    assert obj["public_upload_performed"] is False

def test_release_ready_data_only_candidate_is_frozen_in_receipt():
    receipt = json.loads(
        (ROOT / "reproducibility/CURRENT_RELEASE_STAGING_20260912.json").read_text()
    )
    final = receipt["zenodo_final_candidate"]
    assert final["status"] == "release_ready_not_public"
    assert final["source_main_commit"] == "7483ae8e20fa5df2439ad524312d3ac9a0fbf71f"
    assert final["workflow_run"] == 36288310738
    assert final["github_actions_artifact_id"] == 10921221306
    assert final["github_actions_wrapper_sha256"] == "941c70768664ae9df82cd18f522dda9f4fd352ad723e08ad1e9a8d9e23f0996e"
    assert final["inner_bundle_sha256"] == "be57e9ba80773c97afdc200f9a3982eaa85044cfc6791b62187ccfcd41f54438"
    assert final["durable_drive_file_id"] == "16pildr5nOL08HAb6ATEHepUPei5oEn0t"
    assert final["release_ready"] is True
    assert final["release_gaps"] == []
    assert final["analysis_input_count"] == 5
    assert final["code_included"] is False
    assert final["manuscript_files_included"] is False
    assert final["figures_included"] is False
    assert final["reference_outputs_included"] is False
    assert final["replay_receipts_included"] is False
    assert final["checksum_verification"] == "PASS"
    assert final["public_upload_performed"] is False
    assert len(final["exact_members"]) == 10
    assert receipt["release_metadata_preparation"]["release_approved"] is True
    assert receipt["release_metadata_preparation"]["remaining_author_owned_fields"] == []

def test_zenodo_publication_handoff_matches_final_candidate():
    text = (ROOT / "reproducibility/ZENODO_PUBLICATION_HANDOFF_20260927.md").read_text(
        encoding="utf-8"
    )
    assert "azami_ch1_v3_analysis_inputs.zip" in text
    assert "be57e9ba80773c97afdc200f9a3982eaa85044cfc6791b62187ccfcd41f54438" in text
    assert "7483ae8e20fa5df2439ad524312d3ac9a0fbf71f" in text
    assert "10.5281/zenodo.22295790" in text
    assert "code_included=false" in text
    assert "manuscript_files_included=false" in text
    assert "credential-free public reproduction verified" in text

