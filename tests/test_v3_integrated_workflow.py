import json
import shutil
from pathlib import Path

import pytest

from analysis.v3.integrated_preflight import CONTRACT, INDEX, _check_schedule_replay, _check_source_recovery, _check_protected_replay, validate

ROOT = Path(__file__).resolve().parents[1]


def clone_inputs(tmp_path):
    index = json.loads((ROOT / INDEX).read_text(encoding="utf-8"))
    paths = [CONTRACT, INDEX, "analysis/v3/ecological_analysis_contract.json"]
    paths += ['reproducibility/v3_full_native_environment_qc_selection_20260908.json',
              'analysis/v3/environment_production_contract.json']
    paths += [item["path"] for item in index["inputs"]]
    paths += [item["implementation"] for item in index["requirements"].values() if item.get("implementation")]
    for relative in paths:
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, destination)
    return tmp_path


def test_actual_integrated_preflight_has_no_execution_promotion():
    result = validate(ROOT)
    assert result["status"] == "INTEGRATED_DESIGN_INTEGRITY_VERIFIED_EXECUTION_INCOMPLETE"
    assert [s["stage"] for s in result["stages"]] == ["source", "measurement", "assessability", "ecology", "synthesis"]
    assert [s["optional_for_primary_ecology"] for s in result["stages"]] == [False] * 4 + [True]
    assert result["ecological_fitting_authorized"] is False
    assert result["full_original_stream_authorized"] is False
    assert len(result["verified_public_evidence"]) == 7


def test_pinned_evidence_change_is_rejected_even_if_status_unchanged(tmp_path):
    root = clone_inputs(tmp_path)
    path = root / "analysis/v3/measurement_qualification_decision_20260908.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    value["endpoints"][0]["followup_common_usable"] = 9999
    path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ValueError, match="Pinned evidence changed"):
        validate(root)


def test_status_is_checked_against_receipt_not_repeated_test_literal(tmp_path):
    root = clone_inputs(tmp_path)
    path = root / "analysis/v3/ecological_analysis_contract.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    value["completed_design_gates"]["measurement_freeze_status"] = "FROZEN_BEFORE_ECOLOGY"
    path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ValueError, match="status differs"):
        validate(root)


@pytest.mark.parametrize("field,value", [("status", "production_authorized"), ("stages", [])])
def test_mutated_authority_cannot_enable_execution(tmp_path, field, value):
    root = clone_inputs(tmp_path)
    path = root / CONTRACT
    contract = json.loads(path.read_text(encoding="utf-8"))
    contract[field] = value
    path.write_text(json.dumps(contract), encoding="utf-8")
    with pytest.raises(ValueError):
        validate(root)


def test_secondary_synthesis_cannot_become_primary_prerequisite(tmp_path):
    root = clone_inputs(tmp_path)
    path = root / CONTRACT
    contract = json.loads(path.read_text(encoding="utf-8"))
    contract["hypervolume"]["optional"] = False
    path.write_text(json.dumps(contract), encoding="utf-8")
    with pytest.raises(ValueError, match="Secondary hypervolume"):
        validate(root)


def test_primary_formulation_stays_inside_verified_source_selected_exposures(tmp_path):
    root = clone_inputs(tmp_path)
    path = root / CONTRACT
    contract = json.loads(path.read_text(encoding="utf-8"))
    contract["environment_and_inference"]["retained_variables"].append("tasmax_month")
    path.write_text(json.dumps(contract), encoding="utf-8")
    with pytest.raises(ValueError, match="formulation differs"):
        validate(root)


def test_public_receipt_pins_are_newline_independent(tmp_path):
    root = clone_inputs(tmp_path)
    index = json.loads((root / INDEX).read_text(encoding="utf-8"))
    for item in index["inputs"]:
        path = root / item["path"]
        data = json.loads(path.read_text(encoding="utf-8"))
        path.write_text(json.dumps(data, ensure_ascii=True, indent=4), encoding="utf-8", newline="\r\n")
    assert validate(root)["ecological_fitting_authorized"] is False


def test_two_contributions_are_questions_not_required_positive_results():
    contract = json.loads((ROOT / CONTRACT).read_text(encoding="utf-8"))
    assert contract["source_and_sampling"]["master_deletion"] is False
    assert contract["source_and_sampling"]["primary_range"] == "native_only"
    assert "not 14 independently calibrated" in contract["stages"][1]["output"]
    assert "not claimed as new inventions" in contract["two_contributions"]["methodological"]
    assert "without requiring either v2 candidate to survive" in contract["two_contributions"]["ecological"]
    assert "36 tests" in contract["environment_and_inference"]["primary_probability_family"]
    assert "not evidence of a difference" in contract["environment_and_inference"]["scale_contrast"]


@pytest.mark.parametrize("mode", ["missing", "duplicate_id", "duplicate_path", "unbound_requirement"])
def test_every_historical_evidence_claim_has_a_unique_verified_input(tmp_path, mode):
    root = clone_inputs(tmp_path)
    path = root / INDEX
    index = json.loads(path.read_text(encoding="utf-8"))
    if mode == "missing":
        index["inputs"] = index["inputs"][1:]
    elif mode == "duplicate_id":
        index["inputs"].append(index["inputs"][0])
    elif mode == "duplicate_path":
        index["inputs"][0]["path"] = index["inputs"][1]["path"]
    else:
        index["requirements"]["reconciled_source_links"]["evidence"] = "not_verified.json"
    path.write_text(json.dumps(index), encoding="utf-8")
    with pytest.raises(ValueError, match="Evidence index|not bound"):
        validate(root)


def test_cloud_source_rebuild_stops_before_repeating_lossy_execution():
    workflow = (ROOT / ".github/workflows/ch1-v3-ecological-source-cohort.yml").read_text(encoding="utf-8")
    assert "  push:" not in workflow
    assert workflow.index("exit 1") < workflow.index("actions/checkout")
    assert "exact source, authority, native join and cohort persistence/recovery" in workflow
    assert workflow.count("if: always() && steps.private_archive.outcome == 'success'") == 2


@pytest.mark.parametrize("section,field,value", [
    ("cohort", "historical_csv_sha256", "0" * 64),
    ("cohort", "rows", 319243),
    ("authority", "historical_http_response_identity_verified", True),
    ("authority", "all_three_outputs_byte_identical_in_offline_replay", False),
    ("enrichment", "cohort_membership_changed", True),
    ("private_preservation", "cloud_execution_hold_preserved", False),
])
def test_source_recovery_receipt_cannot_change_identity_or_promote_claims(section, field, value):
    receipt = json.loads((ROOT / "reproducibility/v3_source_cohort_recovery_20260908.json").read_text())
    source = json.loads((ROOT / "reproducibility/v3_source_reconciliation_20260907.json").read_text())
    receipt[section][field] = value
    with pytest.raises(ValueError, match="Source"):
        _check_source_recovery(receipt, source)


def test_tested_or_planned_stage_cannot_claim_the_source_execution_receipt(tmp_path):
    root = clone_inputs(tmp_path)
    path = root / INDEX
    index = json.loads(path.read_text())
    index["requirements"]["joint_covariance_aware_hierarchy"].update(
        state="execution_evidence_recorded", evidence="reproducibility/v3_source_cohort_recovery_20260908.json")
    path.write_text(json.dumps(index))
    with pytest.raises(ValueError, match="Executed evidence requirement"):
        validate(root)


@pytest.mark.parametrize("section,field,value", [
    ("schedule", "all_native_observations_preserved", False),
    ("schedule", "worker_location_taxon_date_fields_present", True),
    ("private_restore", "off_device_private_restore_verified", True),
    ("private_restore", "snapshot_manifest_sha256", "0" * 64),
    ("worker_handoff", "original_and_restored_packets_canonically_identical", False),
    ("partition_audit", "split_known_components", 1),
])
def test_schedule_receipt_cannot_lose_links_or_claim_off_device_validation(section, field, value):
    receipt = json.loads((ROOT / "reproducibility/v3_reconciled_schedule_private_replay_20260908.json").read_text())
    recovery = json.loads((ROOT / "reproducibility/v3_source_cohort_recovery_20260908.json").read_text())
    receipt[section][field] = value
    with pytest.raises(ValueError, match="Schedule replay"):
        _check_schedule_replay(receipt, recovery)


def test_local_replay_cannot_be_relabelled_as_executed_durable_storage(tmp_path):
    root = clone_inputs(tmp_path)
    path = root / INDEX
    index = json.loads(path.read_text())
    index["requirements"]["durable_private_numerical_replay"]["evidence"] = "reproducibility/v3_reconciled_schedule_private_replay_20260908.json"
    path.write_text(json.dumps(index))
    with pytest.raises(ValueError, match="Executed evidence requirement"):
        validate(root)


@pytest.mark.parametrize("section,field,value", [
    ("cloud", "draft_verified", False),
    ("cloud", "anonymous_release_and_asset_requests", "200"),
    ("cloud", "source_files_restored", 871),
    ("cloud", "selected_observations", 127),
    ("local_verification", "selected_observations", 127),
    ("local_verification", "local_and_cloud_packet_canonical_sha256", "0" * 64),
    ("local_verification", "full_original_stream_authorized", True),
])
def test_cloud_recovery_cannot_change_source_or_promote_scientific_readiness(section, field, value):
    receipt = json.loads((ROOT / "reproducibility/v3_protected_numerical_replay_20260908.json").read_text())
    previous = json.loads((ROOT / "reproducibility/v3_reconciled_schedule_private_replay_20260908.json").read_text())
    receipt[section][field] = value
    with pytest.raises(ValueError, match="Protected replay"):
        _check_protected_replay(receipt, previous)
