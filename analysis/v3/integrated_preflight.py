"""Read the active v3 design and verify pinned public evidence without promotion.

This command has no image retrieval, trait join or fitting path. Successful exit
means the design/evidence index is internally consistent, not production readiness.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from .workflow import ROOT, canonical_digest
from .environment_model import definition as environment_definition

CONTRACT = "analysis/v3/integrated_workflow_contract.json"
INDEX = "analysis/v3/integrated_evidence_index.json"
SOURCE_RECOVERY_STATUS = "EXACT_SOURCE_COHORT_RECOVERED_ENRICHED_AND_LOCAL_AUTHORITY_REPLAY_VERIFIED"


def _check_protected_replay(receipt: dict, previous: dict) -> None:
    """Bind cloud recovery to the prior local source without model promotion."""
    cloud, local = receipt["cloud"], receipt["local_verification"]
    source, handoff = previous["private_restore"], previous["worker_handoff"]
    if (receipt["status"] != "PROTECTED_ACTIONS_NUMERICAL_REPLAY_AND_LOCAL_RECOVERY_VERIFIED"
            or cloud["status"] != "GITHUB_DRAFT_NUMERICAL_RESTORE_AND_OUTPUT_ROUNDTRIP_VERIFIED"
            or local["status"] != "ACTIONS_SOURCE_RESTORE_AND_LOCAL_RETURNED_PACKET_VERIFIED"
            or cloud["draft_verified"] is not True
            or cloud["anonymous_release_and_asset_requests"] != "404"
            or cloud["repository"] != "zuizui0223/azami"):
        raise ValueError("Protected replay is not a verified private roundtrip")
    if (cloud["run_id"] != local["cloud_run_id"] or cloud["commit"] != local["cloud_commit"]
            or cloud["contract_sha256"] != local["contract_sha256"]
            or cloud["output_asset"] != local["returned_asset"]
            or cloud["source_files_restored"] != local["source_files_restored_on_actions"]
            or cloud["source_files_restored"] != source["files"]
            or cloud["source_bytes_restored"] != local["source_bytes_restored_on_actions"]
            or cloud["source_bytes_restored"] != source["restored_bytes"]
            or local["source_asset"]["manifest_sha256"] != source["snapshot_manifest_sha256"]
            or local["local_and_cloud_packet_canonical_sha256"] != handoff["packet_canonical_sha256"]
            or not local["selected_observations"] == cloud["selected_observations"] == handoff["selected_observations"]
            or local["request_candidates_not_executed"] != handoff["request_candidates"]):
        raise ValueError("Protected replay source, cloud or local identity differs")
    if (any(cloud[key] != 0 for key in ("source_images_persisted", "image_requests_executed", "ecological_models_executed"))
            or any(local[key] != 0 for key in ("source_files_deleted", "original_images_persisted", "ecological_models_executed"))
            or cloud["production_image_execution_authorized"] is not False
            or local["full_original_stream_authorized"] is not False
            or cloud["ecological_fitting_authorized"] is not False
            or local["ecological_fitting_authorized"] is not False):
        raise ValueError("Protected replay cannot authorize production or ecology")


def _check_schedule_replay(receipt: dict, recovery: dict) -> None:
    """Validate source conservation and local replay without a durability claim."""
    if receipt.get("status") != "RECONCILED_SCHEDULE_AND_LOCAL_PRIVATE_REPLAY_VERIFIED_OFF_DEVICE_PENDING":
        raise ValueError("Schedule replay status is not the bounded local result")
    schedule, restored = receipt["schedule"], receipt["private_restore"]
    packed, handoff = receipt["private_snapshot"], receipt["worker_handoff"]
    counts, enrichment = schedule["counts"], recovery["enrichment"]
    if (schedule["input_sha256"]["enriched"] != enrichment["output_csv_sha256"]
            or schedule["input_sha256"]["reconciliation"] != enrichment["input_sha256"]["reconciliation"]
            or counts["native_observations"] != recovery["cohort"]["rows"]
            or schedule["native_known_components"] != enrichment["dependence"]["cohort_known_components"]
            or not counts["photo_versions"] >= counts["known_photo_links"] >= counts["native_links"] >= counts["photo_jobs"] > 0
            or sum(schedule["photo_states"].values()) != counts["photo_jobs"]):
        raise ValueError("Schedule replay source identity or count conservation differs")
    if (any(schedule[key] is not True for key in ("all_native_observations_preserved", "all_known_target_photo_links_preserved", "shared_photos_partitioned_once"))
            or schedule["worker_location_taxon_date_fields_present"] is not False
            or schedule["status"] != "RECONCILED_NATIVE_PHOTO_SCHEDULE_VERIFIED_NO_IMAGE_EXECUTION"):
        raise ValueError("Schedule replay changed preservation or worker boundaries")
    partition = receipt["partition_audit"]
    if (sum(partition["shard_counts"].values()) != counts["photo_jobs"]
            or len(partition["shard_counts"]) != schedule["shard_count"]
            or partition["split_known_components"] != 0 or partition["all_jobs_accounted_for"] is not True
            or handoff["original_and_restored_packets_canonically_identical"] is not True
            or not 0 < handoff["selected_observations"] <= handoff["maximum_observations"] <= 128
            or handoff["request_candidates"] != handoff["photo_states"].get("request_candidate_not_authorized", 0)
            or sum(handoff["photo_states"].values()) != handoff["selected_photo_jobs"]):
        raise ValueError("Schedule replay partition or worker handoff differs")
    if (packed["status"] != "PRIVATE_LOCAL_SNAPSHOT_WRITTEN_AND_HASH_VERIFIED"
            or restored["status"] != "PRIVATE_LOCAL_RESTORE_ALL_MEMBERS_BYTE_VERIFIED"
            or packed["snapshot_manifest_sha256"] != restored["snapshot_manifest_sha256"]
            or not packed["files"] == restored["files"] == receipt["private_bundle_scope"]["files"] > 0
            or packed["unique_blob_bytes"] > restored["restored_bytes"]):
        raise ValueError("Schedule replay private restoration evidence differs")
    for value in (receipt, packed, restored):
        if (value["off_device_private_restore_verified"] is not False
                or value["production_image_execution_authorized"] is not False
                or value["ecological_fitting_authorized"] is not False):
            raise ValueError("Schedule replay cannot promote off-device or production claims")
    if (receipt["source_files_deleted"] != 0 or receipt["image_requests_executed"] != 0
            or receipt["trait_values_read"] != 0 or receipt["ecological_models_executed"] != 0
            or schedule["image_requests_executed"] != 0 or handoff["images_fetched"] != 0
            or schedule["production_image_execution_authorized"] is not False
            or handoff["production_image_execution_authorized"] is not False):
        raise ValueError("Schedule replay is source-only, not image/ecology execution")


def _check_source_recovery(receipt: dict, source: dict) -> None:
    """Check a public executed receipt without implying access to private files."""
    if receipt.get("status") != SOURCE_RECOVERY_STATUS:
        raise ValueError("Source recovery receipt status is not verified")
    cohort = receipt["cohort"]
    enrichment = receipt["enrichment"]
    authority = receipt["authority"]
    exact = "b52503cd891313daaa83dacddc5be92bd781bdbbb2e96b3b75da740239816314"
    if (cohort["historical_csv_sha256"] != exact or enrichment["input_sha256"]["cohort"] != exact
            or cohort["exact_historical_csv_identity_verified"] is not True
            or cohort["rows"] != 319244 or cohort["accepted_taxa"] != 355
            or enrichment["cohort_rows"] != cohort["rows"]
            or enrichment["accepted_taxa"] != cohort["accepted_taxa"]):
        raise ValueError("Source recovery does not preserve the exact historical cohort")
    if (receipt["source"]["serialization_recovery"]["output_sha256"] !=
            "5632f532a63c8babdc20b023bd0d3c47424d69461df5de9793028e93e819f959"
            or enrichment["dependence"]["full_source_observations"] != source["counts"]["retained_unique_observations"]
            or enrichment["dependence"]["reconciled_source_photo_links"] != source["counts"]["retained_unique_observation_photo_links"]):
        raise ValueError("Source recovery differs from the source identity/count chain")
    if (authority["all_three_outputs_byte_identical_in_offline_replay"] is not True
            or authority["historical_http_response_identity_verified"] is not False
            or authority["additional_search_pages_affect_eligible_source_ranks"] is not False
            or set(authority["offline_replay_outputs_sha256"]) != {
                "native_range_join.csv", "taxon_resolution.csv", "wcvp_distribution_records.csv"}):
        raise ValueError("Source authority replay or historical evidence boundary changed")
    if (enrichment["cohort_membership_changed"] is not False or enrichment["source_rows_deleted"] != 0
            or receipt["independent_saved_csv_audit"]["exact_membership_and_all_12_inherited_fields_preserved"] is not True
            or enrichment["calendar"]["exact_date_rows"] != cohort["rows"]
            or sum(enrichment["calendar"]["hemisphere_counts"].values()) != cohort["rows"]):
        raise ValueError("Source enrichment membership or calendar accounting changed")
    if (receipt["ecological_fitting_authorized"] is not False
            or receipt["production_image_execution_authorized"] is not False
            or enrichment["ecological_fitting_authorized"] is not False
            or receipt["ecological_models_executed"] != 0
            or receipt["trait_files_read"] != 0 or receipt["environment_values_read"] != 0
            or receipt["private_preservation"]["cloud_execution_hold_preserved"] is not True):
        raise ValueError("Source recovery cannot promote measurement/ecological execution")


def _read(root: Path, relative: str) -> dict:
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("Evidence path leaves repository")
    return json.loads(path.read_text(encoding="utf-8"))


def validate(root: Path = ROOT) -> dict:
    contract = _read(root, CONTRACT)
    evidence = _read(root, INDEX)
    if contract.get("status") != "active_design_with_execution_requirements_unfinished":
        raise ValueError("This preflight cannot promote a changed execution authority")
    stage_ids = [stage["id"] for stage in contract["stages"]]
    if stage_ids != ["source", "measurement_qc", "ecology", "output"]:
        raise ValueError("Active stage ordering changed")
    if contract["hypervolume"]["optional"] is not True:
        raise ValueError("Noncanonical hypervolume must remain optional")

    loaded = {}
    checked = []
    ids = [item["id"] for item in evidence["inputs"]]
    paths = [item["path"] for item in evidence["inputs"]]
    required_ids = {"source_receipt", "environment_receipt", "measurement_receipt", "measurement_decision", "source_recovery_receipt", "schedule_replay_receipt", "protected_replay_receipt"}
    if set(ids) != required_ids or len(ids) != len(set(ids)) or len(paths) != len(set(paths)):
        raise ValueError("Evidence index has missing or duplicate required identities/paths")
    for item in evidence["inputs"]:
        data = _read(root, item["path"])
        actual = canonical_digest(data)
        if actual != item["canonical_json_sha256"]:
            raise ValueError(f"Pinned evidence changed: {item['path']}")
        loaded[item["id"]] = data
        checked.append({"id": item["id"], "path": item["path"], "canonical_json_sha256": actual})

    _check_source_recovery(loaded["source_recovery_receipt"], loaded["source_receipt"])
    _check_schedule_replay(loaded["schedule_replay_receipt"], loaded["source_recovery_receipt"])
    _check_protected_replay(loaded["protected_replay_receipt"], loaded["schedule_replay_receipt"])

    measurement = loaded["measurement_receipt"]
    environment = loaded["environment_receipt"]
    decision = loaded["measurement_decision"]
    design = _read(root, "analysis/v3/ecological_analysis_contract.json")
    gates = design["completed_design_gates"]
    if gates["measurement_freeze_status"] != measurement["status"]:
        raise ValueError("Measurement status differs from the actual receipt")
    if gates["environment_freeze_status"] != environment["status"]:
        raise ValueError("Environment status differs from the actual receipt")
    if decision["followup_provenance"]["report_sha256"] != measurement["followup_provenance"]["report_sha256"]:
        raise ValueError("Measurement decision and receipt refer to different evidence")
    routes = Counter(row["ecological_route"] for row in decision["endpoints"])
    if dict(routes) != measurement["ecological_route_counts"]:
        raise ValueError("Measurement route counts disagree")
    if len(decision["endpoints"]) != measurement["endpoint_count"]:
        raise ValueError("Endpoint denominator disagrees")

    active_environment = environment_definition(root)
    active = contract["environment_and_inference"]
    if (active["retained_variables"] != active_environment["variables"]
            or active["process_blocks"] != active_environment["processes"]
            or active["environment_selection_receipt"] != active_environment["receipt"]
            or active["environment_selection_receipt_canonical_sha256"] != active_environment["receipt_canonical_sha256"]
            or "drying" in active or "thermal" in active):
        raise ValueError("Integrated exposure formulation differs from the verified source-QC process model")

    requirements = evidence["requirements"]
    expected = {key for stage in contract["stages"] for key in stage["required"]}
    if not expected.issubset(set(requirements)):
        raise ValueError("Readiness index is missing an active four-stage requirement")

    completed_states = {"historical_evidence_recorded", "execution_evidence_recorded"}
    stage_reports = []
    prior_complete = True
    next_stage = None
    for stage in contract["stages"]:
        items = []
        own_complete = True
        for key in stage["required"]:
            spec = requirements[key]
            if spec["state"] not in {"historical_evidence_recorded", "execution_evidence_recorded", "local_execution_recorded_off_device_pending", "implementation_tested_execution_pending", "not_executed"}:
                raise ValueError("Readiness cannot be asserted by an arbitrary status string")
            if spec["state"] == "historical_evidence_recorded" and spec.get("evidence") not in paths:
                raise ValueError("Historical evidence requirement is not bound to a verified input")
            bound = {"native_authority_input_chain": "source_recovery_receipt", "enriched_source_cohort": "source_recovery_receipt",
                     "reconciled_stream_schedule": "schedule_replay_receipt", "durable_private_numerical_replay": "protected_replay_receipt"}
            if spec["state"] in {"execution_evidence_recorded", "local_execution_recorded_off_device_pending"}:
                expected_state = "execution_evidence_recorded"
                if (key not in bound or spec["state"] != expected_state
                        or spec.get("evidence") != next(item["path"] for item in evidence["inputs"] if item["id"] == bound[key])):
                    raise ValueError("Executed evidence requirement is not bound to its verified bounded receipt")
            if spec.get("implementation") and not (root / spec["implementation"]).is_file():
                raise ValueError(f"Indexed implementation absent: {key}")
            if spec["state"] not in completed_states:
                own_complete = False
            items.append({"id": key, **spec})

        if stage["id"] == "output":
            own_complete = prior_complete
        stage_state = "complete" if prior_complete and own_complete else ("blocked" if not prior_complete else "incomplete")
        if next_stage is None and stage_state != "complete":
            next_stage = stage["id"]
        stage_reports.append({
            "stage": stage["id"],
            "state": stage_state,
            "requirements": items,
        })
        prior_complete = prior_complete and own_complete

    if next_stage is None:
        next_stage = "done"

    return {
        "schema_version": 1,
        "status": "FOUR_STAGE_WORKFLOW_INTEGRITY_VERIFIED_EXECUTION_INCOMPLETE" if next_stage != "done" else "FOUR_STAGE_WORKFLOW_COMPLETE",
        "contract_sha256_canonical_json": canonical_digest(contract),
        "evidence_index_sha256_canonical_json": canonical_digest(evidence),
        "verified_public_evidence": checked,
        "active_environment_model": active_environment,
        "canonical_flow": ["source", "measurement_qc", "ecology", "output"],
        "stages": stage_reports,
        "next_stage": next_stage,
        "noncanonical_requirements_retained_for_history": sorted(set(requirements) - expected),
        "ecological_fitting_authorized": False,
        "full_original_stream_authorized": False,
        "trait_values_read": 0,
        "ecological_models_executed": 0,
        "limits": [
            "Pinned public receipts establish their identity and agreement, not present availability or re-verification of every private numerical input.",
            "Implementation tests and design declarations are not executed full-source coverage or ecological evidence.",
            "Hypervolume/breadth synthesis is outside the canonical four-stage path and cannot block output.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = validate()
    if args.out:
        if args.out.exists():
            raise ValueError("Preserve prior preflight output")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
