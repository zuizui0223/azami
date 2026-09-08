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

CONTRACT = "analysis/v3/integrated_workflow_contract.json"
INDEX = "analysis/v3/integrated_evidence_index.json"
SOURCE_RECOVERY_STATUS = "EXACT_SOURCE_COHORT_RECOVERED_ENRICHED_AND_LOCAL_AUTHORITY_REPLAY_VERIFIED"


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
    if stage_ids != ["source", "measurement", "assessability", "ecology", "synthesis"]:
        raise ValueError("Active stage ordering changed")
    if contract["hypervolume"]["optional"] is not True:
        raise ValueError("Secondary hypervolume must not block primary ecology")
    loaded = {}
    checked = []
    ids = [item["id"] for item in evidence["inputs"]]
    paths = [item["path"] for item in evidence["inputs"]]
    required_ids = {"source_receipt", "environment_receipt", "measurement_receipt", "measurement_decision", "source_recovery_receipt"}
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
    for name in ("drying", "thermal"):
        variables = contract["environment_and_inference"][name]
        if variables != environment["integrated_backbones"][name]:
            raise ValueError("Integrated exposure formulation differs from the frozen environment representation")
        if "vpd_month" in variables and "tasmax_month" in variables:
            raise ValueError("Collinear alternatives cannot be co-entered")
    requirements = evidence["requirements"]
    expected = {key for stage in contract["stages"] for key in stage["required"]}
    if expected != set(requirements):
        raise ValueError("Readiness index does not cover exactly the active requirements")
    stage_reports = []
    for stage in contract["stages"]:
        items = []
        for key in stage["required"]:
            spec = requirements[key]
            if spec["state"] not in {"historical_evidence_recorded", "execution_evidence_recorded", "implementation_tested_execution_pending", "not_executed"}:
                raise ValueError("Readiness cannot be asserted by an arbitrary status string")
            if spec["state"] == "historical_evidence_recorded" and spec.get("evidence") not in paths:
                raise ValueError("Historical evidence requirement is not bound to a verified input")
            if spec["state"] == "execution_evidence_recorded" and (
                    key not in {"native_authority_input_chain", "enriched_source_cohort"}
                    or spec.get("evidence") != next(item["path"] for item in evidence["inputs"] if item["id"] == "source_recovery_receipt")):
                raise ValueError("Executed evidence requirement is not bound to its verified source receipt")
            if spec.get("implementation") and not (root / spec["implementation"]).is_file():
                raise ValueError(f"Indexed implementation absent: {key}")
            items.append({"id": key, **spec})
        stage_reports.append({"stage": stage["id"], "optional_for_primary_ecology": stage["id"] == "synthesis", "requirements": items})
    return {
        "schema_version": 1,
        "status": "INTEGRATED_DESIGN_INTEGRITY_VERIFIED_EXECUTION_INCOMPLETE",
        "contract_sha256_canonical_json": canonical_digest(contract),
        "evidence_index_sha256_canonical_json": canonical_digest(evidence),
        "verified_public_evidence": checked,
        "stages": stage_reports,
        "next_source_gate": "Use the recovered/enriched exact 319244-row native cohort in reconciled photo scheduling; verify off-device private preservation and freeze calendar/support handling before production. Historical HTTP bytes remain unavailable, distinct from the verified new offline replay.",
        "ecological_fitting_authorized": False,
        "full_original_stream_authorized": False,
        "trait_values_read": 0, "ecological_models_executed": 0,
        "limits": [
            "Pinned public receipts establish their identity and agreement, not present availability or re-verification of every private numerical input.",
            "Implementation tests and design declarations are not executed full-source coverage or ecological evidence.",
            "The optional synthesis stage is not a primary-ecology prerequisite.",
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
