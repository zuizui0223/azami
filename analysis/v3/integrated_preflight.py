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
    required_ids = {"source_receipt", "environment_receipt", "measurement_receipt", "measurement_decision"}
    if set(ids) != required_ids or len(ids) != len(set(ids)) or len(paths) != len(set(paths)):
        raise ValueError("Evidence index has missing or duplicate required identities/paths")
    for item in evidence["inputs"]:
        data = _read(root, item["path"])
        actual = canonical_digest(data)
        if actual != item["canonical_json_sha256"]:
            raise ValueError(f"Pinned evidence changed: {item['path']}")
        loaded[item["id"]] = data
        checked.append({"id": item["id"], "path": item["path"], "canonical_json_sha256": actual})
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
            if spec["state"] not in {"historical_evidence_recorded", "implementation_tested_execution_pending", "not_executed"}:
                raise ValueError("Readiness cannot be asserted by an arbitrary status string")
            if spec["state"] == "historical_evidence_recorded" and spec.get("evidence") not in paths:
                raise ValueError("Historical evidence requirement is not bound to a verified input")
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
        "next_source_gate": "Recover the exact 319244-row native cohort CSV and authority chain, then execute the enriched source builder; counts alone are insufficient.",
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
