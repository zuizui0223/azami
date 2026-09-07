"""Freeze the v3 environment representation before any trait join.

This module consumes only aggregate phenotype-blind environment diagnostics plus
an explicit decision JSON. It never reads image, endpoint or ecological-result
files. A freeze receipt is written only when the decision completely partitions
the diagnosed candidate set and the diagnostics came from the native-range,
equal-taxon-weighted, non-pilot environment design.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


FREEZE_STATUS = "V3_ENVIRONMENT_REPRESENTATION_FROZEN_BEFORE_TRAIT_JOIN"
DECISION_STATUS = "phenotype_blind_environment_representation_decision"
ALLOWED_SELECTION_BASIS = {
    "biological_proximity",
    "coverage",
    "environment_only_redundancy",
}
REQUIRED_DIAGNOSTIC_FILES = (
    "environment_diagnostics.json",
    "environment_coverage.csv",
    "environment_correlations.csv",
    "environment_vif.csv",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _variable_rows(items: object, label: str) -> list[dict]:
    if not isinstance(items, list) or not items:
        raise ValueError(f"Decision {label} must be a non-empty list")
    rows: list[dict] = []
    for item in items:
        if not isinstance(item, dict):
            raise ValueError(f"Decision {label} entries must be objects")
        variable = str(item.get("variable") or "").strip()
        rationale = str(item.get("rationale") or "").strip()
        if not variable or not rationale:
            raise ValueError(f"Decision {label} entries require variable and rationale")
        rows.append({**item, "variable": variable, "rationale": rationale})
    variables = [row["variable"] for row in rows]
    if len(variables) != len(set(variables)):
        raise ValueError(f"Decision {label} contains duplicate variables")
    return rows


def load_decision(path: Path) -> dict:
    decision = json.loads(path.read_text(encoding="utf-8"))
    if decision.get("status") != DECISION_STATUS:
        raise ValueError("Environment decision is not in phenotype-blind decision state")
    if decision.get("trait_results_inspected") is not False:
        raise ValueError("Environment decision must explicitly state trait_results_inspected=false")
    basis = set(decision.get("selection_basis") or [])
    if basis != ALLOWED_SELECTION_BASIS:
        raise ValueError(
            "Environment decision selection_basis must be exactly biological_proximity, coverage and environment_only_redundancy"
        )
    try:
        workflow_run_id = int(decision.get("diagnostics_workflow_run_id"))
        artifact_id = int(decision.get("diagnostics_artifact_id"))
    except (TypeError, ValueError):
        raise ValueError("Environment decision requires numeric diagnostics workflow/artifact IDs") from None
    artifact_name = str(decision.get("diagnostics_artifact_name") or "").strip()
    artifact_sha256 = str(decision.get("diagnostics_artifact_sha256") or "").strip().lower()
    if workflow_run_id <= 0 or artifact_id <= 0 or not artifact_name:
        raise ValueError("Environment decision requires positive diagnostic IDs and artifact name")
    if len(artifact_sha256) != 64 or any(ch not in "0123456789abcdef" for ch in artifact_sha256):
        raise ValueError("Environment decision requires a 64-hex diagnostics artifact SHA-256")
    selected = _variable_rows(decision.get("selected"), "selected")
    rejected = _variable_rows(decision.get("rejected"), "rejected")
    selected_names = {row["variable"] for row in selected}
    rejected_names = {row["variable"] for row in rejected}
    if selected_names & rejected_names:
        raise ValueError("Selected and rejected environment variables overlap")
    return {
        **decision,
        "diagnostics_workflow_run_id": workflow_run_id,
        "diagnostics_artifact_id": artifact_id,
        "diagnostics_artifact_name": artifact_name,
        "diagnostics_artifact_sha256": artifact_sha256,
        "selected": selected,
        "rejected": rejected,
    }


def validate_diagnostics(diagnostics_dir: Path) -> tuple[dict, pd.DataFrame]:
    for name in REQUIRED_DIAGNOSTIC_FILES:
        if not (diagnostics_dir / name).is_file():
            raise ValueError(f"Missing environment diagnostic file: {name}")
    report = json.loads((diagnostics_dir / "environment_diagnostics.json").read_text(encoding="utf-8"))
    if report.get("status") != "ENVIRONMENT_ONLY_DIAGNOSTICS_COMPLETE":
        raise ValueError("Environment diagnostics are not complete")
    if int(report.get("trait_columns_read", -1)) != 0:
        raise ValueError("Environment diagnostics are not phenotype-blind")
    if report.get("native_only") is not True:
        raise ValueError("Environment freeze requires native-only diagnostics")
    if report.get("engineering_pilot_only") is not False:
        raise ValueError("Engineering-pilot diagnostics cannot freeze the environment representation")
    if report.get("representation_may_be_frozen_from_this_run") is not True:
        raise ValueError("Diagnostic run does not authorize representation freeze")
    if report.get("weight_column") != "equal_taxon_weight":
        raise ValueError("Environment freeze requires equal-taxon-weighted diagnostics")
    variables = [str(value) for value in report.get("variables") or []]
    if not variables or len(variables) != len(set(variables)):
        raise ValueError("Environment diagnostic candidate variables are empty or duplicated")
    coverage = pd.read_csv(diagnostics_dir / "environment_coverage.csv")
    if "variable" not in coverage.columns:
        raise ValueError("Environment coverage table lacks variable column")
    coverage_variables = set(coverage["variable"].astype(str))
    if coverage_variables != set(variables):
        raise ValueError("Environment coverage variables differ from diagnostic candidate set")
    return report, coverage


def freeze(
    diagnostics_dir: Path,
    decision_path: Path,
    contract_path: Path,
    out_path: Path,
) -> dict:
    report, coverage = validate_diagnostics(diagnostics_dir)
    decision = load_decision(decision_path)
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract.get("status") != "environment_only_design_before_trait_join":
        raise ValueError("Environment exposure contract is not in design-before-trait state")
    contract_id = str(contract.get("contract_id") or "")
    if not contract_id or report.get("contract_id") != contract_id or decision.get("contract_id") != contract_id:
        raise ValueError("Environment contract IDs do not agree across diagnostics and decision")

    candidate_set = set(map(str, report["variables"]))
    selected_set = {row["variable"] for row in decision["selected"]}
    rejected_set = {row["variable"] for row in decision["rejected"]}
    if selected_set | rejected_set != candidate_set:
        missing = sorted(candidate_set - selected_set - rejected_set)
        extra = sorted((selected_set | rejected_set) - candidate_set)
        raise ValueError(f"Environment decision does not exactly partition candidates; missing={missing}, extra={extra}")

    coverage_lookup = coverage.set_index(coverage["variable"].astype(str))["weighted_coverage"].to_dict()
    diagnostic_hashes = {name: sha256_file(diagnostics_dir / name) for name in REQUIRED_DIAGNOSTIC_FILES}
    receipt = {
        "status": FREEZE_STATUS,
        "contract_id": contract_id,
        "decision_status": decision["status"],
        "decision_scope": decision.get("decision_scope", "native_range_taxon_balanced_environment_only"),
        "selection_basis": sorted(ALLOWED_SELECTION_BASIS),
        "diagnostics_provenance": {
            "workflow_run_id": decision["diagnostics_workflow_run_id"],
            "artifact_id": decision["diagnostics_artifact_id"],
            "artifact_name": decision["diagnostics_artifact_name"],
            "artifact_sha256": decision["diagnostics_artifact_sha256"],
        },
        "selected": decision["selected"],
        "rejected": decision["rejected"],
        "selected_variables": [row["variable"] for row in decision["selected"]],
        "rejected_variables": [row["variable"] for row in decision["rejected"]],
        "candidate_variables": list(report["variables"]),
        "candidate_partition_complete": True,
        "weighted_coverage": {variable: float(coverage_lookup[variable]) for variable in report["variables"]},
        "matrix": report.get("matrix"),
        "redundancy_components": report.get("redundancy_components"),
        "absolute_correlation_redundancy_threshold": report.get("absolute_correlation_redundancy_threshold"),
        "diagnostic_file_sha256": diagnostic_hashes,
        "decision_sha256": sha256_file(decision_path),
        "contract_sha256": sha256_file(contract_path),
        "trait_files_read": 0,
        "image_files_read": 0,
        "ecological_models_executed": 0,
        "trait_join_authorized_after_this_receipt": True,
        "claim_boundary": "This receipt freezes exposure representation only; it contains no capitulum-trait association result.",
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        raise ValueError("Freeze receipt already exists; preserve prior decision")
    out_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics-dir", type=Path, required=True)
    parser.add_argument("--decision", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=Path("analysis/v3/environment_exposure_contract.json"))
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = freeze(args.diagnostics_dir, args.decision, args.contract, args.out)
    print(json.dumps({"status": receipt["status"], "selected_variables": receipt["selected_variables"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
