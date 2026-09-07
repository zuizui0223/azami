import json
from pathlib import Path

import pandas as pd
import pytest

from analysis.v3.freeze_environment_representation import FREEZE_STATUS, freeze


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "environment_exposure_contract.json"
ARTIFACT_SHA = "a" * 64


def write_diagnostics(root: Path, variables: list[str]) -> str:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    report = {
        "status": "ENVIRONMENT_ONLY_DIAGNOSTICS_COMPLETE",
        "contract_id": contract["contract_id"],
        "trait_columns_read": 0,
        "native_only": True,
        "engineering_pilot_only": False,
        "representation_may_be_frozen_from_this_run": True,
        "weight_column": "equal_taxon_weight",
        "variables": variables,
        "matrix": {"matrix_rank": len(variables), "condition_number": 2.0},
        "redundancy_components": [variables],
        "absolute_correlation_redundancy_threshold": 0.8,
    }
    (root / "environment_diagnostics.json").write_text(json.dumps(report), encoding="utf-8")
    pd.DataFrame(
        {"variable": variables, "weighted_coverage": [1.0] * len(variables)}
    ).to_csv(root / "environment_coverage.csv", index=False)
    pd.DataFrame(
        {"method": [], "variable_a": [], "variable_b": [], "correlation": []}
    ).to_csv(root / "environment_correlations.csv", index=False)
    pd.DataFrame({"variable": variables, "vif": [1.0] * len(variables)}).to_csv(
        root / "environment_vif.csv", index=False
    )
    return contract["contract_id"]


def decision_base(contract_id: str) -> dict:
    return {
        "status": "phenotype_blind_environment_representation_decision",
        "contract_id": contract_id,
        "decision_scope": "native_range_taxon_balanced_environment_only",
        "trait_results_inspected": False,
        "diagnostics_workflow_run_id": 123456,
        "diagnostics_artifact_id": 789012,
        "diagnostics_artifact_name": "ch1-v3-phenotype-blind-environment-design-123456",
        "diagnostics_artifact_sha256": ARTIFACT_SHA,
        "selection_basis": [
            "biological_proximity",
            "coverage",
            "environment_only_redundancy",
        ],
    }


def test_freeze_requires_complete_phenotype_blind_candidate_partition(tmp_path: Path):
    diagnostics = tmp_path / "diagnostics"
    diagnostics.mkdir()
    contract_id = write_diagnostics(diagnostics, ["pr_month", "BIO12"])
    decision = {
        **decision_base(contract_id),
        "selected": [
            {
                "variable": "pr_month",
                "rationale": "Direct observation-month wetting exposure retained over a broader redundant summary.",
            }
        ],
        "rejected": [
            {
                "variable": "BIO12",
                "rationale": "Broader annual precipitation representation rejected as redundant with the direct exposure.",
            }
        ],
    }
    decision_path = tmp_path / "decision.json"
    decision_path.write_text(json.dumps(decision), encoding="utf-8")
    out = tmp_path / "freeze.json"

    receipt = freeze(diagnostics, decision_path, CONTRACT, out)

    assert receipt["status"] == FREEZE_STATUS
    assert receipt["selected_variables"] == ["pr_month"]
    assert receipt["rejected_variables"] == ["BIO12"]
    assert receipt["candidate_partition_complete"] is True
    assert receipt["diagnostics_provenance"]["workflow_run_id"] == 123456
    assert receipt["diagnostics_provenance"]["artifact_id"] == 789012
    assert receipt["diagnostics_provenance"]["artifact_sha256"] == ARTIFACT_SHA
    assert receipt["trait_files_read"] == 0
    assert receipt["ecological_models_executed"] == 0
    assert out.is_file()


def test_freeze_rejects_incomplete_candidate_partition(tmp_path: Path):
    diagnostics = tmp_path / "diagnostics"
    diagnostics.mkdir()
    contract_id = write_diagnostics(diagnostics, ["pr_month", "BIO12", "BIO18"])
    decision = {
        **decision_base(contract_id),
        "selected": [{"variable": "pr_month", "rationale": "Direct wetting exposure."}],
        "rejected": [{"variable": "BIO12", "rationale": "Broader redundant summary."}],
    }
    decision_path = tmp_path / "decision.json"
    decision_path.write_text(json.dumps(decision), encoding="utf-8")

    with pytest.raises(ValueError, match="does not exactly partition candidates"):
        freeze(diagnostics, decision_path, CONTRACT, tmp_path / "freeze.json")


def test_freeze_rejects_missing_diagnostic_artifact_provenance(tmp_path: Path):
    diagnostics = tmp_path / "diagnostics"
    diagnostics.mkdir()
    contract_id = write_diagnostics(diagnostics, ["pr_month", "BIO12"])
    decision = {
        **decision_base(contract_id),
        "selected": [{"variable": "pr_month", "rationale": "Direct wetting exposure."}],
        "rejected": [{"variable": "BIO12", "rationale": "Broader redundant summary."}],
    }
    decision.pop("diagnostics_artifact_sha256")
    decision_path = tmp_path / "decision.json"
    decision_path.write_text(json.dumps(decision), encoding="utf-8")

    with pytest.raises(ValueError, match="64-hex diagnostics artifact SHA-256"):
        freeze(diagnostics, decision_path, CONTRACT, tmp_path / "freeze.json")
