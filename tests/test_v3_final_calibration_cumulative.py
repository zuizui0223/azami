import json
from pathlib import Path

import pytest

from analysis.v3.summarize_final_calibration_cumulative import (
    CONTINUE_STATUS,
    MAX_FAIL_STATUS,
    QUALIFIED_STATUS,
    OUTERS_REQUIRED_STATUS,
    aggregate,
)
from analysis.v3.workflow import ROOT, canonical_digest


def _contract():
    return json.loads((ROOT / "analysis/v3/final_module_calibration_contract.json").read_text(encoding="utf-8"))


def _write_batch(tmp_path: Path, batch_index: int, *, false_outer: int = 0,
                 covered: int = 8775, total: int = 8775):
    contract = _contract()
    batch_size = contract["sequential_outer_rule"]["batch_size_per_scenario"]
    start = batch_index * batch_size
    stop = start + batch_size
    scenarios = contract["scenarios"]
    grids = contract["grid_degrees"]
    report = {
        "schema_version": 1,
        "status": "FINAL_MULTICOORDINATE_CALIBRATION_BATCH_COMPLETE_CONTINUE_SEQUENTIAL_RULE",
        "contract_id": contract["id"],
        "contract_canonical_sha256": canonical_digest(contract),
        "batch_index": batch_index,
        "outer_start_inclusive": start,
        "outer_stop_exclusive": stop,
        "expected_outer_reports": len(scenarios) * batch_size,
        "recorded_outer_reports": len(scenarios) * batch_size,
        "missing_outer_reports": [],
        "expected_grid_cases": len(scenarios) * batch_size * len(grids),
        "recorded_grid_cases": len(scenarios) * batch_size * len(grids),
        "complete_grid_cases": len(scenarios) * batch_size * len(grids),
        "scenario_grid_summary": [],
        "input_manifest": [],
        # Deliberately untrusted by the cumulative evaluator. A batch-local
        # decision must never substitute for the cumulative sequential rule.
        "batch_can_authorize_stopping": True,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    for scenario in scenarios:
        for grid in grids:
            report["scenario_grid_summary"].append({
                "scenario": scenario,
                "grid_degrees": grid,
                "expected_outer_reports": batch_size,
                "recorded_outer_reports": batch_size,
                "complete_outer_reports": batch_size,
                "outer_with_any_false_holm_rejection": false_outer,
                "false_holm_rejections": false_outer,
                "true_holm_rejections": 20,
                "coefficient_coverage": {
                    "covered": covered,
                    "total": total,
                    "fraction": covered / total,
                },
                "sequential_precision_decision": {
                    "precision_satisfied": True,
                    "admission_satisfied": True,
                    "stop_for_precision_and_admission": True,
                },
            })
        for outer in range(start, stop):
            report["input_manifest"].append({
                "scenario": scenario,
                "outer_replicate": outer,
                "path": f"synthetic/{scenario}-{outer}.json",
                "sha256": f"{outer + 100000 * scenarios.index(scenario):064x}"[-64:],
            })
    path = tmp_path / f"batch-{batch_index}.json"
    path.write_text(json.dumps(report, sort_keys=True) + "\n", encoding="utf-8")
    return path, report


def test_counts_alone_cannot_authorize_stopping_or_next_batch(tmp_path):
    paths = [_write_batch(tmp_path, i)[0] for i in range(3)]
    result = aggregate(paths)
    assert result["status"] == OUTERS_REQUIRED_STATUS
    assert result["cumulative_outer_replicates_per_scenario"] == 75
    assert result["minimum_reached"] is False
    assert result["all_scenario_grid_precision_satisfied"] is False
    assert result["synthetic_calibration_qualified_for_realized_design_gate"] is False
    assert result["batch_local_stopping_decisions_used"] is False
    assert result["next_frozen_batch"] is None
    assert result["ecological_fitting_authorized"] is False


def test_former_count_only_qualification_is_now_explicitly_disabled(tmp_path):
    paths = [_write_batch(tmp_path, i)[0] for i in range(4)]
    result = aggregate(paths)
    assert result["status"] == OUTERS_REQUIRED_STATUS
    assert result["included_batch_indices"] == [0, 1, 2, 3]
    assert result["cumulative_outer_replicates_per_scenario"] == 100
    assert result["minimum_reached"] is True
    assert result["all_scenario_grid_precision_satisfied"] is False
    assert result["all_scenario_grid_admission_satisfied"] is False
    assert result["synthetic_calibration_qualified_for_realized_design_gate"] is False
    assert result["next_frozen_batch"] is None
    assert len(result["scenario_grid_cumulative_summary"]) == 8
    assert all(row["cumulative_outer_reports"] == 100 for row in result["scenario_grid_cumulative_summary"])
    assert result["ecological_fitting_authorized"] is False


def test_missing_intermediate_batch_is_rejected_not_silently_treated_as_cumulative(tmp_path):
    batch0 = _write_batch(tmp_path, 0)[0]
    batch2 = _write_batch(tmp_path, 2)[0]
    with pytest.raises(ValueError, match="contiguous from batch 0"):
        aggregate([batch0, batch2])


def test_duplicate_batch_index_is_rejected(tmp_path):
    batch0, report = _write_batch(tmp_path, 0)
    duplicate = tmp_path / "batch-0-copy.json"
    duplicate.write_text(json.dumps(report, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate batch index 0"):
        aggregate([batch0, duplicate])


def test_completed_batch_must_carry_exact_outer_manifest_and_all_grid_cases(tmp_path):
    _, report = _write_batch(tmp_path, 0)
    report["input_manifest"].pop()
    broken_manifest = tmp_path / "broken-manifest.json"
    broken_manifest.write_text(json.dumps(report, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="exact frozen outer range"):
        aggregate([broken_manifest])

    _, report = _write_batch(tmp_path, 0)
    report["complete_grid_cases"] -= 1
    broken_grid = tmp_path / "broken-grid.json"
    broken_grid.write_text(json.dumps(report, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="incomplete grid cases"):
        aggregate([broken_grid])


def _mock_verified_outers(monkeypatch, *, count, false_per_batch=0):
    # Test aggregation independently of disk validation, covered in the outer
    # verification tests. Fixed final geometry: 351 coefficients per outer.
    values = {(s, i): {g: {"false_family": int(i % 25 < false_per_batch), "covered": 351, "total": 351,
                          "coverage_fraction": 1., "false_holm_rejections": int(i % 25 < false_per_batch),
                          "true_holm_rejections": int(i % 25 < 20)} for g in (2, 5)}
              for s in _contract()["scenarios"] for i in range(count)}
    monkeypatch.setattr("analysis.v3.summarize_final_calibration_cumulative.verify_inventory",
                        lambda *a, **kw: (values, []))


def test_maximum_n400_with_failed_admission_stops_without_rescue_or_ecology(tmp_path, monkeypatch):
    paths = [_write_batch(tmp_path, i, false_outer=5)[0] for i in range(16)]
    _mock_verified_outers(monkeypatch, count=400, false_per_batch=5)
    result = aggregate(paths, outer_report_paths=["mock-only"])
    assert result["status"] == MAX_FAIL_STATUS
    assert result["cumulative_outer_replicates_per_scenario"] == 400
    assert result["all_scenario_grid_admission_satisfied"] is False
    assert result["synthetic_calibration_qualified_for_realized_design_gate"] is False
    assert result["next_frozen_batch"] is None
    assert result["thresholds_changed"] is False
    assert result["seed_changed"] is False
    assert result["grid_changed"] is False
    assert result["family_changed"] is False
    assert result["synthetic_truth_changed"] is False
    assert result["ecological_fitting_authorized"] is False


def test_verified_outers_can_qualify_only_realized_design_gate(tmp_path, monkeypatch):
    paths = [_write_batch(tmp_path, i)[0] for i in range(16)]
    _mock_verified_outers(monkeypatch, count=400)
    result = aggregate(paths, outer_report_paths=["mock-only"])
    assert result["status"] == QUALIFIED_STATUS
    assert result["synthetic_calibration_qualified_for_realized_design_gate"] is True
    assert result["ecological_fitting_authorized"] is False
    assert result["next_frozen_batch"] is None
    assert all(r["cumulative_precision_decision"]["coverage_confidence_sequence"]["independent_outer_datasets"] == 400
               for r in result["scenario_grid_cumulative_summary"])


def test_verified_small_batch_continues_but_count_disagreement_stops(tmp_path, monkeypatch):
    path, report = _write_batch(tmp_path, 0)
    _mock_verified_outers(monkeypatch, count=25)
    result = aggregate([path], outer_report_paths=["mock-only"])
    assert result["status"] == CONTINUE_STATUS
    assert result["next_frozen_batch"]["outer_start_inclusive"] == 25
    report["scenario_grid_summary"][0]["coefficient_coverage"]["covered"] -= 1
    path.write_text(json.dumps(report), encoding="utf-8")
    with pytest.raises(ValueError, match="batch counts differ"):
        aggregate([path], outer_report_paths=["mock-only"])
