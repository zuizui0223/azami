import json

import numpy as np
import pytest

from analysis.v3.module_ecology import summarize_draws
from analysis.v3.simulate_full_family_tail_pilot import _slot_rows
from analysis.v3.verify_calibration_outers import IMPLEMENTATIONS, verify_inventory, verify_outer
from analysis.v3.workflow import ROOT, canonical_digest, digest, text_digest


def _json(path, value):
    path.write_text(json.dumps(value, allow_nan=False), encoding="utf-8")


@pytest.fixture
def outer(tmp_path):
    spec = json.loads((ROOT / "analysis/v3/final_module_calibration_contract.json").read_text())
    blocks = json.loads((ROOT / "analysis/v3/crossed_bootstrap_simulation_contract.json").read_text())["process_indices"]
    execution = {"specification_canonical_sha256": canonical_digest(spec), "scenario": "iid_null",
                 "outer_replicate": 0, "outer_seed": spec["seed"],
                 "implementation_sha256_text_lf": {p: text_digest(ROOT / p) for p in IMPLEMENTATIONS},
                 "empirical_trait_environment_values_read": 0, "ecological_models_executed": 0,
                 "ecological_fitting_authorized": False}
    _json(tmp_path / "execution_contract.json", execution)
    rng = np.random.default_rng(9912)
    cases = []
    for grid in (2, 5):
        records = [{"replicate": i, "seed": spec["seed"], "grid_degrees": grid, "status": "estimated"}
                   for i in range(999)]
        records_path = tmp_path / f"shared_draw_records_{grid}deg.jsonl"
        records_path.write_text("\n".join(json.dumps(r) for r in records) + "\n", encoding="utf-8")
        arrays, summaries, truths, coverage = {}, {}, {}, {}
        for module, dim in spec["module_dimensions"].items():
            point = truth = np.zeros((3, dim, 9))
            draws = rng.normal(size=(999, 3, dim, 9))
            draws[:, 2] = draws[:, 1] - draws[:, 0]
            arrays.update({f"{module}_point": point, f"{module}_draws": draws, f"{module}_truth": truth})
            summaries[module] = summarize_draws(point, draws, blocks)
            truths[module] = truth
            coverage[module] = {"covered": truth.size, "total": truth.size, "rate": 1.}
        arrays_path = tmp_path / f"shared_draws_{grid}deg.npz"
        np.savez_compressed(arrays_path, **arrays)
        cases.append({"scenario": "iid_null", "outer_replicate": 0, "outer_seed": spec["seed"],
                      "grid_degrees": grid, "status": "FINAL_MULTICOORDINATE_OUTER_COMPLETE_NOT_AGGREGATED",
                      "planned_bootstrap_replicates": 999, "estimable_shared_bootstrap_replicates": 999,
                      "shared_draw_records_sha256": digest(records_path), "shared_draws_sha256": digest(arrays_path),
                      "family_slots": _slot_rows(summaries, truths, blocks), "interval_coverage": coverage,
                      "any_false_holm_rejection": False, "false_holm_rejections": 0, "true_holm_rejections": 0,
                      "ecological_fitting_authorized": False})
    report = {"execution_contract": execution, "cases": cases, "empirical_trait_environment_values_read": 0,
              "ecological_models_executed": 0, "ecological_fitting_authorized": False}
    path = tmp_path / "public_report.json"
    _json(path, report)
    return path, report, spec


def test_saved_arrays_recompute_351_coverage_slots_and_all_36_tests(outer):
    path, _, spec = outer
    key, values = verify_outer(path, spec)
    assert key == ("iid_null", 0)
    assert set(values) == {2, 5}
    assert all(v["covered"] == v["total"] == 351 and v["false_family"] == 0 for v in values.values())


@pytest.mark.parametrize("mutation,match", [
    ("unestimated_slot", "probability differs"), ("coverage", "coverage differs"),
    ("missing_grid", "Missing frozen grid"), ("failed_draw", "Incomplete outer"),
    ("wrong_seed", "Case identity differs"),
])
def test_a_complete_label_never_overrides_missing_or_inconsistent_numerics(outer, mutation, match):
    path, report, spec = outer
    case = report["cases"][0]
    if mutation == "unestimated_slot":
        case["family_slots"][0]["candidate_probability"] = None
    elif mutation == "coverage":
        case["interval_coverage"]["orientation"]["covered"] -= 1
    elif mutation == "missing_grid":
        report["cases"].pop()
    elif mutation == "failed_draw":
        case["estimable_shared_bootstrap_replicates"] = 998
    elif mutation == "wrong_seed":
        case["outer_seed"] += 1
    _json(path, report)
    with pytest.raises(ValueError, match=match):
        verify_outer(path, spec)


def test_non_estimable_saved_rank_is_not_hidden_by_retaining_36_slots(outer):
    path, report, spec = outer
    target = path.parent / "shared_draws_2deg.npz"
    with np.load(target) as archive:
        arrays = {k: archive[k] for k in archive.files}
    arrays["orientation_draws"][:] = 0
    np.savez_compressed(target, **arrays)
    report["cases"][0]["shared_draws_sha256"] = digest(target)
    _json(path, report)
    with pytest.raises(ValueError, match="Non-estimable family slot"):
        verify_outer(path, spec)


def test_draw_identity_is_checked_even_when_report_hash_is_updated(outer):
    path, report, spec = outer
    target = path.parent / "shared_draw_records_2deg.jsonl"
    records = target.read_text().splitlines()
    records[-1] = records[0]
    target.write_text("\n".join(records) + "\n", encoding="utf-8")
    report["cases"][0]["shared_draw_records_sha256"] = digest(target)
    _json(path, report)
    with pytest.raises(ValueError, match="reordered or failed"):
        verify_outer(path, spec)


def test_inventory_rejects_missing_duplicate_and_changed_outer_reports(outer):
    path, _, spec = outer
    manifest = [{"scenario": "iid_null", "outer_replicate": 0, "sha256": digest(path)}]
    values, _ = verify_inventory([path], manifest, spec)
    assert len(values) == 1
    with pytest.raises(ValueError, match="every planned outer"):
        verify_inventory([], manifest, spec)
    with pytest.raises(ValueError, match="duplicate outer"):
        verify_inventory([path, path], manifest, spec)
    manifest[0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="report SHA differs"):
        verify_inventory([path], manifest, spec)
