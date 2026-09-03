from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "ch1_global" / "v2" / "ANALYSIS_MANIFEST.tsv"
CATALOG = ROOT / "reproducibility" / "actions_artifact_catalog.json"
RECOVERY_INVENTORY = ROOT / "reproducibility" / "recovery_inventory.json"
PYPROJECT = ROOT / "pyproject.toml"
FROZEN_REQUIREMENTS = ROOT / "reproducibility" / "frozen-numerical-rebuild-requirements.txt"
ATLAS_VALIDATION = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27" / "v2_full27_environment_validation.json"
SAMPLING_VALIDATION = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27" / "sampling" / "v2_full27_sampling_composition_validation.json"
SENSITIVITY_VALIDATION = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27" / "v2_full27_sensitivities_validation.json"
FIGURE_REPORT = ROOT / "reproducibility" / "figures" / "figure_build_report.json"


def load_rows() -> list[dict[str, str]]:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_manifest_has_only_live_repository_paths() -> None:
    rows = load_rows()
    assert rows
    for row in rows:
        relative = row["file"]
        assert relative, row
        assert (ROOT / relative).exists(), relative


def test_manifest_contains_current_reproducibility_entrypoints() -> None:
    rows = load_rows()
    by_stage = {row["stage"]: row for row in rows}
    required = {
        "canonical_rebuild_runner": "analysis/rebuild_frozen_analysis.py",
        "broad_region_lookup": "analysis/rebuild_frozen_broad_region_lookup.py",
        "broad_region_assignment_builder": "analysis/build_naturalearth_broad_region_lookup.py",
        "native_status_rebuild": "analysis/rebuild_frozen_native_status.py",
        "native_status_contract": "analysis/ch1/native_range_sensitivity_contract.json",
        "frozen_python_environment": "reproducibility/frozen-numerical-rebuild-requirements.txt",
        "full27_atlas": "analysis/run_geb_v2_full27_environment_atlas.py",
        "full27_validation": "analysis/validate_geb_v2_full27_environment_atlas.py",
        "figure_rebuild": ".github/workflows/rebuild-frozen-figures.yml",
        "integrity_ci": ".github/workflows/reproducibility-integrity.yml",
        "source_artifact_catalog": "reproducibility/actions_artifact_catalog.json",
        "frozen_outputs": "analysis_outputs/README.md",
    }
    for stage, path in required.items():
        assert by_stage[stage]["file"] == path


def test_removed_submission_workflows_are_not_marked_active() -> None:
    text = MANIFEST.read_text(encoding="utf-8")
    stale = (
        ".github/workflows/ch1-lightweight-spde-inla.yml",
        ".github/workflows/ch1-submission-ci.yml",
        ".github/workflows/ch1-build-submission-bundle.yml",
    )
    for path in stale:
        assert path not in text


def test_full_extra_covers_recovered_pipeline_dependencies() -> None:
    payload = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    full = payload["project"]["optional-dependencies"]["full"]
    names = {requirement.split("<", 1)[0].split(">", 1)[0].split("=", 1)[0] for requirement in full}
    assert {"rasterio", "pyproj", "geopandas", "pyogrio", "ultralytics"}.issubset(names)


def test_frozen_numerical_environment_matches_recovered_artifact_versions() -> None:
    lines = {
        line.strip()
        for line in FROZEN_REQUIREMENTS.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    required = {
        "numpy==2.4.6",
        "pandas==3.0.5",
        "scipy==1.17.1",
        "statsmodels==0.14.6",
        "requests==2.34.2",
        "rasterio==1.4.4",
        "pyproj==3.7.2",
        "shapely==2.1.2",
        "geopandas==1.1.4",
        "pyogrio==0.13.0",
    }
    assert required.issubset(lines)


def test_frozen_validation_reports_are_still_complete() -> None:
    expected = {
        ATLAS_VALIDATION: ("PASS", 24, 24),
        SAMPLING_VALIDATION: ("PASS", 7, 7),
        SENSITIVITY_VALIDATION: ("PASS", 11, 11),
    }
    for path, (status, n_checks, n_passed) in expected.items():
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["status"] == status
        assert payload["n_checks"] == n_checks
        assert payload["n_passed"] == n_passed

    figure = json.loads(FIGURE_REPORT.read_text(encoding="utf-8"))
    assert figure["status"] == "ok"
    assert len(figure["main_figure_stems"]) == 5
    assert len(figure["supporting_figure_stems"]) == 7


def test_recovery_inventory_closes_active_v2_and_preserves_history() -> None:
    payload = json.loads(RECOVERY_INVENTORY.read_text(encoding="utf-8"))
    assert payload["active_v2"]["status"] == "complete_for_frozen_chapter1_v2_reanalysis"
    assert payload["active_v2"]["baseline_ci"]["reproducibility_integrity_conclusion"] == "success"
    assert payload["active_v2"]["baseline_ci"]["frozen_figure_conclusion"] == "success"

    completion = payload["completion"]
    for key in (
        "active_v2_code",
        "active_v2_frozen_outputs",
        "active_v2_direct_inputs",
        "active_v2_figures",
        "active_v2_ci",
        "precleanup_scientific_history",
    ):
        assert completion[key] == "complete"
    assert completion["full_upstream_binary_durability"] == "tracked_external_gate_issue_85"

    recovery = payload["historical_recovery"]
    assert recovery["immutable_tag"] == "azami-ch1-v2-2026-08-27"
    assert recovery["tagged_commit"] == "584e4a9863ab0c133a050f15b4bb862730db7faf"
    assert recovery["browseable_archive_branch"] == "archive/ch1-precleanup-20260827"
    assert recovery["status"] == "complete"

    direct = payload["direct_reanalysis_inputs"]
    assert direct["continuous_source"]["durability"] == "owner_side_archive_verified"
    assert direct["multilevel_process_source"]["durability"] == "owner_side_archive_verified"
    assert direct["native_status"]["durability"] == "git_retained_and_rebuildable"
    assert direct["broad_region"]["durability"] == "deterministically_rebuildable_fail_closed"
    assert payload["ephemeral_upstream_artifacts"]["tracking_issue"] == 85
    assert payload["ephemeral_upstream_artifacts"]["not_required_for_current_frozen_v2_numerical_reanalysis"] is True


def test_artifact_catalog_retains_reconstruction_checkpoints() -> None:
    payload = json.loads(CATALOG.read_text(encoding="utf-8"))
    assert payload["catalog_version"] >= 4
    assert payload["external_archive_tracking_issue"] == 85
    by_id = {int(row["artifact_id"]): row for row in payload["artifacts"]}
    required_ids = {
        8076736948,
        8269246732,
        8227254443,
        8983877726,
        9612943217,
        9632715852,
        8521924881,
        8521925441,
        8521926057,
    }
    assert required_ids.issubset(by_id)
    assert by_id[8076736948]["model_weight_sha256"] == (
        "4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed"
    )
    for artifact_id in (9612943217, 9632715852):
        row = by_id[artifact_id]
        assert row["github_digest"].removeprefix("sha256:") == row["local_archive_sha256"]
    historical = by_id[8227254443]
    assert historical["github_digest"] == (
        "sha256:499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc"
    )
    assert historical["tree_files"]["historical_trees/gbotb_lcvp_scenario2_randomized.trees"] == (
        "sha256:82655f79297e44a6630a599d8b0a1dc6f85e792812b8b01f2234e567a478e3af"
    )
    spatial = by_id[8983877726]
    assert spatial["inner_files"]["spatial_regions/broad_region_lookup.csv"] == (
        "sha256:085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff"
    )
    continuous = by_id[9612943217]
    assert continuous["inner_files"]["environment/strict_spatial_chelsa.csv"] == (
        "sha256:2172e3570f684770d0f919ecd81265c8460574e287bc4fb057db4f719cab7bb0"
    )
    qc_companion = by_id[8521925441]
    assert qc_companion["github_digest"] == (
        "sha256:4c39fe51b155d0eb09038e470088ababce4367f740a476dc38127d801d2e0087"
    )
    native = payload["git_frozen_inputs"]["sampling_native_status"]
    assert native["immutable_tag"] == "azami-ch1-v2-2026-08-27"
    assert native["git_blob_sha"] == "b98af47482fd86b1353546573492519659cda848"
    assert native["sha256_recorded_by_frozen_sampling_report"] == (
        "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a"
    )
    assert native["source_cohort_sha256"] == (
        "2172e3570f684770d0f919ecd81265c8460574e287bc4fb057db4f719cab7bb0"
    )
    assert payload["historical_recovery"]["immutable_tag"] == "azami-ch1-v2-2026-08-27"
    assert payload["historical_recovery"]["archive_branch"] == "archive/ch1-precleanup-20260827"
