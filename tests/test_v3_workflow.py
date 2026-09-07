"""Full-source information-retention tests; synthetic records only."""
import csv
import importlib
import json
from pathlib import Path
import sqlite3
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
workflow = importlib.import_module("analysis.v3.workflow")


@pytest.fixture
def contract():
    return json.loads(workflow.CONTRACT.read_text(encoding="utf-8"))


@pytest.fixture
def source(tmp_path, contract):
    rows = [
        ["1", "11", "0", "A", "species", "false", "false", "", "chunk1"],
        ["1", "12", "1", "A", "species", "false", "false", "cc-by", "chunk1"],
        ["2", "21", "0", "Genus", "genus", "unknown", "true", "cc-by-nc", "chunk1"],
        ["3", "31", "0", "B", "species", "true", "false", "cc0", "chunk2"],
        ["3", "31", "0", "B", "species", "true", "false", "cc0", "chunk2"],
        ["4", "11", "0", "C", "species", "true", "", "cc-by", "chunk2"],
    ]
    path = tmp_path / "metadata.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank",
                         "coordinate_usable_for_environment", "captive", "photo_license_code", "metadata_chunk_source"])
        writer.writerows(rows)
    contract["source"].update(member_sha256=workflow.digest(path), expected_photo_rows=6,
                               expected_observations_with_photos=4)
    return path


def test_full_source_contract_preserves_data(contract):
    workflow.validate_contract(contract)
    assert contract["source"]["expected_photo_rows"] == 1122854
    assert contract["source"]["expected_observations_with_photos"] == 665115


@pytest.mark.parametrize("field", ["drop_rows_for_qc", "thin_at_acquisition", "first_photo_only",
                                   "native_only_at_acquisition", "require_coordinates_for_image_measurement",
                                   "delete_detector_negative_photos"])
def test_cannot_move_analysis_filters_to_acquisition(contract, field):
    contract["preservation"][field] = True
    with pytest.raises(ValueError, match="acquisition universe"):
        workflow.validate_contract(contract)


@pytest.mark.parametrize("field", ["retain_every_source_row", "retain_all_photos_per_observation"])
def test_cannot_silently_drop_records(contract, field):
    contract["preservation"][field] = False
    with pytest.raises(ValueError):
        workflow.validate_contract(contract)


@pytest.mark.parametrize("field", ["no_required_result_direction", "no_requirement_that_candidates_survive"])
def test_no_preselected_conclusion(contract, field):
    contract["reporting"][field] = False
    with pytest.raises(ValueError):
        workflow.validate_contract(contract)


def test_all_photos_and_restricted_rows_retained(source, tmp_path, contract):
    out = tmp_path / "run"
    report = workflow.inventory_records(source, out, contract)
    counts = report["counts"]
    assert counts["retained_photo_record_rows"] == 6
    assert counts["rows_removed_in_v3_inventory"] == 0
    assert counts["observations_with_photos"] == 4
    assert counts["observations_with_multiple_photos"] == 1
    assert counts["unique_photo_ids"] == 4
    assert counts["duplicate_observation_photo_rows"] == 1
    assert counts["photo_ids_linked_to_multiple_observations"] == 1
    assert report["photo_record_counts_by_status"]["coordinate_status"] == {"false": 2, "true": 3, "unknown": 1}
    assert report["photo_record_counts_by_status"]["captive_status"]["true"] == 1
    with sqlite3.connect(out / "source_ledger.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM photo_records WHERE obs_id='1'").fetchone()[0] == 2
        assert db.execute("SELECT COUNT(*) FROM photo_records WHERE taxon_rank='genus'").fetchone()[0] == 1


def test_legacy_subset_is_annotation_not_filter(source, tmp_path, contract):
    report = workflow.inventory_records(source, tmp_path / "run", contract,
                                         {"1": ("A", "native"), "2": ("Genus", "introduced"), "99": ("Z", "native")})
    assert report["counts"]["retained_photo_record_rows"] == 6
    assert report["legacy_v2_comparison"]["overlapping_observations"] == 2
    assert report["legacy_v2_comparison"]["missing_from_photo_snapshot"] == 1
    assert report["photo_record_counts_by_status"]["native_status_v2"]["not_assessed_outside_v2"] == 3


def test_changed_taxon_does_not_inherit_native_status(source, tmp_path, contract):
    report = workflow.inventory_records(source, tmp_path / "run", contract, {"1": ("Wrong", "native")})
    assert report["photo_record_counts_by_status"]["native_status_v2"]["not_assessed_taxon_changed"] == 2


def test_missing_legacy_input_is_unavailable_not_zero(source, tmp_path, contract):
    report = workflow.inventory_records(source, tmp_path / "run", contract)
    assert report["legacy_v2_comparison"]["status"] == "NOT_SUPPLIED"
    assert report["legacy_v2_comparison"]["overlapping_observations"] is None


def test_inventory_does_not_certify_image_availability(source, tmp_path, contract):
    report = workflow.inventory_records(source, tmp_path / "run", contract)
    assert report["image_download_completeness"] == "NOT_ASSESSED_FROM_METADATA"
    assert report["detector_negative_count"] is None
    assert not report["image_operations_performed"]
    assert not report["ecological_models_executed"]
    assert report["stages"]["recover"] == "NOT_EXECUTED"


def test_rerun_is_byte_identical_and_source_unchanged(source, tmp_path, contract):
    before = workflow.digest(source)
    a = workflow.inventory_records(source, tmp_path / "a", contract)
    b = workflow.inventory_records(source, tmp_path / "b", contract)
    assert a == b
    assert workflow.digest(source) == before
    assert workflow.digest(tmp_path / "a/source_inventory_report.json") == workflow.digest(tmp_path / "b/source_inventory_report.json")


def test_bad_hash_stops_before_output_creation(source, tmp_path, contract):
    contract["source"]["member_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA-256"):
        workflow.inventory_records(source, tmp_path / "out", contract)
    assert not (tmp_path / "out").exists()


def test_wrong_denominator_leaves_explicit_incomplete_run(source, tmp_path, contract):
    contract["source"]["expected_photo_rows"] = 5
    out = tmp_path / "run"
    with pytest.raises(ValueError, match="counts differ"):
        workflow.inventory_records(source, out, contract)
    assert not (out / "source_inventory_report.json").exists()
    assert json.loads((out / "incomplete_run.json").read_text())["status"] == "INCOMPLETE_DO_NOT_USE"


def test_missing_identifier_is_retained_not_dropped(source, tmp_path, contract):
    text = source.read_text().replace("4,11,0,C", ",11,0,C")
    source.write_text(text)
    contract["source"]["member_sha256"] = workflow.digest(source)
    contract["source"]["expected_observations_with_photos"] = 3
    report = workflow.inventory_records(source, tmp_path / "run", contract)
    assert report["counts"]["rows_with_missing_identifiers"] == 1
    assert report["counts"]["retained_photo_record_rows"] == 6


def test_existing_run_cannot_be_overwritten(source, tmp_path, contract):
    with pytest.raises(ValueError, match="already exists"):
        workflow.inventory_records(source, tmp_path, contract)


def test_frozen_repository_outputs_are_protected(source, contract):
    with pytest.raises(ValueError, match="external directory"):
        workflow.inventory_records(source, ROOT / "analysis_outputs/do-not-create", contract)


def test_plan_command_needs_no_data():
    run = subprocess.run([sys.executable, "-m", "analysis.v3.workflow", "plan"], cwd=ROOT,
                         text=True, capture_output=True, check=True)
    report = json.loads(run.stdout)
    assert report["status"] == "DESIGN_CHECKED_NOT_EXECUTED"
    assert report["source"]["expected_photo_rows"] > 1000000
