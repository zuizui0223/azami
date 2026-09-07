"""Historical record reconciliation with small synthetic archives."""
import csv
import importlib
import io
import json
from pathlib import Path
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
audit = importlib.import_module("analysis.v3.audit_processing_history")


def csv_text(fields, rows):
    handle = io.StringIO(newline="")
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(fields.split())
    writer.writerows(rows)
    return handle.getvalue()


def inputs(tmp_path, fault=None):
    source = tmp_path / "metadata.csv"
    source.write_text(csv_text("obs_id photo_id", [[1, 11], [1, 12], [2, 21], [3, 31]]))
    counts = {"n_queue_photos": 3, "n_queue_observations": 2, "n_queue_species": 1,
              "n_detected_photos": 1, "n_no_detection_photos": 1, "n_failed_screening_photos": 1,
              "n_observations_with_detected_heads": 1, "n_detected_heads": 1,
              "n_colour_usable": 1, "n_shape_usable": 0, "n_orientation_usable": 1}
    if fault == "saved_count":
        counts["n_detected_heads"] = 2
    members = {
        "exhaustive_photo_queue_merged.csv": csv_text("queue_id photo_id obs_id taxon_name", [["q1", 11, 9 if fault == "source_link" else 1, "A"], ["q2", 12, 1, "A"], ["q3", 21, 2, "A"]]),
        "exhaustive_screening_results.csv": csv_text("queue_id screen_status n_detections", [["q1", "detected", 2 if fault == "head_count" else 1], ["q2", "unexpected" if fault == "state" else "no_detection", 0], ["q3", "error", 0]]),
        "exhaustive_continuous_head_level.csv": csv_text("annotation_unit_id queue_id photo_id obs_id det_index colour_status shape_status orientation_status", [["h1", "q1", 11, 1, 1, "usable", "unusable", "usable"]]),
        "exhaustive_yolo_crop_metadata.csv": csv_text("queue_id det_index photo_id obs_id", [["q1", 1, 12 if fault == "crop_identity" else 11, 1]]),
        "exhaustive_continuous_merge_report.json": json.dumps(counts),
    }
    archive = tmp_path / "source.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        for name, data in members.items():
            zipped.writestr("exhaustive_merged/" + name, data)
        zipped.writestr("postprediction_cohorts/strict_spatial_thinned_observations.csv", csv_text("obs_id taxon_name", [[9 if fault == "cohort" else 1, "A"]]))
    spec = {"artifact_id": 1, "archive_bytes": archive.stat().st_size, "archive_sha256": audit.digest(archive)}
    contract = {"source": {"member_sha256": audit.digest(source), "expected_photo_rows": 4, "expected_observations_with_photos": 3}}
    return archive, source, spec, contract


def test_counts_and_missing_states_remain_distinct(tmp_path):
    archive, source, spec, contract = inputs(tmp_path)
    report = audit.audit_history(archive, source, tmp_path / "out", spec, contract)
    assert report["counts"]["n_queue_photos"] == 3
    assert report["counts"]["source_photos_not_in_historical_queue"] == 1
    assert report["counts"]["queued_observations_with_multiple_photos"] == 1
    assert report["screen_status_counts"] == {"detected": 1, "error": 1, "no_detection": 1}
    assert report["postprediction_views"][0]["observations"] == 1
    assert not any(report["identity_checks"].values())


def test_recount_is_not_new_measurement_or_current_image_availability(tmp_path):
    archive, source, spec, contract = inputs(tmp_path)
    report = audit.audit_history(archive, source, tmp_path / "out", spec, contract)
    assert report["historical_report_counts_matched"]
    assert report["current_image_cache_availability"] == "NOT_ASSESSED"
    assert report["archive_image_members"] == 0
    assert report["new_image_operations"] is False
    assert report["ecological_models_executed"] is False


@pytest.mark.parametrize("fault", ["saved_count", "source_link", "head_count", "state", "crop_identity", "cohort"])
def test_inconsistent_processing_does_not_pass(tmp_path, fault):
    archive, source, spec, contract = inputs(tmp_path, fault)
    with pytest.raises(ValueError):
        audit.audit_history(archive, source, tmp_path / "out", spec, contract)
    assert (tmp_path / "out/incomplete_run.json").exists()
    assert not (tmp_path / "out/historical_processing_report.json").exists()


def test_wrong_archive_fails_before_output_creation(tmp_path):
    archive, source, spec, contract = inputs(tmp_path)
    spec["archive_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="archive identity"):
        audit.audit_history(archive, source, tmp_path / "out", spec, contract)
    assert not (tmp_path / "out").exists()
