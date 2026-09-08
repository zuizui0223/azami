import csv
import json

import pytest

from analysis.v3.verify_original_stream import verify
from analysis.v3.workflow import digest
from test_v3_original_stream_traits import DECISION, offline_run  # noqa: F401
from analysis.v3 import stream_original_traits as stream


def rewrite(out, name, mutate):
    path = out / name
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    mutate(rows)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    report_path = out / "original_stream_report.json"
    report = json.loads(report_path.read_text())
    report["numerical_file_sha256"][name] = digest(path)
    report_path.write_text(json.dumps(report), encoding="utf-8")


def test_recomputes_all27_without_admitting_held_routes(offline_run):
    out, _ = offline_run.invoke()
    result = verify(out, DECISION)
    assert result["raw_endpoint_slots"] == 54
    assert result["bbox_slots"] == 50
    assert result["technical_hold_endpoint_count"] == 13
    assert result["measurement_eligible_heads"]["bract_projection_roughness"] == 0
    assert result["endpoint_status_counts"]["bract_projection_roughness"] == {"usable": 2}
    assert result["endpoint_coverage"]["bract_projection_roughness"]["raw_usable_observations"] == 1
    assert result["endpoint_coverage"]["bract_projection_roughness"]["measurement_eligible_observations"] == 0
    assert result["ecological_fitting_authorized"] is False


def test_rejects_file_tampering(offline_run):
    out, _ = offline_run.invoke()
    with (out / "endpoint_measurements_private.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="hash changed"):
        verify(out, DECISION)


@pytest.mark.parametrize("failure", ["transfer", "no_detection", "head_measurement", "bbox_measurement"])
def test_preserves_failures_without_relabelling_as_usable(offline_run, monkeypatch, failure):
    def fail(*args):
        raise ValueError("deliberate offline failure")
    if failure == "transfer":
        monkeypatch.setattr(stream, "_download", fail)
    elif failure == "no_detection":
        monkeypatch.setattr(stream, "_detector_records", lambda *_: [])
    elif failure == "head_measurement":
        monkeypatch.setattr(stream, "measure_head", fail)
    else:
        original = stream.features.measure
        calls = {"n": 0}
        def fail_shift(*args):
            calls["n"] += 1
            if calls["n"] == 2:
                fail()
            return original(*args)
        monkeypatch.setattr(stream.features, "measure", fail_shift)
    out, _ = offline_run.invoke()
    result = verify(out, DECISION)
    if failure in {"transfer", "no_detection"}:
        assert result["detected_heads"] == result["raw_endpoint_slots"] == 0
    else:
        assert result["raw_endpoint_slots"] == 54
    expected = 1 if failure == "bbox_measurement" else 0
    assert result["measurement_eligible_heads"]["orientation_image_vertical_angle"] == expected


@pytest.mark.parametrize("name,mutate,message", [
    ("endpoint_measurements_private.csv", lambda rows: rows.pop(), "27 raw endpoint slots"),
    ("endpoint_measurements_private.csv", lambda rows: rows.append(rows[0]), "Duplicate row identity"),
    ("bbox_uncertainty_private.csv", lambda rows: rows[0].update(bbox_shift_median_abs_change="999"), "BBox arithmetic"),
    ("bbox_measurements_private.csv", lambda rows: rows.pop(), "BBox condition slots"),
    ("transfer_private.csv", lambda rows: rows.clear(), "Scheduled transfer"),
    ("endpoint_measurements_private.csv", lambda rows: rows[0].update(primary_measurement_eligible="True"), "Measurement eligibility"),
])
def test_detects_arithmetic_and_denominator_errors_even_after_rehash(offline_run, name, mutate, message):
    out, _ = offline_run.invoke()
    rewrite(out, name, mutate)
    with pytest.raises(ValueError, match=message):
        verify(out, DECISION)
