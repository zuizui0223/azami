import csv
import importlib
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
recovery = importlib.import_module("analysis.v3.recover_display_composition")


def source(tmp_path):
    path = tmp_path / "heads.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["annotation_unit_id", "photo_id", "obs_id", "colour_status"] + recovery.FIELDS)
        writer.writerows([
            ["h1", "p1", "o1", "usable", .2, 1, 0, 0, 0],
            ["h2", "p1", "o1", "usable", .4, 0, 1, 0, 0],
            ["h3", "p1", "o1", "usable", .6, 0, 0, 1, 0],
            ["h4", "p2", "o1", "usable", .8, 1, 0, 0, 0],
            ["h5", "p3", "o1", "low_colour_quality", .1, 0, 0, 0, 1],
            ["h6", "p4", "o2", "unassessable", "", "", "", "", ""],
        ])
    return path


def run(tmp_path):
    path = source(tmp_path)
    return recovery.recover(path, tmp_path / "out", recovery.digest(path), 6)


def test_equal_photo_weight_not_extra_weight_for_photo_with_many_heads(tmp_path):
    report = run(tmp_path)
    with sqlite3.connect(tmp_path / "out/display_composition.sqlite") as db:
        value, visible = db.execute("SELECT corolla_white_fraction,corolla_visible_fraction FROM observation_values WHERE obs_id='o1'").fetchone()
    assert value == pytest.approx(2 / 3)
    assert visible == pytest.approx(.6)
    assert report["composition_max_abs_sum_error"] < 1e-12


def test_bad_and_missing_head_values_remain_in_source_ledger(tmp_path):
    report = run(tmp_path)
    assert report["counts"]["retained_head_records"] == 6
    assert report["counts"]["source_photos"] == 4
    assert report["counts"]["source_observations"] == 2
    with sqlite3.connect(tmp_path / "out/display_composition.sqlite") as db:
        assert db.execute("SELECT corolla_visible_fraction FROM observation_values WHERE obs_id='o2'").fetchone()[0] is None
        assert db.execute("SELECT n_photos_total FROM observation_values WHERE obs_id='o1'").fetchone()[0] == 3


def test_invalid_composition_is_not_silently_normalized():
    row = dict(zip(recovery.FIELDS, [.5, .4, .4, .4, .4]), colour_status="usable")
    values, visible, composition, conflict = recovery.classify(row)
    assert visible == 1
    assert composition == 0
    assert conflict == 1
    assert sum(values[1:]) == 1.6


def test_visibility_can_be_available_without_complete_composition():
    row = dict(zip(recovery.FIELDS, [.5, "", .4, .3, .3]), colour_status="usable")
    assert recovery.classify(row)[1:] == (1, 0, 1)


def test_bounded_fraction_zero_is_valid_not_missing():
    row = dict(zip(recovery.FIELDS, [0, 0, 0, 0, 1]), colour_status="usable")
    assert recovery.classify(row)[1:] == (1, 1, 0)


def test_recovery_does_not_claim_new_image_or_ecological_execution(tmp_path):
    report = run(tmp_path)
    assert report["new_image_operations"] is False
    assert report["ecological_models_executed"] is False
    assert report["frozen_v2_results_changed"] is False
    assert len(report["endpoints"]) == 5


def test_wrong_source_stops_before_outputs(tmp_path):
    path = source(tmp_path)
    with pytest.raises(ValueError, match="identity"):
        recovery.recover(path, tmp_path / "out", "0" * 64, 6)
    assert not (tmp_path / "out").exists()


def test_wrong_count_fails_closed(tmp_path):
    path = source(tmp_path)
    with pytest.raises(ValueError, match="denominator"):
        recovery.recover(path, tmp_path / "out", recovery.digest(path), 7)
    assert (tmp_path / "out/incomplete_run.json").exists()
    assert not (tmp_path / "out/display_composition_report.json").exists()


def test_repeated_execution_is_identical(tmp_path):
    path = source(tmp_path)
    first = recovery.recover(path, tmp_path / "a", recovery.digest(path), 6)
    second = recovery.recover(path, tmp_path / "b", recovery.digest(path), 6)
    assert first == second
