import csv
import importlib
import json
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
grouping = importlib.import_module("analysis.v3.build_dependence_groups")


def fixture(tmp_path):
    root = tmp_path / "workspace"
    root.mkdir()
    with sqlite3.connect(root / "image_workspace.sqlite") as db:
        db.executescript("""
            CREATE TABLE photos(photo_id TEXT PRIMARY KEY);
            CREATE TABLE source_links(obs_id TEXT,photo_id TEXT);
            CREATE TABLE photo_versions(photo_id TEXT,sha256 TEXT);
            CREATE TABLE objects(sha256 TEXT,decoded_rgb_sha256 TEXT);
            CREATE TABLE cache_records(record_id INTEGER,photo_id TEXT,source_path TEXT,pool TEXT);
        """)
        db.executemany("INSERT INTO photos VALUES (?)", [(str(i),) for i in range(1, 8)])
        db.executemany("INSERT INTO source_links VALUES (?,?)", [("a","1"),("b","1"),("b","2"),("c","3"),("d","4"),("e","5"),("f","6"),("g","7")])
        db.executemany("INSERT INTO photo_versions VALUES (?,?)", [("2","bytes1"),("3","bytes1"),("4","bytes2"),("5","bytes3"),("6","bytes4")])
        db.executemany("INSERT INTO objects VALUES (?,?)", [("bytes1","pixels1"),("bytes2","pixels1"),("bytes3","pixels3"),("bytes4","pixels4")])
        db.executemany("INSERT INTO cache_records VALUES (?,?,?,?)", [(1,"2","cache/one.jpg",grouping.DEVELOPMENT_POOL),(2,"4","cache/two.jpg",grouping.DEVELOPMENT_POOL),(3,"3","cache/three.jpg","historical_audit")])
    report = {"workspace_database_sha256": grouping.digest(root/"image_workspace.sqlite"), "counts": {"source_photo_ids":7,"source_observation_ids":7,"source_observation_photo_links":8}}
    (root/"image_workspace_report.json").write_text(json.dumps(report))
    manifest = tmp_path/"training.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["queue_id","split","source_image"])
        writer.writerows([("q1","train","old/images/one.jpg"),("q2","val","old/images/two.jpg")])
    return root, manifest


def run(tmp_path):
    root, manifest = fixture(tmp_path)
    return grouping.build(root, manifest, tmp_path/"out", grouping.digest(manifest))


def test_transitive_photo_bytes_pixels_and_source_retention(tmp_path):
    report = run(tmp_path)
    c = report["counts"]
    assert c["observations_retained"] == c["photos_retained"] == 7
    assert c["source_links_retained"] == 8
    assert c["components"] == 4
    assert c["maximum_observations_in_component"] == 4
    assert not any(report["integrity_checks"].values())
    with sqlite3.connect(tmp_path/"out/dependence_groups.sqlite") as db:
        assert db.execute("SELECT DISTINCT component_id FROM observations WHERE obs_id IN ('a','b','c','d')").fetchall() == [("obs:a",)]


def test_exposure_propagates_across_entire_component(tmp_path):
    c = run(tmp_path)["counts"]
    assert c["components_shared_by_historical_train_and_validation"] == 1
    assert c["cached_objects_in_development_exposed_components"] == 2
    assert c["components_with_multiple_historical_usage_pools"] == 1


def test_union_order_does_not_change_identity():
    one, two = grouping.Components(), grouping.Components()
    for a,b in [("c","b"),("a","d"),("d","c")]:
        one.union(a,b)
    for a,b in [("d","c"),("a","d"),("c","b")]:
        two.union(b,a)
    assert [one.find(k) for k in "abcd"] == [two.find(k) for k in "abcd"] == ["a"]*4


def test_empty_identity_stops():
    with pytest.raises(ValueError, match="Missing"):
        grouping.Components().find("")


def test_stable_fold():
    assert grouping.fold("obs:123") == grouping.fold("obs:123")
    assert 0 <= grouping.fold("obs:123") < 5


def test_changed_source_stops(tmp_path):
    root, manifest = fixture(tmp_path)
    with sqlite3.connect(root/"image_workspace.sqlite") as db:
        db.execute("DELETE FROM source_links WHERE photo_id='7'")
    with pytest.raises(ValueError, match="workspace identity"):
        grouping.build(root, manifest, tmp_path/"out", grouping.digest(manifest))


def test_changed_training_provenance_stops(tmp_path):
    root, manifest = fixture(tmp_path)
    with pytest.raises(ValueError, match="training-manifest identity"):
        grouping.build(root, manifest, tmp_path/"out", "0"*64)


def test_orphan_training_filename_stops(tmp_path):
    root, manifest = fixture(tmp_path)
    manifest.write_text("queue_id,split,source_image\nx,train,unknown.jpg\n")
    with pytest.raises(ValueError, match="not linked"):
        grouping.build(root, manifest, tmp_path/"out", grouping.digest(manifest))
    assert (tmp_path/"out/incomplete_run.json").exists()


def test_repeat_build_byte_identical(tmp_path):
    root, manifest = fixture(tmp_path)
    first = grouping.build(root, manifest, tmp_path/"a", grouping.digest(manifest))
    second = grouping.build(root, manifest, tmp_path/"b", grouping.digest(manifest))
    assert first == second
