import csv
import importlib
from pathlib import Path
import sqlite3
import sys

from PIL import Image, PngImagePlugin
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
cache = importlib.import_module("analysis.v3.build_image_workspace")


def fixture(tmp_path):
    metadata = tmp_path / "source.csv"
    with metadata.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["photo_id", "medium_image_url", "large_image_url", "photo_license_code"])
        writer.writerows([[str(i), "https://example.invalid/m.jpg", "https://example.invalid/l.jpg", ""] for i in range(1, 4)])
    links = tmp_path / "links.sqlite"
    with sqlite3.connect(links) as db:
        db.execute("CREATE TABLE source_links (obs_id TEXT, photo_id TEXT)")
        db.executemany("INSERT INTO source_links VALUES (?,?)", [("o1", "1"), ("o2", "1"), ("o3", "2"), ("o4", "3")])
    one = tmp_path / "one.png"
    two = tmp_path / "two.png"
    image = Image.new("RGB", (20, 30), "purple")
    image.save(one)
    info = PngImagePlugin.PngInfo()
    info.add_text("description", "same pixels, different encoding")
    image.save(two, pnginfo=info)
    declarations = [dict(photo_id=p, path=str(f), pool=pool, expected_sha256=cache.digest(f), source_manifest="synthetic.csv", source_row=i)
                    for i, (p, f, pool) in enumerate([("1", one, "development"), ("2", one, "audit"), ("1", two, "audit")], 1)]
    contract = {"source": {"member_sha256": cache.digest(metadata), "expected_photo_rows": 3}}
    return metadata, links, declarations, contract


def run(tmp_path, modifier=None):
    metadata, links, declarations, contract = fixture(tmp_path)
    if modifier:
        modifier(declarations)
    return cache.build_workspace(metadata, links, tmp_path / "out", declarations, [], contract, cache.digest(links))


def test_all_photo_ids_and_shared_observation_links_retained(tmp_path):
    report = run(tmp_path)
    c = report["counts"]
    assert c["source_photo_ids"] == 3
    assert c["source_observation_photo_links"] == 4
    assert c["source_observation_ids"] == 4
    assert c["photo_ids_without_verified_local_image"] == 1
    assert c["source_photo_rows_deleted"] == 0


def test_content_and_photo_versions_are_distinct(tmp_path):
    c = run(tmp_path)["counts"]
    assert c["verified_image_objects"] == 2
    assert c["retained_photo_content_versions"] == 3
    assert c["objects_linked_to_multiple_photo_ids"] == 1
    assert c["photos_with_multiple_cached_versions"] == 1
    assert c["decoded_pixel_groups_with_multiple_encodings"] == 1


def test_original_bytes_and_pools_preserved(tmp_path):
    run(tmp_path)
    with sqlite3.connect(tmp_path / "out/image_workspace.sqlite") as db:
        for sha, path in db.execute("SELECT sha256,relative_path FROM objects"):
            assert cache.digest(cache.readable(tmp_path / "out" / path)) == sha
        assert db.execute("SELECT COUNT(*) FROM cache_records").fetchone()[0] == 3
        assert db.execute("SELECT COUNT(*) FROM image_jobs WHERE status='pending_detection'").fetchone()[0] == 2


def test_changed_source_is_not_imported_as_valid(tmp_path):
    report = run(tmp_path, lambda rows: rows[0].update(expected_sha256="0" * 64))
    assert report["cache_record_states"]["source_hash_mismatch"] == 1
    assert report["counts"]["source_photo_ids"] == 3


def test_missing_cache_is_not_detection_negative(tmp_path):
    report = run(tmp_path, lambda rows: rows[0].update(path=str(tmp_path / "absent.png")))
    assert report["cache_record_states"]["missing_local_file"] == 1
    assert report["detector_executed"] is False
    assert report["new_source_photo_downloads"] == 0


def test_orphan_photo_stops_run(tmp_path):
    with pytest.raises(ValueError, match="outside"):
        run(tmp_path, lambda rows: rows[0].update(photo_id="unknown"))
    assert (tmp_path / "out/incomplete_run.json").exists()
    assert not (tmp_path / "out/image_workspace_report.json").exists()


def test_bad_source_hash_stops_before_outputs(tmp_path):
    metadata, links, rows, contract = fixture(tmp_path)
    with pytest.raises(ValueError, match="identity"):
        cache.build_workspace(metadata, links, tmp_path / "out", rows, [], contract, "0" * 64)
    assert not (tmp_path / "out").exists()


def test_corrupt_image_gets_explicit_state(tmp_path):
    source = tmp_path / "bad.jpg"
    source.write_bytes(b"not an image")
    assert cache.inspect_and_copy(source, tmp_path / "out", None)["status"] == "decode_failed"


def test_exif_transform_is_recorded_and_original_bytes_unchanged(tmp_path):
    source = tmp_path / "rotated.jpg"
    exif = Image.Exif()
    exif[274] = 6
    Image.new("RGB", (20, 30), "green").save(source, exif=exif)
    before = cache.digest(source)
    record = cache.inspect_and_copy(source, tmp_path / "out", before)
    assert (record["width"], record["height"]) == (30, 20)
    assert record["exif_orientation"] == 6
    assert cache.digest(source) == before


@pytest.mark.parametrize("name", ["../photo.jpg", "a/b.jpg", "C:\\a.jpg", "..", ""])
def test_manifest_path_traversal_rejected(name):
    with pytest.raises(ValueError):
        cache.safe_filename(name)
