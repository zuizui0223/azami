"""Source-link conservation with synthetic archives, never live API records."""
import csv
import gzip
import importlib
import io
import json
from pathlib import Path
import sqlite3
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
audit = importlib.import_module("analysis.v3.reconcile_sources")


def csv_bytes(rows, merged=False):
    handle = io.StringIO(newline="")
    fields = ["obs_id", "photo_id", "taxon_name"] + (["metadata_chunk_source"] if merged else [])
    writer = csv.DictWriter(handle, fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode("utf-8-sig")


def row(obs, photo, origin=None):
    result = {"obs_id": str(obs), "photo_id": str(photo), "taxon_name": "synthetic"}
    if origin is not None:
        result["metadata_chunk_source"] = str(origin)
    return result


def fixture(tmp_path, compressed=False):
    archives = tmp_path / "archives"
    specs = []
    for origin, rows, seen, start, end, raw in [
        (1, [row(1, 10)], 2, 0, 1, None),
        (2, [row(2, 10), row(3, 30)], 5, 1, 6, [
            {"id": 2, "photos": [{"id": 10}]},
            {"id": 3, "photos": [{"id": 30}, {"id": 10}]},
            {"id": 4, "photos": [{"id": 10}]},
            {"id": 5, "photos": []},
            {"id": 6, "photos": [{"id": None}]},
        ]),
    ]:
        folder = archives / str(origin)
        folder.mkdir(parents=True)
        source = folder / "source.zip"
        metadata = csv_bytes(rows)
        with zipfile.ZipFile(source, "w") as zipped:
            zipped.writestr("photo_metadata.csv", metadata)
            zipped.writestr("collection_provenance.json", json.dumps({
                "resume_state_at_start": {"last_obs_id": start}, "parameters": {"has": "photos"},
                "started_at_utc": "test", "completed_at_utc": "test"}))
            zipped.writestr("checkpoint.json", json.dumps({"observations_seen": seen, "photos_written": len(rows), "last_obs_id": end}))
            if raw is not None:
                raw_bytes = "".join(json.dumps(r) + "\n" for r in raw).encode()
                zipped.writestr("observation_raw.ndjson.gz" if compressed else "observation_raw.ndjson",
                                gzip.compress(raw_bytes, mtime=0) if compressed else raw_bytes)
        specs.append({"artifact_id": origin, "role": "raw_metadata_chunk", "archive_bytes": source.stat().st_size,
                      "archive_sha256": audit.digest(source), "required_member": "photo_metadata.csv",
                      "required_member_sha256": audit.hashlib.sha256(metadata).hexdigest()})
    merged = tmp_path / "merged.csv"
    merged.write_bytes(csv_bytes([row(1, 10, 1), row(3, 30, 2)], merged=True))
    contract = {"source": {"member_sha256": audit.digest(merged), "expected_photo_rows": 2, "expected_observations_with_photos": 2}}
    return archives, merged, {"archives": specs}, contract


def run(tmp_path):
    archives, merged, manifest, contract = fixture(tmp_path)
    return audit.reconcile(archives, merged, tmp_path / "result", manifest, contract)


def test_all_original_rows_and_api_links_preserved(tmp_path):
    report = run(tmp_path)
    counts = report["counts"]
    assert counts["original_metadata_rows"] == 3
    assert counts["merged_metadata_rows"] == 2
    assert counts["retained_unique_observation_photo_links"] == 5
    assert counts["retained_unique_photo_ids"] == 2
    assert counts["source_links_absent_from_merged_metadata"] == 3
    assert report["source_rows_removed_in_this_reconciliation"] == 0
    with sqlite3.connect(tmp_path / "result/source_reconciliation.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM original_photo_rows").fetchone()[0] == 3
        assert db.execute("SELECT COUNT(*) FROM api_observations").fetchone()[0] == 5
        assert db.execute("SELECT COUNT(*) FROM api_photo_links").fetchone()[0] == 5


def test_merge_loss_and_collection_loss_separate(tmp_path):
    report = run(tmp_path)
    assert report["counts"]["original_versions_absent_from_merge"] == 1
    assert report["counts"]["original_unique_links_absent_from_merge"] == 1
    assert report["counts"]["raw_api_unique_links_absent_from_chunk_metadata"] == 2


def test_photo_shared_between_observations_is_not_new_photo(tmp_path):
    report = run(tmp_path)
    assert report["counts"]["photo_ids_with_multiple_observation_links"] == 1
    assert report["counts"]["retained_unique_observation_photo_links"] > report["counts"]["retained_unique_photo_ids"]


def test_no_photo_observations_and_missing_ids_retained(tmp_path):
    report = run(tmp_path)
    assert report["counts"]["retained_unique_observations"] == 6
    assert report["counts"]["source_observations_without_valid_photo_link"] == 2
    assert report["chunks"][1]["raw_observations_without_photos"] == 1
    assert report["chunks"][1]["raw_photo_entries_missing_id"] == 1


def test_missing_raw_archive_is_explicit_gap(tmp_path):
    report = run(tmp_path)
    assert report["raw_api_missing_archives"] == [1]
    assert report["chunks"][0]["raw_api_observations"] is None
    assert report["chunks"][0]["seen_minus_metadata_observations"] == 1
    assert report["image_bytes_verified"] is False
    assert report["ecological_models_executed"] is False


def test_bad_merged_hash_fails_before_output(tmp_path):
    archives, merged, manifest, contract = fixture(tmp_path)
    contract["source"]["member_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="Merged source identity"):
        audit.reconcile(archives, merged, tmp_path / "result", manifest, contract)
    assert not (tmp_path / "result").exists()


def test_unknown_merged_version_fails_closed(tmp_path):
    archives, merged, manifest, contract = fixture(tmp_path)
    changed = row(3, 30, 2)
    changed["taxon_name"] = "changed"
    merged.write_bytes(csv_bytes([row(1, 10, 1), changed], merged=True))
    contract["source"]["member_sha256"] = audit.digest(merged)
    with pytest.raises(ValueError, match="exact original source versions"):
        audit.reconcile(archives, merged, tmp_path / "result", manifest, contract)
    assert (tmp_path / "result/incomplete_run.json").exists()
    assert not (tmp_path / "result/source_reconciliation_report.json").exists()


def test_wrong_archive_stops_reconciliation(tmp_path):
    archives, merged, manifest, contract = fixture(tmp_path)
    manifest["archives"][0]["archive_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="Original archive identity"):
        audit.reconcile(archives, merged, tmp_path / "result", manifest, contract)


def test_existing_outputs_preserved(tmp_path):
    archives, merged, manifest, contract = fixture(tmp_path)
    result = tmp_path / "result"
    result.mkdir()
    marker = result / "keep.txt"
    marker.write_text("keep")
    with pytest.raises(ValueError, match="Output exists"):
        audit.reconcile(archives, merged, result, manifest, contract)
    assert marker.read_text() == "keep"


def test_compressed_api_source_is_enumerated_not_marked_missing(tmp_path):
    archives, merged, manifest, contract = fixture(tmp_path, compressed=True)
    report = audit.reconcile(archives, merged, tmp_path / "result", manifest, contract)
    assert report["raw_api_missing_archives"] == [1]
    assert report["chunks"][1]["raw_api_member"] == "observation_raw.ndjson.gz"
    assert report["counts"]["retained_unique_observation_photo_links"] == 5


def test_repeated_run_has_identical_ledger_and_report(tmp_path):
    archives, merged, manifest, contract = fixture(tmp_path)
    first = audit.reconcile(archives, merged, tmp_path / "first", manifest, contract)
    second = audit.reconcile(archives, merged, tmp_path / "second", manifest, contract)
    assert first == second
    assert audit.digest(tmp_path / "first/source_reconciliation_report.json") == audit.digest(tmp_path / "second/source_reconciliation_report.json")


@pytest.mark.parametrize("text", ["obs_id,photo_id\n1\n", "obs_id,photo_id\n1,2,3\n", "foo,bar\n1,2\n"])
def test_malformed_csv_stops(text):
    with pytest.raises(ValueError):
        list(audit.csv_records(io.StringIO(text)))
