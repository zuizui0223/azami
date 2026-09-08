import csv
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import zipfile

import pytest

from analysis.v3 import reconciled_photo_schedule as schedule
from analysis.v3.enriched_source_cohort import COHORT_FIELDS, EXTRA_FIELDS
from analysis.v3.workflow import canonical_digest, digest


def fixture(tmp_path):
    source = tmp_path / "source.sqlite"
    metadata_rows = []
    for obs, photo, license_code in [("1", "10", "cc-by"), ("2", "20", "cc-by"),
                                      ("3", "30", ""), ("99", "10", "cc0")]:
        metadata_rows.append({"obs_id": obs, "photo_id": photo, "photo_license_code": license_code,
            "photo_attribution": "Synthetic", "medium_image_url": f"https://static.inaturalist.org/photos/{photo}/medium.jpg",
            "large_image_url": f"https://static.inaturalist.org/photos/{photo}/large.jpg", "latitude": "SECRET_LOCATION"})
    api = {"id": 99, "photos": [{"id": 20, "license_code": "cc-by", "attribution": "Synthetic",
                                 "url": "https://static.inaturalist.org/photos/20/square.jpg"}],
           "private_latitude": "NEVER_COPY"}
    raw = (json.dumps(api) + "\n").encode()
    with sqlite3.connect(source) as db:
        db.executescript("""
            CREATE TABLE source_links(obs_id TEXT,photo_id TEXT,PRIMARY KEY(obs_id,photo_id));
            CREATE TABLE original_photo_rows(origin TEXT,source_row INTEGER,obs_id TEXT,photo_id TEXT,payload_sha256 TEXT);
            CREATE TABLE api_observations(origin TEXT,source_line INTEGER,obs_id TEXT,photo_array_count INTEGER,raw_line_sha256 TEXT);
            CREATE TABLE api_photo_links(origin TEXT,source_line INTEGER,photo_index INTEGER,obs_id TEXT,photo_id TEXT);
        """)
        db.executemany("INSERT INTO source_links VALUES (?,?)", [(r["obs_id"], r["photo_id"]) for r in metadata_rows] + [("99", "20")])
        db.executemany("INSERT INTO original_photo_rows VALUES (?,?,?,?,?)", [
            ("1", i, r["obs_id"], r["photo_id"], canonical_digest(r)) for i, r in enumerate(metadata_rows, 1)])
        db.execute("INSERT INTO api_observations VALUES ('1',1,'99',1,?)", (hashlib.sha256(raw).hexdigest(),))
        db.execute("INSERT INTO api_photo_links VALUES ('1',1,1,'99','20')")
    folder = tmp_path / "archives/1"
    folder.mkdir(parents=True)
    archive = folder / "source.zip"
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(metadata_rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(metadata_rows)
    csv_bytes = buffer.getvalue().encode()
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("photo_metadata.csv", csv_bytes)
        zipped.writestr("observation_raw.ndjson", raw)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"trait_files_in_scope": 0, "archives": [{"artifact_id": 1,
        "role": "raw_metadata_chunk", "archive_bytes": archive.stat().st_size, "archive_sha256": digest(archive),
        "required_member": "photo_metadata.csv", "required_member_sha256": hashlib.sha256(csv_bytes).hexdigest()}]}))
    enriched = tmp_path / "enriched.csv"
    rows = []
    for obs, component in [("1", "component-a"), ("2", "component-a"), ("3", "component-b")]:
        row = {key: "" for key in COHORT_FIELDS + EXTRA_FIELDS}
        row.update(obs_id=obs, dependence_component_id=component, native_range_status="native",
                   source_taxon_rank="species", reconciled_source_photo_count="1")
        rows.append(row)
    with enriched.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COHORT_FIELDS + EXTRA_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    receipt = {"output_csv_sha256": digest(enriched), "input_sha256": {"reconciliation": digest(source)},
               "cohort_rows": 3, "dependence": {"cohort_known_components": 2}}
    (tmp_path / "enriched_source_cohort_report.json").write_text(json.dumps(receipt))
    return dict(enriched=enriched, reconciliation=source, archives=tmp_path / "archives", out=tmp_path / "output",
        expected_enriched_sha256=digest(enriched), expected_reconciliation_sha256=digest(source), manifest_path=manifest)


def test_source_links_versions_blocked_states_and_blind_worker_database(tmp_path):
    kwargs = fixture(tmp_path)
    report = schedule.build(**kwargs)
    assert report["counts"] == {"native_observations": 3, "native_links": 3, "known_photo_links": 5,
                                 "photo_versions": 5, "photo_jobs": 3}
    assert report["photo_states"] == {"license_conflict": 1, "license_unavailable": 1, "request_candidate_not_authorized": 1}
    assert report["native_observations_with_request_candidates"] == 1
    assert report["production_image_execution_authorized"] is False
    path = kwargs["out"] / "reconciled_photo_schedule_private.sqlite"
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(DISTINCT shard_index) FROM photo_jobs WHERE component_id='component-a'").fetchone()[0] == 1
        dump = "\n".join(db.iterdump())
        assert "SECRET_LOCATION" not in dump and "NEVER_COPY" not in dump
        assert '"latitude"' not in dump and '"source_taxon_rank"' not in dump
        assert db.execute("SELECT original_url FROM photo_jobs WHERE photo_id='20'").fetchone()[0].endswith("/original.jpg")
    assert digest(path) == report["output_sqlite_sha256"]


@pytest.mark.parametrize("url", [
    "https://example.test/photos/10/large.jpg", "file:///photos/10/large.jpg",
    "https://static.inaturalist.org/photos/99/large.jpg", "https://user:pass@static.inaturalist.org/photos/10/large.jpg",
    "https://static.inaturalist.org:443/photos/10/large.jpg", "https://static.inaturalist.org/photos/10/large.jpg?x=1",
])
def test_bad_or_mismatched_urls_are_never_request_candidates(url):
    assert schedule.original_identity("10", {"urls": [url]})[0] == "url_invalid_or_identity_conflict"


def test_url_sizes_agree_but_missing_or_different_sources_are_not_imputed():
    urls = [f"https://static.inaturalist.org/photos/10/{size}.jpg" for size in ("square", "medium", "large")]
    assert schedule.original_identity("10", {"urls": urls}) == ("valid", "https://static.inaturalist.org/photos/10/original.jpg")
    assert schedule.original_identity("10", {"urls": [""]}) == ("url_unavailable", "")
    urls.append("https://inaturalist-open-data.s3.amazonaws.com/photos/10/large.jpg")
    assert schedule.original_identity("10", {"urls": urls}) == ("url_conflict", "")


def test_input_identity_change_stops_before_output(tmp_path):
    kwargs = fixture(tmp_path)
    kwargs["expected_enriched_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="input identity"):
        schedule.build(**kwargs)
    assert not kwargs["out"].exists()


def test_archive_damage_is_not_silently_skipped(tmp_path):
    kwargs = fixture(tmp_path)
    archive = kwargs["archives"] / "1/source.zip"
    archive.write_bytes(archive.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="Archive byte count"):
        schedule.build(**kwargs)
    assert not (kwargs["out"] / "reconciled_photo_schedule_report.json").exists()
    assert json.loads((kwargs["out"] / "incomplete_run.json").read_text())["status"] == "INCOMPLETE_DO_NOT_USE"


def test_missing_reconciled_link_and_changed_source_locator_fail(tmp_path):
    kwargs = fixture(tmp_path)
    with sqlite3.connect(kwargs["reconciliation"]) as db:
        db.execute("UPDATE original_photo_rows SET payload_sha256='wrong' WHERE source_row=1")
    kwargs["expected_reconciliation_sha256"] = digest(kwargs["reconciliation"])
    receipt_path = kwargs["enriched"].with_name("enriched_source_cohort_report.json")
    receipt = json.loads(receipt_path.read_text())
    receipt["input_sha256"]["reconciliation"] = kwargs["expected_reconciliation_sha256"]
    receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="locator/hash"):
        schedule.build(**kwargs)


def test_existing_schedule_is_never_overwritten(tmp_path):
    kwargs = fixture(tmp_path)
    schedule.build(**kwargs)
    original = digest(kwargs["out"] / "reconciled_photo_schedule_private.sqlite")
    with pytest.raises(ValueError, match="Fresh output"):
        schedule.build(**kwargs)
    assert original == digest(kwargs["out"] / "reconciled_photo_schedule_private.sqlite")


def test_worker_packet_keeps_whole_components_and_all_blocked_links(tmp_path):
    from analysis.v3.reconciled_stream_input import pilot_input
    kwargs = fixture(tmp_path)
    report = schedule.build(**kwargs)
    path = kwargs["out"] / "reconciled_photo_schedule_private.sqlite"
    packet = pilot_input(path, report["output_sqlite_sha256"], 3)
    assert set(packet["selected"]) == {"1", "2", "3"}
    assert packet["selection_scores"]["1"] == packet["selection_scores"]["2"]
    assert len(packet["links"]) == 3 and len(packet["queue"]) == 1
    assert packet["queue"][0]["known_metadata_obs_ids"] == ["2", "99"]
    assert {v["kind"] for v in packet["queue"][0]["source_versions"]} == {"metadata", "api"}
    assert "SECRET_LOCATION" not in json.dumps(packet) and "NEVER_COPY" not in json.dumps(packet)
    assert digest(path) == report["output_sqlite_sha256"]
    contract = json.loads((schedule.ROOT / "analysis/v3/reconciled_photo_schedule_contract.json").read_text())
    first = min((schedule.component_score(c, contract["ordering_salt"]), c) for c in ("component-a", "component-b"))[1]
    if first == "component-a":
        with pytest.raises(ValueError, match="First complete component"):
            pilot_input(path, report["output_sqlite_sha256"], 1)
    else:
        bounded = pilot_input(path, report["output_sqlite_sha256"], 2)
        assert bounded["selected"] == ["3"]  # Never fill the remaining slot by splitting a component.
        assert bounded["queue"] == []  # No replacement of blocked photos with available ones.


def test_worker_packet_requires_exact_pin_and_valid_execution_contract(tmp_path):
    from analysis.v3.reconciled_stream_input import pilot_input
    kwargs = fixture(tmp_path)
    report = schedule.build(**kwargs)
    path = kwargs["out"] / "reconciled_photo_schedule_private.sqlite"
    with pytest.raises(ValueError, match="identity mismatch"):
        pilot_input(path, "0" * 64, 3)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE execution SET value_json='true' WHERE key='production_image_execution_authorized'")
    with pytest.raises(ValueError, match="execution contract"):
        pilot_input(path, digest(path), 3)
