"""Reconcile original photo rows, archived API links and the historical merge.

Read pinned ZIP archives, not the live API. Keep source records and a unique
observation-photo link view separately; a repeated photo is not a new replicate.
Raw IDs and record locators stay local. Only the aggregate report is public.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import zipfile

from .workflow import CONTRACT, ROOT, canonical_digest, digest, text_digest


def csv_records(handle):
    reader = csv.DictReader(handle)
    if not {"obs_id", "photo_id"}.issubset(reader.fieldnames or []):
        raise ValueError("Source requires observation and photo identifiers")
    for number, row in enumerate(reader, 1):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed CSV row {number}")
        yield number, row


def id_text(value):
    if value is None or isinstance(value, bool):
        return ""
    return str(value)


def ingest_archive(db, archive: Path, spec: dict) -> dict:
    origin = str(spec["artifact_id"])
    if archive.stat().st_size != spec["archive_bytes"] or digest(archive) != spec["archive_sha256"]:
        raise ValueError("Original archive identity mismatch: " + origin)
    report = {"artifact_id": spec["artifact_id"], "archive_sha256": spec["archive_sha256"]}
    with zipfile.ZipFile(archive) as zipped:
        names = zipped.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate archive members")
        member = spec["required_member"]
        with zipped.open(member) as handle:
            actual = hashlib.file_digest(handle, "sha256").hexdigest()
        if actual != spec["required_member_sha256"]:
            raise ValueError("Original photo metadata identity mismatch")
        with zipped.open(member) as binary, io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as handle:
            rows = 0
            batch = []
            for rows, row in csv_records(handle):
                batch.append((origin, rows, row["obs_id"], row["photo_id"], canonical_digest(row)))
                if len(batch) == 20000:
                    db.executemany("INSERT INTO original_photo_rows VALUES (?,?,?,?,?)", batch)
                    batch.clear()
            db.executemany("INSERT INTO original_photo_rows VALUES (?,?,?,?,?)", batch)
        report["metadata_rows"] = rows
        report["metadata_sha256"] = actual
        provenance = json.loads(zipped.read("collection_provenance.json"))
        checkpoint = json.loads(zipped.read("checkpoint.json"))
        report["collection"] = {
            "requested_start_obs_id": provenance.get("chunk_requested_start_obs_id", provenance["resume_state_at_start"]["last_obs_id"]),
            "last_obs_id": checkpoint["last_obs_id"],
            "observations_seen": checkpoint["observations_seen"],
            "photos_written": checkpoint["photos_written"],
            "started_at_utc": provenance["started_at_utc"],
            "completed_at_utc": provenance["completed_at_utc"],
            "query_parameters": provenance["parameters"],
        }
        if checkpoint["photos_written"] != rows:
            raise ValueError("Collection photo count does not reconcile")
        raw_members = [n for n in names if n in ("observation_raw.ndjson", "observation_raw.ndjson.gz")]
        if len(raw_members) > 1:
            raise ValueError("Ambiguous raw API source versions")
        raw_name = raw_members[0] if raw_members else None
        report["raw_api_member"] = raw_name
        if raw_name is None:
            report["raw_api_status"] = "NOT_PRESENT_IN_ARCHIVE"
            report["raw_api_observations"] = None
            return report
        raw_hash = hashlib.sha256()
        observations = []
        links = []
        absent_id = no_photos = photo_id_missing = observation_count = link_count = 0
        with zipped.open(raw_name) as binary:
            handle = gzip.GzipFile(fileobj=binary) if raw_name.endswith(".gz") else binary
            for line, payload in enumerate(handle, 1):
                raw_hash.update(payload)
                record = json.loads(payload)
                if not isinstance(record, dict):
                    raise ValueError("Invalid raw API record")
                obs = id_text(record.get("id"))
                photos = record.get("photos") or []
                if not isinstance(photos, list) or any(not isinstance(p, dict) for p in photos):
                    raise ValueError("Invalid raw API photo array")
                observations.append((origin, line, obs, len(photos), hashlib.sha256(payload).hexdigest()))
                observation_count += 1
                absent_id += not bool(obs)
                no_photos += not bool(photos)
                for index, photo in enumerate(photos, 1):
                    photo_id = id_text(photo.get("id"))
                    photo_id_missing += not bool(photo_id)
                    links.append((origin, line, index, obs, photo_id))
                    link_count += 1
                if len(observations) == 5000:
                    db.executemany("INSERT INTO api_observations VALUES (?,?,?,?,?)", observations)
                    db.executemany("INSERT INTO api_photo_links VALUES (?,?,?,?,?)", links)
                    observations.clear()
                    links.clear()
        db.executemany("INSERT INTO api_observations VALUES (?,?,?,?,?)", observations)
        db.executemany("INSERT INTO api_photo_links VALUES (?,?,?,?,?)", links)
        if observation_count != checkpoint["observations_seen"]:
            raise ValueError("Raw API count does not reconcile with the collection checkpoint")
        report.update(raw_api_status="ARCHIVED_RECORDS_ENUMERATED", raw_api_sha256_uncompressed=raw_hash.hexdigest(),
                      raw_api_observations=observation_count, raw_photo_array_entries=link_count,
                      raw_observations_missing_id=absent_id, raw_observations_without_photos=no_photos,
                      raw_photo_entries_missing_id=photo_id_missing)
    return report


def reconcile(archive_root: Path, metadata: Path, out: Path, manifest: dict, contract: dict) -> dict:
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Output exists; preserve earlier runs")
    if digest(metadata) != contract["source"]["member_sha256"]:
        raise ValueError("Merged source identity mismatch")
    out.mkdir(parents=True)
    try:
        with sqlite3.connect(out / "source_reconciliation.sqlite") as db:
            db.executescript("""
                CREATE TABLE original_photo_rows (origin TEXT, source_row INTEGER, obs_id TEXT, photo_id TEXT, payload_sha256 TEXT, PRIMARY KEY(origin, source_row));
                CREATE TABLE merged_photo_rows (source_row INTEGER PRIMARY KEY, origin TEXT, obs_id TEXT, photo_id TEXT, payload_sha256 TEXT);
                CREATE TABLE api_observations (origin TEXT, source_line INTEGER, obs_id TEXT, photo_array_count INTEGER, raw_line_sha256 TEXT, PRIMARY KEY(origin, source_line));
                CREATE TABLE api_photo_links (origin TEXT, source_line INTEGER, photo_index INTEGER, obs_id TEXT, photo_id TEXT, PRIMARY KEY(origin, source_line, photo_index));
            """)
            chunks = []
            for spec in manifest["archives"]:
                if spec["role"] == "raw_metadata_chunk":
                    chunks.append(ingest_archive(db, archive_root / str(spec["artifact_id"]) / "source.zip", spec))
                    db.commit()
            with metadata.open(encoding="utf-8-sig", newline="") as handle:
                batch = []
                for number, row in csv_records(handle):
                    origin = row.pop("metadata_chunk_source")
                    batch.append((number, origin, row["obs_id"], row["photo_id"], canonical_digest(row)))
                    if len(batch) == 20000:
                        db.executemany("INSERT INTO merged_photo_rows VALUES (?,?,?,?,?)", batch)
                        batch.clear()
                db.executemany("INSERT INTO merged_photo_rows VALUES (?,?,?,?,?)", batch)
            db.executescript("""
                CREATE INDEX original_link ON original_photo_rows(obs_id, photo_id);
                CREATE INDEX merged_link ON merged_photo_rows(obs_id, photo_id);
                CREATE INDEX api_link ON api_photo_links(obs_id, photo_id);
                CREATE TABLE source_links AS SELECT obs_id, photo_id FROM original_photo_rows WHERE obs_id<>'' AND photo_id<>'' UNION SELECT obs_id, photo_id FROM api_photo_links WHERE obs_id<>'' AND photo_id<>'';
                CREATE TABLE source_observations AS SELECT obs_id FROM original_photo_rows WHERE obs_id<>'' UNION SELECT obs_id FROM api_observations WHERE obs_id<>'';
            """)
            scalar = lambda sql: db.execute(sql).fetchone()[0]
            counts = {
                "original_metadata_rows": scalar("SELECT COUNT(*) FROM original_photo_rows"),
                "merged_metadata_rows": scalar("SELECT COUNT(*) FROM merged_photo_rows"),
                "original_metadata_observations": scalar("SELECT COUNT(DISTINCT obs_id) FROM original_photo_rows WHERE obs_id<>''"),
                "merged_metadata_observations": scalar("SELECT COUNT(DISTINCT obs_id) FROM merged_photo_rows WHERE obs_id<>''"),
                "original_versions_absent_from_merge": scalar("SELECT COUNT(*) FROM original_photo_rows o WHERE NOT EXISTS (SELECT 1 FROM merged_photo_rows m WHERE o.origin=m.origin AND o.payload_sha256=m.payload_sha256)"),
                "original_unique_links_absent_from_merge": scalar("SELECT COUNT(*) FROM (SELECT DISTINCT obs_id,photo_id FROM original_photo_rows o WHERE NOT EXISTS (SELECT 1 FROM merged_photo_rows m WHERE o.obs_id=m.obs_id AND o.photo_id=m.photo_id))"),
                "raw_api_records": scalar("SELECT COUNT(*) FROM api_observations"),
                "raw_api_unique_links_absent_from_chunk_metadata": scalar("SELECT COUNT(*) FROM (SELECT DISTINCT a.obs_id,a.photo_id FROM api_photo_links a WHERE a.obs_id<>'' AND a.photo_id<>'' AND NOT EXISTS (SELECT 1 FROM original_photo_rows o WHERE o.obs_id=a.obs_id AND o.photo_id=a.photo_id))"),
                "retained_unique_observation_photo_links": scalar("SELECT COUNT(*) FROM source_links"),
                "retained_unique_photo_ids": scalar("SELECT COUNT(DISTINCT photo_id) FROM source_links"),
                "retained_unique_observations": scalar("SELECT COUNT(*) FROM source_observations"),
                "photo_ids_with_multiple_observation_links": scalar("SELECT COUNT(*) FROM (SELECT photo_id FROM source_links GROUP BY photo_id HAVING COUNT(*)>1)"),
                "source_links_absent_from_merged_metadata": scalar("SELECT COUNT(*) FROM source_links s WHERE NOT EXISTS (SELECT 1 FROM merged_photo_rows m WHERE s.obs_id=m.obs_id AND s.photo_id=m.photo_id)"),
                "source_observations_without_valid_photo_link": scalar("SELECT COUNT(*) FROM source_observations o WHERE NOT EXISTS (SELECT 1 FROM source_links s WHERE o.obs_id=s.obs_id)"),
            }
            if counts["merged_metadata_rows"] != contract["source"]["expected_photo_rows"] or counts["merged_metadata_observations"] != contract["source"]["expected_observations_with_photos"]:
                raise ValueError("Merged counts differ from pinned snapshot")
        report = {
            "status": "ARCHIVED_SOURCE_LINKS_RECONCILED_WITH_EXPLICIT_COVERAGE_GAPS",
            "counts": counts,
            "chunks": chunks,
            "raw_api_missing_archives": [c["artifact_id"] for c in chunks if c["raw_api_member"] is None],
            "manifest_sha256_canonical_json": canonical_digest(manifest),
            "contract_sha256_canonical_json": canonical_digest(contract),
            "merged_metadata_sha256": contract["source"]["member_sha256"],
            "implementation_sha256_text_lf": text_digest(Path(__file__)),
            "ledger_sha256": digest(out / "source_reconciliation.sqlite"),
            "source_rows_removed_in_this_reconciliation": 0,
            "image_bytes_verified": False,
            "ecological_models_executed": False,
        }
        (out / "source_reconciliation_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        return report
    except (Exception, KeyboardInterrupt) as error:
        (out / "incomplete_run.json").write_text(json.dumps({"status":"INCOMPLETE_DO_NOT_USE","error_type":type(error).__name__,"reason":str(error)}, indent=2) + "\n", encoding="utf-8")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=ROOT / "analysis/v3/upstream_sources.json")
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    result = reconcile(args.archives, args.metadata, args.out_dir, json.loads(args.manifest.read_text()), json.loads(args.contract.read_text()))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
