"""Prepare phenotype-blind full-source observations for v3 environment design.

Only the six immutable iNaturalist metadata/API source archives are read. No
image, detector, trait or fitted-result artifact is admitted. One reversible row
is produced per source observation with taxon/date/public-location provenance.
"""
from __future__ import annotations

import argparse
import csv
from datetime import date
import gzip
import hashlib
import io
import json
from itertools import groupby
from pathlib import Path
import sqlite3
import zipfile

import pandas as pd

from .prepare_observation_annotations import api_fields, csv_fields, consensus, annotate
from .workflow import ROOT, canonical_digest, digest

DEFAULT_MANIFEST = ROOT / "analysis" / "v3" / "environment_source_archives.json"


def verify_archive(path: Path, spec: dict) -> zipfile.ZipFile:
    if path.stat().st_size != int(spec["archive_bytes"]):
        raise ValueError(f"Archive byte count mismatch: {spec['artifact_id']}")
    if digest(path) != spec["archive_sha256"]:
        raise ValueError(f"Archive SHA-256 mismatch: {spec['artifact_id']}")
    zipped = zipfile.ZipFile(path)
    names = zipped.namelist()
    if len(names) != len(set(names)):
        zipped.close()
        raise ValueError("Duplicate archive members")
    member = spec["required_member"]
    with zipped.open(member) as handle:
        member_sha = hashlib.file_digest(handle, "sha256").hexdigest()
    if member_sha != spec["required_member_sha256"]:
        zipped.close()
        raise ValueError(f"Metadata member SHA-256 mismatch: {spec['artifact_id']}")
    return zipped


def save_versions(db: sqlite3.Connection, rows: list[tuple]) -> None:
    db.executemany(
        "INSERT OR IGNORE INTO versions(obs_id,kind,origin,source_row,version_hash,fields_json) VALUES (?,?,?,?,?,?)",
        rows,
    )


def ingest_archive(db: sqlite3.Connection, archive: Path, spec: dict) -> dict:
    zipped = verify_archive(archive, spec)
    origin = str(spec["artifact_id"])
    counts = {"metadata": 0, "api": 0}
    try:
        batch: list[tuple] = []
        with zipped.open(spec["required_member"]) as binary, io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if "obs_id" not in (reader.fieldnames or []):
                raise ValueError("Metadata member lacks obs_id")
            for line, row in enumerate(reader, 1):
                if None in row or any(value is None for value in row.values()) or not row.get("obs_id"):
                    raise ValueError(f"Malformed metadata row {line} in {origin}")
                fields = csv_fields(row)
                batch.append((str(row["obs_id"]), "metadata", origin, line, canonical_digest(fields), json.dumps(fields, ensure_ascii=False, sort_keys=True)))
                counts["metadata"] += 1
                if len(batch) >= 10000:
                    save_versions(db, batch)
                    batch.clear()
            save_versions(db, batch)

        raw_members = [name for name in zipped.namelist() if name in {"observation_raw.ndjson", "observation_raw.ndjson.gz"}]
        if len(raw_members) > 1:
            raise ValueError(f"Ambiguous raw API member in {origin}")
        if raw_members:
            raw_name = raw_members[0]
            with zipped.open(raw_name) as binary:
                stream = gzip.GzipFile(fileobj=binary) if raw_name.endswith(".gz") else binary
                batch = []
                for line, payload in enumerate(stream, 1):
                    record = json.loads(payload)
                    obs = str(record.get("id") or "").strip()
                    if not obs:
                        raise ValueError(f"Raw API row without observation id in {origin}:{line}")
                    fields = api_fields(record)
                    batch.append((obs, "api", origin, line, canonical_digest(fields), json.dumps(fields, ensure_ascii=False, sort_keys=True)))
                    counts["api"] += 1
                    if len(batch) >= 10000:
                        save_versions(db, batch)
                        batch.clear()
                save_versions(db, batch)
    finally:
        zipped.close()
    return {"artifact_id": int(spec["artifact_id"]), **counts}


def source_observation_rows(db: sqlite3.Connection):
    cursor = db.execute(
        "SELECT obs_id,kind,fields_json FROM versions ORDER BY obs_id,kind,version_hash"
    )
    for obs_id, records in groupby(cursor, key=lambda row: row[0]):
        records = list(records)
        preferred, n_versions, fields, conflicts = consensus(
            [(row[1], json.loads(row[2])) for row in records]
        )
        derived = annotate(fields, conflicts, preferred)
        observed_on = str(fields.get("observed_on") or "")
        observation_month = None
        if derived["date_status"] == "exact_day":
            observation_month = date.fromisoformat(observed_on).month
        yield {
            "obs_id": obs_id,
            "preferred_source_kind": preferred,
            "n_candidate_versions": int(n_versions),
            "conflicted_fields_json": json.dumps(conflicts, ensure_ascii=False),
            "source_taxon_name": derived["source_taxon_name"],
            "source_taxon_rank": derived["source_taxon_rank"],
            "observed_on": observed_on,
            "observation_month": observation_month,
            "date_status": derived["date_status"],
            "analysis_latitude": derived["analysis_latitude"],
            "analysis_longitude": derived["analysis_longitude"],
            "coordinate_status": derived["coordinate_status"],
            "position_accuracy_status": derived["position_accuracy_status"],
            "position_accuracy_m": derived["position_accuracy_m"],
            "captive_state": derived["captive_state"],
        }


def build(archives: Path, out_dir: Path, manifest_path: Path = DEFAULT_MANIFEST) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = manifest["expected"]
    if int(manifest.get("trait_files_in_scope", -1)) != 0:
        raise ValueError("Environment source manifest must admit zero trait files")
    if out_dir.exists():
        raise ValueError("Output exists; preserve prior source-only run")
    out_dir.mkdir(parents=True)
    try:
        db_path = out_dir / "environment_source_observations.sqlite"
        with sqlite3.connect(db_path) as db:
            db.execute(
                "CREATE TABLE versions(obs_id TEXT,kind TEXT,origin TEXT,source_row INTEGER,version_hash TEXT,fields_json TEXT,PRIMARY KEY(obs_id,kind,origin,source_row))"
            )
            ingested = []
            for spec in manifest["archives"]:
                if spec.get("role") != "raw_metadata_chunk":
                    raise ValueError("Non-source archive admitted to environment source manifest")
                archive = archives / str(spec["artifact_id"]) / "source.zip"
                ingested.append(ingest_archive(db, archive, spec))
                db.commit()
                print(json.dumps(ingested[-1]), flush=True)

            metadata_rows = int(db.execute("SELECT COUNT(*) FROM versions WHERE kind='metadata'").fetchone()[0])
            api_rows = int(db.execute("SELECT COUNT(*) FROM versions WHERE kind='api'").fetchone()[0])
            unique_observations = int(db.execute("SELECT COUNT(DISTINCT obs_id) FROM versions").fetchone()[0])
            if metadata_rows != int(expected["metadata_rows"]):
                raise ValueError(f"Metadata denominator changed: {metadata_rows}")
            if api_rows != int(expected["api_rows"]):
                raise ValueError(f"API denominator changed: {api_rows}")
            if unique_observations != int(expected["unique_observations"]):
                raise ValueError(f"Observation denominator changed: {unique_observations}")

            rows = list(source_observation_rows(db))
            if len(rows) != unique_observations:
                raise ValueError("Consensus observation count differs from source universe")
            frame = pd.DataFrame(rows)
            frame.to_sql("observations", db, index=False, if_exists="replace")
            db.execute("CREATE UNIQUE INDEX observation_id ON observations(obs_id)")
            db.commit()

        csv_path = out_dir / "environment_source_observations.csv"
        frame.to_csv(csv_path, index=False)
        states = {
            "date_status": {str(k): int(v) for k, v in frame["date_status"].value_counts(dropna=False).items()},
            "coordinate_status": {str(k): int(v) for k, v in frame["coordinate_status"].value_counts(dropna=False).items()},
            "captive_state": {str(k): int(v) for k, v in frame["captive_state"].value_counts(dropna=False).items()},
            "source_taxon_rank": {str(k): int(v) for k, v in frame["source_taxon_rank"].value_counts(dropna=False).items()},
        }
        report = {
            "status": "PHENOTYPE_BLIND_FULL_SOURCE_OBSERVATIONS_PREPARED",
            "manifest_sha256_canonical_json": canonical_digest(manifest),
            "counts": {
                "metadata_rows": metadata_rows,
                "api_rows": api_rows,
                "unique_observations": unique_observations,
                "exact_date_rows": int(frame["date_status"].eq("exact_day").sum()),
                "public_coordinate_rows": int(frame["coordinate_status"].eq("public_location_present_precision_not_gated").sum()),
            },
            "states": states,
            "source_database_sha256": digest(db_path),
            "source_csv_sha256": digest(csv_path),
            "trait_files_read": 0,
            "image_files_read": 0,
            "ecological_models_executed": 0,
            "limits": [
                "The first source chunk has no archived raw API stream; its flattened metadata remains an explicit fallback.",
                "Source taxon assignment is not independently verified from the image.",
                "Public coordinates are retained without yet deciding model-resolution positional-accuracy eligibility.",
                "Observation date identifies the source date, not independently verified flowering stage.",
            ],
        }
        (out_dir / "environment_source_observations_report.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        return report
    except Exception as error:
        (out_dir / "incomplete_run.json").write_text(
            json.dumps({"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__, "reason": str(error)}, indent=2) + "\n",
            encoding="utf-8",
        )
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build(args.archives.resolve(), args.out_dir.resolve(), args.manifest)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
