"""Source-only, location-blind photo jobs from the full reconciled link ledger."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sqlite3
from urllib.parse import urlsplit

from .enriched_source_cohort import COHORT_FIELDS, EXTRA_FIELDS, pinned_hash
from .prepare_environment_source_observations import verify_archive
from .recover_native_source_authority import private_directory, write_new_json
from .resolution_stream_gate import LICENSES, sized_url, text
from .workflow import ROOT, canonical_digest, digest, text_digest

CONTRACT = ROOT / "analysis/v3/reconciled_photo_schedule_contract.json"
STATUS = "RECONCILED_NATIVE_PHOTO_SCHEDULE_VERIFIED_NO_IMAGE_EXECUTION"
HOSTS = {"static.inaturalist.org", "inaturalist-open-data.s3.amazonaws.com"}


def component_score(component: str, salt: str) -> str:
    return hashlib.sha256((salt + "|" + component).encode()).hexdigest()


def photo_fields(row: dict, kind: str) -> dict:
    if kind == "metadata":
        urls = [text(row.get(k)) for k in ("raw_image_url", "small_image_url", "medium_image_url", "large_image_url")]
        return {"license_code": text(row.get("photo_license_code")).casefold(),
                "attribution": text(row.get("photo_attribution")), "urls": urls}
    return {"license_code": text(row.get("license_code")).casefold(),
            "attribution": text(row.get("attribution")),
            "urls": [text(row.get(k)) for k in ("url", "small_url", "medium_url", "large_url", "original_url")]}


def original_identity(photo: str, fields: dict) -> tuple[str, str]:
    urls = [url for url in fields["urls"] if url]
    if not urls:
        return "url_unavailable", ""
    normalized = set()
    for url in urls:
        try:
            parsed = urlsplit(url)
            if (parsed.scheme not in {"http", "https"} or parsed.hostname not in HOSTS
                    or parsed.username or parsed.password or parsed.port or parsed.query or parsed.fragment
                    or not re.fullmatch(r"/photos/" + re.escape(photo) + r"/(square|small|medium|large|original)\.[A-Za-z0-9]+", parsed.path)):
                return "url_invalid_or_identity_conflict", ""
            normalized.add(sized_url(url, "original"))
        except ValueError:
            return "url_invalid_or_identity_conflict", ""
    if len(normalized) != 1:
        return "url_conflict", ""
    return "valid", next(iter(normalized))


def read_native_projection(path: Path) -> list[tuple]:
    rows = []
    seen = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or []) != set(COHORT_FIELDS + EXTRA_FIELDS):
            raise ValueError("Enriched source schema differs; no arbitrary trait input is allowed")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Duplicate enriched source fields")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed enriched source row")
            obs, component = row["obs_id"], row["dependence_component_id"]
            if (not obs or obs in seen or not component or row["native_range_status"] != "native"
                    or row["source_taxon_rank"] not in {"species", "subspecies", "variety"}):
                raise ValueError("Invalid or duplicate native source membership")
            count = int(row["reconciled_source_photo_count"])
            if count < 1:
                raise ValueError("Native source row has no reconciled photo support")
            seen.add(obs)
            rows.append((obs, component, count))
    if not rows:
        raise ValueError("Empty native source cohort")
    return sorted(rows)


def ingest_versions(db: sqlite3.Connection, source: sqlite3.Connection, archives: Path,
                    specs: list[dict], target_photos: set[str]) -> list[dict]:
    reports = []
    for spec in specs:
        origin = str(spec["artifact_id"])
        archive = archives / origin / "source.zip"
        counts = Counter()
        batch = []

        def save(kind, number, index, obs, photo, record_hash, raw_photo):
            if photo not in target_photos:
                return
            fields = photo_fields(raw_photo, kind)
            state, url = original_identity(photo, fields)
            batch.append((photo, obs, kind, origin, number, index, record_hash,
                          canonical_digest(raw_photo), fields["license_code"], state, url,
                          json.dumps(fields, sort_keys=True, ensure_ascii=True)))
            counts["retained_" + kind + "_photo_versions"] += 1
            if len(batch) >= 10000:
                db.executemany("INSERT INTO photo_versions VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", batch)
                batch.clear()

        with verify_archive(archive, spec) as zipped:
            expected = source.execute("SELECT source_row,obs_id,photo_id,payload_sha256 FROM original_photo_rows WHERE origin=? ORDER BY source_row", (origin,))
            with zipped.open(spec["required_member"]) as binary, io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as handle:
                for number, row in enumerate(csv.DictReader(handle), 1):
                    if None in row or any(v is None for v in row.values()):
                        raise ValueError("Malformed archived photo metadata")
                    record_hash = canonical_digest(row)
                    if expected.fetchone() != (number, row["obs_id"], row["photo_id"], record_hash):
                        raise ValueError("Archived metadata differs from reconciled source locator/hash")
                    counts["metadata_rows_checked"] += 1
                    save("metadata", number, 0, row["obs_id"], row["photo_id"], record_hash, row)
            if expected.fetchone() is not None:
                raise ValueError("Archived metadata ended before reconciled rows")
            raw_members = [n for n in zipped.namelist() if n in {"observation_raw.ndjson", "observation_raw.ndjson.gz"}]
            if len(raw_members) > 1:
                raise ValueError("Ambiguous archived API stream")
            expected_obs = source.execute("SELECT source_line,obs_id,photo_array_count,raw_line_sha256 FROM api_observations WHERE origin=? ORDER BY source_line", (origin,))
            expected_links = source.execute("SELECT source_line,photo_index,obs_id,photo_id FROM api_photo_links WHERE origin=? ORDER BY source_line,photo_index", (origin,))
            if raw_members:
                member = raw_members[0]
                with zipped.open(member) as binary:
                    stream = gzip.GzipFile(fileobj=binary) if member.endswith(".gz") else binary
                    for number, payload in enumerate(stream, 1):
                        row = json.loads(payload)
                        obs = text(row.get("id"))
                        photos = row.get("photos") or []
                        record_hash = hashlib.sha256(payload).hexdigest()
                        if expected_obs.fetchone() != (number, obs, len(photos), record_hash):
                            raise ValueError("Archived API observation differs from reconciled locator/hash")
                        counts["api_observations_checked"] += 1
                        for index, photo in enumerate(photos, 1):
                            photo_id = text(photo.get("id"))
                            if expected_links.fetchone() != (number, index, obs, photo_id):
                                raise ValueError("Archived API photo link differs from reconciliation")
                            counts["api_photo_links_checked"] += 1
                            save("api", number, index, obs, photo_id, record_hash, photo)
            if expected_obs.fetchone() is not None or expected_links.fetchone() is not None:
                raise ValueError("Archived API ended before reconciled records")
        db.executemany("INSERT INTO photo_versions VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", batch)
        db.commit()
        reports.append({"artifact_id": int(origin), **counts})
        print(json.dumps({"stage": "schedule_source_versions_checked", **reports[-1]}), flush=True)
    return reports


def finalize_jobs(db: sqlite3.Connection, salt: str, shard_count: int) -> None:
    for photo, component in db.execute("SELECT l.photo_id,MIN(o.component_id) FROM native_links l JOIN native_observations o USING(obs_id) GROUP BY l.photo_id ORDER BY l.photo_id"):
        licenses = set()
        url_states, urls = set(), set()
        n_versions = 0
        for license_code, url_state, url in db.execute("SELECT license_code,url_state,original_url FROM photo_versions WHERE photo_id=?", (photo,)):
            licenses.add(license_code)
            url_states.add(url_state)
            urls.add(url)
            n_versions += 1
        if not n_versions:
            state = "source_photo_version_missing"
        elif len(licenses) != 1:
            state = "license_conflict"
        elif not licenses.issubset(LICENSES):
            state = "license_unavailable"
        elif url_states != {"valid"}:
            state = "url_unavailable_or_invalid"
        elif len(urls) != 1:
            state = "url_conflict"
        else:
            state = "request_candidate_not_authorized"
        score = component_score(component, salt)
        db.execute("INSERT INTO photo_jobs VALUES (?,?,?,?,?,?,?,?)", (
            photo, component, score, int(score[:16], 16) % shard_count, state,
            next(iter(licenses)) if len(licenses) == 1 else "",
            next(iter(urls)) if state == "request_candidate_not_authorized" else "", n_versions))
    db.commit()


def build(enriched: Path, reconciliation: Path, archives: Path, out: Path, *,
          expected_enriched_sha256: str, expected_reconciliation_sha256: str,
          manifest_path: Path = ROOT / "analysis/v3/environment_source_archives.json",
          contract_path: Path = CONTRACT, shard_count: int = 64) -> dict:
    out = private_directory(out)
    if out.exists() or shard_count < 1:
        raise ValueError("Fresh output and a positive scheduling shard count are required")
    hashes = {"enriched": pinned_hash(enriched, expected_enriched_sha256, "enriched cohort"),
              "reconciliation": pinned_hash(reconciliation, expected_reconciliation_sha256, "reconciliation")}
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract["production_image_execution_authorized"] is not False or set(contract["license_codes"]) != LICENSES:
        raise ValueError("Photo schedule contract cannot promote production or alter inherited license codes")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("trait_files_in_scope") != 0 or any(s.get("role") != "raw_metadata_chunk" for s in manifest["archives"]):
        raise ValueError("Schedule source manifest must contain raw metadata archives only")
    sibling = enriched.with_name("enriched_source_cohort_report.json")
    if not sibling.exists():
        raise ValueError("Enriched source receipt required")
    enrichment = json.loads(sibling.read_text(encoding="utf-8"))
    if enrichment["output_csv_sha256"] != hashes["enriched"] or enrichment["input_sha256"]["reconciliation"] != hashes["reconciliation"]:
        raise ValueError("Enriched receipt input/output identity mismatch")
    projection = read_native_projection(enriched)
    out.mkdir(parents=True)
    db_path = out / "reconciled_photo_schedule_private.sqlite"
    try:
        with sqlite3.connect(reconciliation.resolve().as_uri() + "?mode=ro", uri=True) as source, sqlite3.connect(db_path, uri=True) as db:
            db.execute("ATTACH DATABASE ? AS source", (reconciliation.resolve().as_uri() + "?mode=ro",))
            db.executescript("""
                CREATE TABLE native_observations (obs_id TEXT PRIMARY KEY,component_id TEXT NOT NULL,expected_photo_count INTEGER NOT NULL);
                CREATE TABLE photo_versions (photo_id TEXT,obs_id TEXT,kind TEXT,origin TEXT,source_row INTEGER,photo_index INTEGER,source_record_sha256 TEXT,photo_payload_sha256 TEXT,license_code TEXT,url_state TEXT,original_url TEXT,photo_fields_json TEXT,PRIMARY KEY(kind,origin,source_row,photo_index));
                CREATE TABLE photo_jobs (photo_id TEXT PRIMARY KEY,component_id TEXT,component_score TEXT,shard_index INTEGER,state TEXT,license_code TEXT,original_url TEXT,source_version_count INTEGER);
            """)
            db.executemany("INSERT INTO native_observations VALUES (?,?,?)", projection)
            db.executescript("""
                CREATE TABLE native_links AS SELECT s.obs_id,s.photo_id FROM source.source_links s JOIN native_observations o USING(obs_id);
                CREATE UNIQUE INDEX native_link ON native_links(obs_id,photo_id);
                CREATE INDEX native_photo ON native_links(photo_id);
                CREATE TABLE known_photo_links AS SELECT s.obs_id,s.photo_id FROM source.source_links s WHERE s.photo_id IN (SELECT photo_id FROM native_links);
                CREATE UNIQUE INDEX known_link ON known_photo_links(obs_id,photo_id);
                CREATE INDEX known_photo ON known_photo_links(photo_id);
            """)
            mismatches = db.execute("SELECT COUNT(*) FROM native_observations o WHERE expected_photo_count != (SELECT COUNT(*) FROM native_links l WHERE l.obs_id=o.obs_id)").fetchone()[0]
            split = db.execute("SELECT COUNT(*) FROM (SELECT photo_id FROM native_links JOIN native_observations USING(obs_id) GROUP BY photo_id HAVING COUNT(DISTINCT component_id)>1)").fetchone()[0]
            if mismatches or split or len(projection) != enrichment["cohort_rows"]:
                raise ValueError("Native observation/photo membership or full-source dependence mismatch")
            targets = {r[0] for r in db.execute("SELECT DISTINCT photo_id FROM native_links")}
            chunks = ingest_versions(db, source, archives, manifest["archives"], targets)
            db.execute("CREATE INDEX photo_version ON photo_versions(photo_id)")
            missing = db.execute("SELECT COUNT(*) FROM known_photo_links l WHERE NOT EXISTS(SELECT 1 FROM photo_versions v WHERE v.photo_id=l.photo_id AND v.obs_id=l.obs_id)").fetchone()[0]
            unexpected = db.execute("SELECT COUNT(*) FROM photo_versions v WHERE NOT EXISTS(SELECT 1 FROM known_photo_links l WHERE l.photo_id=v.photo_id AND l.obs_id=v.obs_id)").fetchone()[0]
            if missing or unexpected:
                raise ValueError("Target-photo version ledger does not conserve all known links")
            finalize_jobs(db, contract["ordering_salt"], shard_count)
            counts = {table: db.execute("SELECT COUNT(*) FROM " + table).fetchone()[0] for table in
                      ("native_observations", "native_links", "known_photo_links", "photo_versions", "photo_jobs")}
            states = dict(db.execute("SELECT state,COUNT(*) FROM photo_jobs GROUP BY state"))
            known_components = db.execute("SELECT COUNT(DISTINCT component_id) FROM native_observations").fetchone()[0]
            candidate_observations = db.execute("SELECT COUNT(DISTINCT l.obs_id) FROM native_links l JOIN photo_jobs j USING(photo_id) WHERE j.state='request_candidate_not_authorized'").fetchone()[0]
            if known_components != enrichment["dependence"]["cohort_known_components"]:
                raise ValueError("Known native dependence component count changed")
            if db.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                raise ValueError("Schedule SQLite integrity failed")
            execution = {"status": STATUS, "input_sha256": hashes, "contract_canonical_sha256": canonical_digest(contract),
                         "archive_manifest_canonical_sha256": canonical_digest(manifest), "shard_count": shard_count,
                         "native_source_rows": len(projection), "production_image_execution_authorized": False}
            db.execute("CREATE TABLE execution (key TEXT PRIMARY KEY,value_json TEXT)")
            db.executemany("INSERT INTO execution VALUES (?,?)", [(k, json.dumps(v, sort_keys=True)) for k, v in execution.items()])
            db.commit()
        report = {"schema_version": 1, **execution, "counts": counts, "photo_states": states,
                  "native_known_components": known_components, "native_observations_with_request_candidates": candidate_observations,
                  "all_native_observations_preserved": True, "all_known_target_photo_links_preserved": True,
                  "shared_photos_partitioned_once": True, "worker_location_taxon_date_fields_present": False,
                  "source_archive_checks": chunks, "output_sqlite_sha256": digest(db_path),
                  "implementation_sha256_text_lf": text_digest(Path(__file__)),
                  "source_rows_deleted": 0, "image_requests_executed": 0, "trait_files_read": 0,
                  "ecological_models_executed": 0, "ecological_fitting_authorized": False,
                  "limits": ["Scheduling state is source support, not measured image availability or endpoint assessability.",
                             "Known dependence components are operational partitions, not independent validation folds.",
                             "Recorded photo-license codes are not independent legal clearance or permission to republish.",
                             "The schedule does not verify an off-device private archive or enable production image execution."]}
        write_new_json(out / "reconciled_photo_schedule_report.json", report)
        return report
    except Exception as error:
        write_new_json(out / "incomplete_run.json", {"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__})
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("enriched", "reconciliation", "archives", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("enriched", "reconciliation"):
        parser.add_argument("--expected-" + name + "-sha256", required=True)
    parser.add_argument("--shard-count", type=int, default=64)
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
