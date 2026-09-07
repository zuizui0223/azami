"""Recount pinned, unthinned historical processing without rerunning models.

The local ledger links source photos, queue jobs, screening, detections and
post-prediction views. Historical processing is not current image-cache coverage.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import sqlite3
import zipfile

from .reconcile_sources import csv_records
from .workflow import CONTRACT, ROOT, canonical_digest, digest, text_digest


def insert_rows(db, table, records, size):
    statement = f"INSERT INTO {table} VALUES ({','.join(['?'] * size)})"
    batch = []
    for row in records:
        batch.append(row)
        if len(batch) == 20000:
            db.executemany(statement, batch)
            batch.clear()
    db.executemany(statement, batch)


def archive_rows(zipped, member):
    with zipped.open(member) as binary, io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed historical CSV: " + member)
            yield row


def audit_history(archive, metadata, out, spec, contract):
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Output exists; preserve earlier runs")
    if archive.stat().st_size != spec["archive_bytes"] or digest(archive) != spec["archive_sha256"]:
        raise ValueError("Processing archive identity mismatch")
    if digest(metadata) != contract["source"]["member_sha256"]:
        raise ValueError("Merged source identity mismatch")
    out.mkdir(parents=True)
    try:
        with sqlite3.connect(out / "historical_processing.sqlite") as db, zipfile.ZipFile(archive) as zipped:
            if len(zipped.namelist()) != len(set(zipped.namelist())):
                raise ValueError("Duplicate archive member names")
            db.executescript("""
                CREATE TABLE source_photos (photo_id TEXT PRIMARY KEY, obs_id TEXT);
                CREATE INDEX source_obs ON source_photos(obs_id);
                CREATE TABLE queue (queue_id TEXT PRIMARY KEY, photo_id TEXT UNIQUE, obs_id TEXT, taxon_name TEXT);
                CREATE INDEX queue_obs ON queue(obs_id);
                CREATE TABLE screens (queue_id TEXT PRIMARY KEY, status TEXT, n_detections INTEGER);
                CREATE TABLE heads (head_id TEXT PRIMARY KEY, queue_id TEXT, photo_id TEXT, obs_id TEXT, det_index TEXT,
                    colour_status TEXT, shape_status TEXT, orientation_status TEXT);
                CREATE TABLE crops (queue_id TEXT, det_index TEXT, photo_id TEXT, obs_id TEXT, PRIMARY KEY(queue_id, det_index));
                CREATE TABLE cohort_members (cohort TEXT, obs_id TEXT, taxon_name TEXT, PRIMARY KEY(cohort, obs_id));
            """)
            with metadata.open(encoding="utf-8-sig", newline="") as handle:
                insert_rows(db, "source_photos", ((r["photo_id"], r["obs_id"]) for _, r in csv_records(handle)), 2)
            prefix = "exhaustive_merged/"
            insert_rows(db, "queue", ((r["queue_id"], r["photo_id"], r["obs_id"], r["taxon_name"])
                                      for r in archive_rows(zipped, prefix + "exhaustive_photo_queue_merged.csv")), 4)
            insert_rows(db, "screens", ((r["queue_id"], r["screen_status"], int(r["n_detections"]))
                                        for r in archive_rows(zipped, prefix + "exhaustive_screening_results.csv")), 3)
            print(json.dumps({"stage": "historical_queue_and_screens_loaded"}), flush=True)
            insert_rows(db, "heads", ((r["annotation_unit_id"], r["queue_id"], r["photo_id"], r["obs_id"], r["det_index"],
                                      r["colour_status"], r["shape_status"], r["orientation_status"])
                                      for r in archive_rows(zipped, prefix + "exhaustive_continuous_head_level.csv")), 8)
            insert_rows(db, "crops", ((r["queue_id"], r["det_index"], r["photo_id"], r["obs_id"])
                                      for r in archive_rows(zipped, prefix + "exhaustive_yolo_crop_metadata.csv")), 4)
            db.executescript("""
                CREATE UNIQUE INDEX head_detection ON heads(queue_id, det_index);
                CREATE INDEX head_obs ON heads(obs_id);
                CREATE INDEX head_photo ON heads(photo_id);
            """)
            scalar = lambda sql: db.execute(sql).fetchone()[0]
            checks = {
                "queue_source_identity_mismatches": "SELECT COUNT(*) FROM queue q LEFT JOIN source_photos s ON q.photo_id=s.photo_id WHERE s.photo_id IS NULL OR q.obs_id<>s.obs_id",
                "queue_without_screen": "SELECT COUNT(*) FROM queue q WHERE NOT EXISTS (SELECT 1 FROM screens s WHERE q.queue_id=s.queue_id)",
                "screen_without_queue": "SELECT COUNT(*) FROM screens s WHERE NOT EXISTS (SELECT 1 FROM queue q WHERE q.queue_id=s.queue_id)",
                "invalid_screen_states": "SELECT COUNT(*) FROM screens WHERE status NOT IN ('detected','no_detection','missing_image','error') OR n_detections<0",
                "head_queue_identity_mismatches": "SELECT COUNT(*) FROM heads h LEFT JOIN queue q ON h.queue_id=q.queue_id WHERE q.queue_id IS NULL OR h.photo_id<>q.photo_id OR h.obs_id<>q.obs_id",
                "head_crop_identity_mismatches": "SELECT COUNT(*) FROM heads h LEFT JOIN crops c ON h.queue_id=c.queue_id AND h.det_index=c.det_index WHERE c.queue_id IS NULL OR h.photo_id<>c.photo_id OR h.obs_id<>c.obs_id",
                "crop_without_head": "SELECT COUNT(*) FROM crops c WHERE NOT EXISTS (SELECT 1 FROM heads h WHERE h.queue_id=c.queue_id AND h.det_index=c.det_index)",
                "screen_head_count_mismatches": "SELECT COUNT(*) FROM screens s LEFT JOIN (SELECT queue_id,COUNT(*) n FROM heads GROUP BY queue_id) h ON s.queue_id=h.queue_id WHERE s.n_detections<>COALESCE(h.n,0) OR (s.status='detected')<>(COALESCE(h.n,0)>0)",
            }
            results = {name: scalar(sql) for name, sql in checks.items()}
            if any(results.values()):
                raise ValueError("Historical identity/state reconciliation failed: " + json.dumps(results))
            counts = {
                "n_queue_photos": scalar("SELECT COUNT(*) FROM queue"),
                "n_queue_observations": scalar("SELECT COUNT(DISTINCT obs_id) FROM queue"),
                "n_queue_species": scalar("SELECT COUNT(DISTINCT taxon_name) FROM queue"),
                "n_detected_photos": scalar("SELECT COUNT(*) FROM screens WHERE status='detected'"),
                "n_no_detection_photos": scalar("SELECT COUNT(*) FROM screens WHERE status='no_detection'"),
                "n_failed_screening_photos": scalar("SELECT COUNT(*) FROM screens WHERE status IN ('missing_image','error')"),
                "n_observations_with_detected_heads": scalar("SELECT COUNT(DISTINCT obs_id) FROM heads"),
                "n_detected_heads": scalar("SELECT COUNT(*) FROM heads"),
                "n_colour_usable": scalar("SELECT COUNT(*) FROM heads WHERE colour_status='usable'"),
                "n_shape_usable": scalar("SELECT COUNT(*) FROM heads WHERE shape_status='usable'"),
                "n_orientation_usable": scalar("SELECT COUNT(*) FROM heads WHERE orientation_status='usable'"),
            }
            saved = json.loads(zipped.read(prefix + "exhaustive_continuous_merge_report.json"))
            if any(saved[k] != value for k, value in counts.items()):
                raise ValueError("Recount differs from historical merge report")
            counts.update(source_photos_not_in_historical_queue=scalar("SELECT COUNT(*) FROM source_photos s WHERE NOT EXISTS (SELECT 1 FROM queue q WHERE s.photo_id=q.photo_id)"),
                          queued_observations_with_multiple_photos=scalar("SELECT COUNT(*) FROM (SELECT obs_id FROM queue GROUP BY obs_id HAVING COUNT(*)>1)"),
                          detected_observations_with_multiple_photos=scalar("SELECT COUNT(*) FROM (SELECT obs_id FROM heads GROUP BY obs_id HAVING COUNT(DISTINCT photo_id)>1)"))
            cohort_report = []
            for member in sorted(zipped.namelist()):
                if member.startswith("postprediction_cohorts/") and member.endswith("observations.csv"):
                    name = Path(member).stem
                    insert_rows(db, "cohort_members", ((name, r["obs_id"], r["taxon_name"]) for r in archive_rows(zipped, member)), 3)
                    count, taxa = db.execute("SELECT COUNT(*),COUNT(DISTINCT taxon_name) FROM cohort_members WHERE cohort=?", (name,)).fetchone()
                    cohort_report.append({"cohort": name, "observations": count, "taxa": taxa})
            if scalar("SELECT COUNT(*) FROM cohort_members c WHERE NOT EXISTS (SELECT 1 FROM heads h WHERE h.obs_id=c.obs_id)"):
                raise ValueError("Post-prediction cohort contains observations without detected heads")
            source_counts = {"photos": scalar("SELECT COUNT(*) FROM source_photos"), "observations": scalar("SELECT COUNT(DISTINCT obs_id) FROM source_photos")}
            if source_counts != {"photos": contract["source"]["expected_photo_rows"], "observations": contract["source"]["expected_observations_with_photos"]}:
                raise ValueError("Merged snapshot count mismatch")
            screen_counts = dict(db.execute("SELECT status,COUNT(*) FROM screens GROUP BY status ORDER BY status"))
            images = sum(Path(name).suffix.lower() in {".jpg", ".jpeg", ".png", ".tif", ".tiff", ".webp"} for name in zipped.namelist())
        report = {
            "status": "HISTORICAL_PROCESSING_RECOUNTED_NOT_V3_REMEASUREMENT",
            "artifact_id": spec["artifact_id"], "archive_sha256": spec["archive_sha256"],
            "merged_source_sha256": contract["source"]["member_sha256"],
            "contract_sha256_canonical_json": canonical_digest(contract),
            "implementation_sha256_text_lf": text_digest(Path(__file__)),
            "counts": counts, "source_counts": source_counts, "screen_status_counts": screen_counts,
            "identity_checks": results, "postprediction_views": cohort_report,
            "historical_report_counts_matched": True,
            "archive_image_members": images, "current_image_cache_availability": "NOT_ASSESSED",
            "new_image_operations": False, "ecological_models_executed": False,
            "ledger_sha256": digest(out / "historical_processing.sqlite"),
            "limits": ["Head counts and technical usability are historical pipeline outputs, not validated biological individuals or accuracy.",
                       "A no-detection result is not an ecological absence; missing/error and not-queued states remain separate.",
                       "Several photos within an observation are not independent samples or automatically repeat measurements of the same head.",
                       "Source links reconstructed from archived API records are tracked separately by reconcile_sources; this comparison uses the exact historical merged snapshot.",
                       "Only the core source/queue/screen/head/crop relations and cohort membership were recounted; trait values and model results were not recomputed."]}
        (out / "historical_processing_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
        return report
    except (Exception, KeyboardInterrupt) as error:
        (out / "incomplete_run.json").write_text(json.dumps({"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__, "reason": str(error)}, indent=2) + "\n", encoding="utf-8")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((ROOT / "analysis/v3/upstream_sources.json").read_text(encoding="utf-8"))
    spec, = [s for s in manifest["archives"] if s["role"] == "historical_all_photo_processing"]
    result = audit_history(args.archive, args.metadata, args.out_dir, spec, json.loads(CONTRACT.read_text(encoding="utf-8")))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
