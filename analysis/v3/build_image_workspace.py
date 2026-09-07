"""Retain the full photo universe while importing verifiable local image bytes.

No source-photo downloads, detector inference or ecological analysis. Image
objects are content-addressed; source photo IDs, versions and usage pools are
separate relationships. Workers need not read geographic or taxonomic metadata.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import sqlite3

from PIL import Image, ImageOps

from .workflow import CONTRACT, ROOT, canonical_digest, digest, text_digest


def readable(path: Path) -> Path:
    """Support existing Windows archive paths longer than MAX_PATH."""
    absolute = str(path.resolve())
    if os.name == "nt" and not absolute.startswith("\\\\?\\"):
        return Path("\\\\?\\UNC\\" + absolute[2:] if absolute.startswith("\\\\") else "\\\\?\\" + absolute)
    return Path(absolute)


def rows(path):
    with readable(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed input manifest")
            yield row


def unique_map(records, key):
    result = {}
    for row in records:
        value = row[key]
        if not value or value in result:
            raise ValueError("Missing or duplicate manifest key: " + key)
        result[value] = row
    return result


def safe_filename(value):
    if not value or value in (".", "..") or Path(value).name != value or any(c in value for c in "/\\:"):
        raise ValueError("Unsafe image filename in source manifest")
    return value


def known_local_sources(materials: Path, perturbation: Path | None):
    """Explicit adapters for the archived local caches; do not guess file joins."""
    sha_manifest = materials / "LOCAL_SHA256_MANIFEST.json"
    saved = json.loads(readable(sha_manifest).read_text(encoding="utf-8"))
    expected = unique_map(saved["files"], "relative_path")
    sources = [sha_manifest]
    declarations = []

    def add(photo_id, path, pool, manifest, source_row):
        relative = path.relative_to(materials).as_posix() if path.is_relative_to(materials) else None
        expected_sha = expected.get(relative, {}).get("sha256") if relative else None
        declarations.append({"photo_id": photo_id, "path": str(path), "pool": pool,
                             "expected_sha256": expected_sha, "source_manifest": str(manifest),
                             "source_row": source_row})

    downloads = materials / "detector_development_qc_images/downloaded_audit/download_results.csv"
    sources.append(downloads)
    for number, row in enumerate(rows(downloads), 1):
        if row["download_status"] != "success":
            raise ValueError("This archived development-cache adapter expects completed downloads")
        basename = safe_filename(PurePosixPath(row["image_local_path"].replace("\\", "/")).name)
        add(row["photo_id"], downloads.parent / "images" / basename, "historical_detector_development_pool", downloads, number)

    key = materials / "detector_audit_private/detector_independent_audit_private_key.csv"
    blinded = materials / "detector_audit_blinded/detector_independent_audit_blinded_manifest.csv"
    sources += [key, blinded]
    by_queue = unique_map(rows(key), "queue_id")
    used = set()
    for number, row in enumerate(rows(blinded), 1):
        private = by_queue[row["queue_id"]]
        if row["audit_id"] != private["audit_id"] or row["queue_id"] in used:
            raise ValueError("Audit-cache identity mismatch")
        used.add(row["queue_id"])
        add(private["photo_id"], blinded.parent / "images" / safe_filename(row["screen_download_filename"]),
            "historical_detector_audit_pool_not_independent_v3_evaluation", blinded, number)
    if used != set(by_queue):
        raise ValueError("Audit key and blinded manifest coverage differ")

    if perturbation is not None:
        sample = perturbation / "perturbation_sample.csv"
        sources.append(sample)
        for number, row in enumerate(rows(sample), 1):
            photo = row["photo_id"]
            add(photo, perturbation / "image_cache" / safe_filename(f"photo_{photo}.jpg"),
                "historical_trait_perturbation_pool", sample, number)
    provenance = [{"path": str(p), "sha256": digest(readable(p))} for p in sources]
    return declarations, provenance


def inspect_and_copy(source: Path, out: Path, expected_sha: str | None):
    if not readable(source).is_file():
        return {"status": "missing_local_file"}
    actual = digest(readable(source))
    if expected_sha and actual != expected_sha:
        return {"status": "source_hash_mismatch", "actual_sha256": actual}
    try:
        with Image.open(readable(source)) as opened:
            opened.verify()
        with Image.open(readable(source)) as opened:
            orientation = opened.getexif().get(274)
            raw_width, raw_height = opened.size
            prepared = ImageOps.exif_transpose(opened).convert("RGB")
            prepared.load()
            width, height = prepared.size
            pixel_hash = hashlib.sha256(width.to_bytes(8, "big") + height.to_bytes(8, "big") + prepared.tobytes()).hexdigest()
    except Exception as error:
        return {"status": "decode_failed", "actual_sha256": actual, "error_type": type(error).__name__}
    relative = Path("objects") / actual[:2] / (actual + ".image")
    target = out / relative
    readable(target.parent).mkdir(parents=True, exist_ok=True)
    if readable(target).exists():
        if digest(readable(target)) != actual:
            raise ValueError("Existing content-addressed object changed")
    else:
        temporary = readable(target.with_suffix(".copying"))
        with readable(source).open("rb") as src, temporary.open("xb") as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
        if digest(temporary) != actual:
            raise ValueError("Image copy integrity failure")
        temporary.rename(readable(target))
    return {"status": "available", "sha256": actual, "relative_path": relative.as_posix(),
            "bytes": readable(target).stat().st_size, "raw_width": raw_width, "raw_height": raw_height,
            "width": width, "height": height, "exif_orientation": orientation,
            "decoded_rgb_sha256": pixel_hash,
            "version_evidence": "matches_archived_local_hash" if expected_sha else "first_hash_inventory_of_existing_local_file"}


def build_workspace(metadata, links, out, declarations, source_provenance, contract, expected_ledger_sha):
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Output exists; do not overwrite earlier workspaces")
    if digest(metadata) != contract["source"]["member_sha256"] or digest(links) != expected_ledger_sha:
        raise ValueError("Source metadata or reconciled-ledger identity mismatch")
    out.mkdir(parents=True)
    try:
        with sqlite3.connect(out / "image_workspace.sqlite") as db:
            db.executescript("""
                CREATE TABLE photos (photo_id TEXT PRIMARY KEY, medium_url TEXT, large_url TEXT, license_code TEXT);
                CREATE TABLE source_links (obs_id TEXT, photo_id TEXT, PRIMARY KEY(obs_id,photo_id));
                CREATE INDEX link_photo ON source_links(photo_id);
                CREATE TABLE objects (sha256 TEXT PRIMARY KEY, relative_path TEXT, bytes INTEGER, raw_width INTEGER, raw_height INTEGER,
                    width INTEGER, height INTEGER, exif_orientation INTEGER, decoded_rgb_sha256 TEXT);
                CREATE TABLE cache_records (record_id INTEGER PRIMARY KEY, photo_id TEXT, pool TEXT, source_path TEXT, expected_sha256 TEXT,
                    source_manifest TEXT, source_row INTEGER, status TEXT, content_sha256 TEXT, version_evidence TEXT);
                CREATE TABLE photo_versions (photo_id TEXT, sha256 TEXT, PRIMARY KEY(photo_id,sha256));
                CREATE INDEX version_content ON photo_versions(sha256);
            """)
            db.executemany("INSERT INTO photos VALUES (?,?,?,?)", ((r["photo_id"], r["medium_image_url"], r["large_image_url"], r["photo_license_code"]) for r in rows(metadata)))
            with sqlite3.connect(links.resolve().as_uri() + "?mode=ro", uri=True) as source:
                db.executemany("INSERT INTO source_links VALUES (?,?)", source.execute("SELECT obs_id,photo_id FROM source_links ORDER BY obs_id,photo_id"))
            if db.execute("SELECT COUNT(*) FROM source_links l WHERE NOT EXISTS (SELECT 1 FROM photos p WHERE p.photo_id=l.photo_id)").fetchone()[0]:
                raise ValueError("Reconciled photo link lacks source metadata")
            for number, record in enumerate(declarations, 1):
                if not db.execute("SELECT 1 FROM photos WHERE photo_id=?", (record["photo_id"],)).fetchone():
                    raise ValueError("Cache declaration photo is outside the recovered source")
                inspected = inspect_and_copy(Path(record["path"]), out, record.get("expected_sha256"))
                if inspected["status"] == "available":
                    db.execute("INSERT OR IGNORE INTO objects VALUES (?,?,?,?,?,?,?,?,?)", tuple(inspected[k] for k in
                               ("sha256", "relative_path", "bytes", "raw_width", "raw_height", "width", "height", "exif_orientation", "decoded_rgb_sha256")))
                    db.execute("INSERT OR IGNORE INTO photo_versions VALUES (?,?)", (record["photo_id"], inspected["sha256"]))
                db.execute("INSERT INTO cache_records VALUES (?,?,?,?,?,?,?,?,?,?)", (number, record["photo_id"], record["pool"],
                           record["path"], record.get("expected_sha256"), record["source_manifest"], record["source_row"],
                           inspected["status"], inspected.get("sha256", inspected.get("actual_sha256")), inspected.get("version_evidence")))
                if number % 250 == 0:
                    print(json.dumps({"local_image_records_checked": number}), flush=True)
            db.execute("CREATE TABLE image_jobs AS SELECT sha256, 'pending_detection' AS status FROM objects ORDER BY sha256")
            db.execute("CREATE UNIQUE INDEX image_job ON image_jobs(sha256)")
            scalar = lambda sql: db.execute(sql).fetchone()[0]
            counts = {"source_photo_ids": scalar("SELECT COUNT(*) FROM photos"),
                      "source_observation_ids": scalar("SELECT COUNT(DISTINCT obs_id) FROM source_links"),
                      "source_observation_photo_links": scalar("SELECT COUNT(*) FROM source_links"),
                      "declared_local_image_records": len(declarations), "verified_image_objects": scalar("SELECT COUNT(*) FROM objects"),
                      "cached_source_photo_ids": scalar("SELECT COUNT(DISTINCT photo_id) FROM photo_versions"),
                      "retained_photo_content_versions": scalar("SELECT COUNT(*) FROM photo_versions"),
                      "photos_with_multiple_cached_versions": scalar("SELECT COUNT(*) FROM (SELECT photo_id FROM photo_versions GROUP BY photo_id HAVING COUNT(*)>1)"),
                      "photo_ids_without_verified_local_image": scalar("SELECT COUNT(*) FROM photos p WHERE NOT EXISTS (SELECT 1 FROM photo_versions v WHERE v.photo_id=p.photo_id)"),
                      "objects_linked_to_multiple_photo_ids": scalar("SELECT COUNT(*) FROM (SELECT sha256 FROM photo_versions GROUP BY sha256 HAVING COUNT(*)>1)"),
                      "decoded_pixel_groups_with_multiple_encodings": scalar("SELECT COUNT(*) FROM (SELECT decoded_rgb_sha256 FROM objects GROUP BY decoded_rgb_sha256 HAVING COUNT(*)>1)"),
                      "objects_with_nonidentity_exif_orientation": scalar("SELECT COUNT(*) FROM objects WHERE exif_orientation NOT IN (1) AND exif_orientation IS NOT NULL"),
                      "source_photo_rows_deleted": 0}
            if counts["source_photo_ids"] != contract["source"]["expected_photo_rows"]:
                raise ValueError("Full photo denominator changed")
            states = dict(db.execute("SELECT status,COUNT(*) FROM cache_records GROUP BY status ORDER BY status"))
            pools = [dict(pool=p, records=n, available=a) for p,n,a in db.execute("SELECT pool,COUNT(*),SUM(status='available') FROM cache_records GROUP BY pool ORDER BY pool")]
        report = {"status": "LOCAL_IMAGE_BYTES_INDEXED_FULL_SOURCE_RETAINED", "counts": counts,
                  "cache_record_states": states, "usage_pools": pools,
                  "metadata_sha256": contract["source"]["member_sha256"], "reconciled_ledger_sha256": expected_ledger_sha,
                  "implementation_sha256_text_lf": text_digest(Path(__file__)),
                  "contract_sha256_canonical_json": canonical_digest(contract),
                  "source_manifest_sha256": [r["sha256"] for r in source_provenance],
                  "workspace_database_sha256": digest(out / "image_workspace.sqlite"),
                  "new_source_photo_downloads": 0, "detector_executed": False, "ecological_models_executed": False,
                  "limits": ["Local cache coverage is not full-source download completeness; unobserved cache locations are outside this audit.",
                             "Historical development/audit/perturbation pools are usage provenance, not new independent v3 evaluation splits.",
                             "All byte versions and links are retained. Shared photos, content or decoded pixels require grouping before splitting or inference.",
                             "Raw bytes are preserved. Display orientation uses recorded EXIF; this is not gravity-referenced biological orientation.",
                             "Cached local copies may be processed under the existing research scope; their presence does not authorize public image redistribution."]}
        (out / "local_cache_provenance.json").write_text(json.dumps(source_provenance, indent=2) + "\n", encoding="utf-8")
        (out / "image_workspace_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
        return report
    except (Exception, KeyboardInterrupt) as error:
        (out / "incomplete_run.json").write_text(json.dumps({"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__, "reason": str(error)}, indent=2) + "\n", encoding="utf-8")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--source-links", type=Path, required=True)
    parser.add_argument("--local-materials-root", type=Path, required=True)
    parser.add_argument("--perturbation-root", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    declarations, sources = known_local_sources(args.local_materials_root.resolve(), args.perturbation_root.resolve() if args.perturbation_root else None)
    receipt = json.loads((ROOT / "reproducibility/v3_source_reconciliation_20260907.json").read_text(encoding="utf-8"))
    report = build_workspace(args.metadata, args.source_links, args.out_dir, declarations, sources,
                             json.loads(CONTRACT.read_text(encoding="utf-8")), receipt["ledger_sha256"])
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
