"""Resumable all-endpoint baseline measurement on a pinned cached detector pass.

No taxon, environment or human-label inputs. This is baseline extraction with
mirror QC and paired context, not the completed technical-evaluation protocol.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import sqlite3

import cv2
import numpy as np
from PIL import Image

from .build_image_workspace import readable
from .detect_cached_images import save_png
from . import image_features as features
from .workflow import ROOT, digest, text_digest


def load_crop(root, crop):
    path = root / crop["path"]
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Crop path escapes detector output")
    if digest(readable(path)) != crop["sha256"]:
        raise ValueError("Crop byte identity mismatch")
    with Image.open(readable(path)) as image:
        array = np.asarray(image.convert("RGB"))[:, :, ::-1].copy()
    x1, y1, x2, y2 = crop["box"]
    if array.shape[:2] != (y2-y1, x2-x1):
        raise ValueError("Crop dimensions differ from detector coordinates")
    return array


def insert_endpoints(db, head_id, endpoints):
    expected = {r["endpoint_id"] for r in features.registry()}
    if len(endpoints) != 27 or {e["endpoint_id"] for e in endpoints} != expected:
        raise ValueError("Every detected head requires the full 27-endpoint denominator")
    db.executemany("INSERT INTO endpoints VALUES (?,?,?,?,?,?,?,?)", [
        (head_id, e["endpoint_id"], e["unit"], e["value"], e["original"], e["mirror"], e["mirror_abs_difference"], e["status"]) for e in endpoints])


def failed_endpoints(reason):
    return [{"endpoint_id": row["endpoint_id"], "unit": row["unit"], "value": None,
             "original": None, "mirror": None, "mirror_abs_difference": None, "status": reason} for row in features.registry()]


def software_versions():
    versions = {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "pillow")}
    versions["cv2"] = cv2.__version__
    for name in ("opencv-python", "opencv-python-headless"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return versions


def run(detection, out, limit=0):
    detection, out = detection.resolve(), out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if limit < 0:
        raise ValueError("Negative work limit")
    source_report = json.loads((detection / "detection_report.json").read_text(encoding="utf-8"))
    if source_report["job_states"].get("pending"):
        raise ValueError("Finish the pinned detector pass before measuring it")
    source_db = detection / "detection.sqlite"
    if digest(source_db) != source_report["detection_database_sha256"]:
        raise ValueError("Detector database identity mismatch")
    context = {"detector_database_sha256": source_report["detection_database_sha256"],
               "detector_execution": source_report["execution_contract"], "feature_specification": features.specification(),
               "worker_sha256_text_lf": text_digest(Path(__file__)),
               "python": platform.python_version(),
               "software": software_versions(),
               "io_helpers_sha256_text_lf": {name: text_digest(Path(__file__).with_name(name)) for name in ("detect_cached_images.py", "build_image_workspace.py")},
               "opencv_threads": 1, "seed": 20260907,
               "mask_outputs": "primary head foreground/floral union and paired context masks; extended engine regions remain reconstructible from crops plus pinned code"}
    out.mkdir(parents=True, exist_ok=True)
    lock = out / "execution_contract.json"
    if lock.exists():
        if json.loads(lock.read_text(encoding="utf-8")) != context:
            raise ValueError("Measurement context changed; use a new versioned run")
    else:
        if any(out.iterdir()):
            raise ValueError("Unrecognized existing output directory")
        lock.write_text(json.dumps(context, indent=2) + "\n", encoding="utf-8", newline="\n")
    cv2.setNumThreads(1)
    cv2.setRNGSeed(20260907)
    with sqlite3.connect(source_db.as_uri() + "?mode=ro", uri=True) as source, sqlite3.connect(out / "measurements.sqlite") as db:
        db.executescript("""
            CREATE TABLE IF NOT EXISTS jobs (head_id TEXT PRIMARY KEY, sha256 TEXT, status TEXT, error_type TEXT, error_message TEXT);
            CREATE TABLE IF NOT EXISTS endpoints (head_id TEXT, endpoint_id TEXT, unit TEXT, value REAL, original REAL, mirror REAL, mirror_abs_difference REAL, status TEXT, PRIMARY KEY(head_id,endpoint_id));
            CREATE TABLE IF NOT EXISTS details (head_id TEXT PRIMARY KEY, raw_and_diagnostics_json TEXT, masks_json TEXT);
        """)
        rows = {r[0]: r for r in source.execute("SELECT head_id,sha256,det_index,roi_status,crops_json FROM detections ORDER BY head_id")}
        all_boxes = {}
        for head_id, sha, index, roi_status, crops in rows.values():
            if roi_status == "roi_ready":
                all_boxes.setdefault(sha, []).append(json.loads(crops)["head"]["box"])
            db.execute("INSERT OR IGNORE INTO jobs(head_id,sha256,status) VALUES (?,?, 'pending')", (head_id, sha))
        db.commit()
        completed = 0
        while not limit or completed < limit:
            db.execute("BEGIN IMMEDIATE")
            candidate = db.execute("SELECT head_id FROM jobs WHERE status='pending' ORDER BY head_id LIMIT 1").fetchone()
            if candidate is None:
                db.commit()
                break
            head_id = candidate[0]
            _, sha, index, roi_status, encoded = rows[head_id]
            try:
                if roi_status != "roi_ready":
                    insert_endpoints(db, head_id, failed_endpoints("invalid_roi"))
                    status = "invalid_roi"
                else:
                    crops = json.loads(encoded)
                    head = load_crop(detection, crops["head"])
                    background = load_crop(detection, crops["context"])
                    result, masks = features.measure(head, background, crops["context"]["box"], all_boxes[sha])
                    saved_masks = {}
                    for name, mask in masks.items():
                        relative = Path("masks") / sha / f"{index:04d}_{name}.png"
                        mask_sha = save_png(Image.fromarray(mask.astype(np.uint8)*255), out / relative)
                        saved_masks[name] = {"path": relative.as_posix(), "sha256": mask_sha, "height": mask.shape[0], "width": mask.shape[1]}
                    insert_endpoints(db, head_id, result["endpoints"])
                    db.execute("INSERT INTO details VALUES (?,?,?)", (head_id, json.dumps(result, sort_keys=True, allow_nan=False), json.dumps(saved_masks, sort_keys=True)))
                    status = "measured_with_engine_errors" if result["engine_errors"] else "measured"
                db.execute("UPDATE jobs SET status=? WHERE head_id=?", (status, head_id))
                db.commit()
            except Exception as error:
                db.rollback()
                insert_endpoints(db, head_id, failed_endpoints("worker_error"))
                db.execute("UPDATE jobs SET status='error',error_type=?,error_message=? WHERE head_id=?", (type(error).__name__, str(error), head_id))
                db.commit()
            completed += 1
            if completed % 50 == 0:
                print(json.dumps({"head_jobs_completed_this_invocation": completed}), flush=True)
        states = dict(db.execute("SELECT status,COUNT(*) FROM jobs GROUP BY status ORDER BY status"))
        counts = {"input_detections": len(rows), "completed_this_invocation": completed,
                  "endpoint_rows_retained": db.execute("SELECT COUNT(*) FROM endpoints").fetchone()[0],
                  "heads_with_any_usable_endpoint": db.execute("SELECT COUNT(DISTINCT head_id) FROM endpoints WHERE status='usable'").fetchone()[0],
                  "heads_with_all_27_usable": db.execute("SELECT COUNT(*) FROM (SELECT head_id FROM endpoints WHERE status='usable' GROUP BY head_id HAVING COUNT(*)=27)").fetchone()[0]}
        coverage = {r["endpoint_id"]: {"finite_raw_means": db.execute("SELECT COUNT(*) FROM endpoints WHERE endpoint_id=? AND value IS NOT NULL", (r["endpoint_id"],)).fetchone()[0],
                     "usable_heads": db.execute("SELECT COUNT(*) FROM endpoints WHERE endpoint_id=? AND status='usable'", (r["endpoint_id"],)).fetchone()[0]} for r in features.registry()}
        if counts["endpoint_rows_retained"] != (len(rows)-states.get("pending", 0))*27:
            raise ValueError("Incomplete per-head endpoint denominator")
    status = ("CACHED_27_BASELINE_PARTIAL" if states.get("pending") else "CACHED_27_BASELINE_COMPLETED_WITH_ERRORS" if states.get("error") or states.get("measured_with_engine_errors") else "CACHED_27_BASELINE_COMPLETED")
    report = {"status": status, "execution_contract": context, "counts": counts, "job_states": states,
              "endpoint_coverage": coverage, "measurement_database_sha256": digest(out / "measurements.sqlite"),
              "ecological_models_executed": False, "technical_perturbation_evaluation_completed": False,
              "independent_accuracy_estimated": False,
              "limits": ["Cached historical images are not a representative full-source sample or a newly independent evaluation set.",
                         "All 27 slots are retained, including nulls, errors, QC failures and finite values below eligibility thresholds.",
                         "Mirror QC and baseline extraction do not replace crop, resolution and photometric perturbation evaluation.",
                         "Paired context masks exclude all detected head boxes, not all true flowers; green pixels are not verified leaf tissue.",
                         "Context colour does not alone separate photographic illumination from background ecology or prove flower-specific biology.",
                         "EXIF vertical, JPEG-derived colour and normalized geometry remain image features, not calibrated physical or botanical quantities.",
                         "The 27-endpoint denominator includes joint hue and a closed four-part composition; it is not 27 independent tests."]}
    (out / "measurement_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(report, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detection", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    run(args.detection, args.out_dir, args.limit)


if __name__ == "__main__":
    main()
