"""Run the pinned detector on all verified cached image objects, resumably.

One job per exact byte object, retaining all photo/observation links in the input
workspace. No development labels or taxonomic/geographic fields are read.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import io
import json
import math
import os
from pathlib import Path
import platform
import sqlite3

import numpy as np
from PIL import Image, ImageOps

from .build_image_workspace import readable
from .workflow import ROOT, digest, text_digest

MODEL_SHA = "4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed"
PARAMETERS = {"conf": 0.25, "iou": 0.70, "imgsz": 640, "max_det": 300,
              "device": "cpu", "augment": False, "verbose": False}


def crop_box(box, width, height, padding):
    if width <= 0 or height <= 0 or not math.isfinite(padding) or padding < 0:
        raise ValueError("Invalid image extent or crop padding")
    if len(box) != 4 or not all(math.isfinite(v) for v in box):
        raise ValueError("Nonfinite detector box")
    x1, y1, x2, y2 = box
    if x2 <= x1 or y2 <= y1:
        raise ValueError("Nonpositive detector box")
    raw = (math.floor(x1 - (x2-x1)*padding), math.floor(y1-(y2-y1)*padding),
           math.ceil(x2+(x2-x1)*padding), math.ceil(y2+(y2-y1)*padding))
    clipped = (max(0, raw[0]), max(0, raw[1]), min(width, raw[2]), min(height, raw[3]))
    if clipped[2] <= clipped[0] or clipped[3] <= clipped[1]:
        raise ValueError("Detector box lies outside the image")
    return clipped, raw != clipped


def save_png(image, path):
    target = readable(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    handle = io.BytesIO()
    image.save(handle, format="PNG", optimize=False)
    payload = handle.getvalue()
    sha = hashlib.sha256(payload).hexdigest()
    if target.exists():
        if digest(target) != sha:
            raise ValueError("Existing crop differs; do not overwrite")
    else:
        temporary = readable(path.with_suffix(".writing"))
        with temporary.open("xb") as sink:
            sink.write(payload)
        temporary.rename(target)
    return sha


def decode_object(workspace, row):
    sha, relative, width, height, pixel_sha = row
    path = workspace / relative
    if not path.resolve().is_relative_to(workspace.resolve()):
        raise ValueError("Image path escapes the source workspace")
    if digest(readable(path)) != sha:
        raise ValueError("Cached image byte identity changed")
    with Image.open(readable(path)) as opened:
        image = ImageOps.exif_transpose(opened).convert("RGB")
        image.load()
    actual_pixels = hashlib.sha256(image.width.to_bytes(8, "big") + image.height.to_bytes(8, "big") + image.tobytes()).hexdigest()
    if image.size != (width, height) or actual_pixels != pixel_sha:
        raise ValueError("Decoded image differs from the cache inventory")
    return image


def write_detections(db, out, sha, image, records):
    invalid = 0
    for index, (box, confidence, class_id) in enumerate(records):
        identifier = f"{sha}:{index}"
        reason = None
        details = {}
        try:
            if class_id != 0 or not math.isfinite(confidence) or not 0 <= confidence <= 1:
                raise ValueError("Unexpected detector class or score")
            extents = {kind: crop_box(box, image.width, image.height, padding)
                       for kind, padding in (("head", .12), ("context", .8))}
        except ValueError as error:
            invalid += 1
            reason = str(error)
        if reason is None:
            # A changed existing crop is an integrity failure, not a bad ROI.
            for kind, (extent, clipped) in extents.items():
                relative = Path("crops") / sha / f"{index:04d}_{kind}.png"
                crop_hash = save_png(image.crop(extent), out / relative)
                details[kind] = {"path": relative.as_posix(), "sha256": crop_hash, "box": extent, "clipped": clipped}
        db.execute("INSERT INTO detections VALUES (?,?,?,?,?,?,?,?,?)", (identifier, sha, index, json.dumps(box), confidence, class_id,
                   "invalid_roi" if reason else "roi_ready", reason, json.dumps(details, sort_keys=True)))
    return invalid


def run(workspace, weights, out, limit=0):
    workspace, out = workspace.resolve(), out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if limit < 0:
        raise ValueError("Batch limit cannot be negative")
    if digest(readable(weights)) != MODEL_SHA:
        raise ValueError("Detector weight identity mismatch")
    report = json.loads((workspace / "image_workspace_report.json").read_text(encoding="utf-8"))
    source_db = workspace / "image_workspace.sqlite"
    if digest(source_db) != report["workspace_database_sha256"]:
        raise ValueError("Image workspace database changed")
    out.mkdir(parents=True, exist_ok=True)
    (out / "runtime_settings").mkdir(exist_ok=True)
    os.environ["YOLO_CONFIG_DIR"] = str(out / "runtime_settings")
    os.environ["YOLO_AUTOINSTALL"] = "false"
    import torch
    from ultralytics import YOLO
    torch.set_num_threads(2)
    torch.manual_seed(20260907)
    torch.use_deterministic_algorithms(True)
    import cv2
    cv2.setNumThreads(1)
    cv2.setRNGSeed(20260907)
    model = YOLO(str(readable(weights)))
    context = {"source_workspace_sha256": report["workspace_database_sha256"], "model_sha256": MODEL_SHA,
               "implementation_sha256_text_lf": text_digest(Path(__file__)), "parameters": PARAMETERS,
               "workspace_helper_sha256_text_lf": text_digest(Path(__file__).with_name("build_image_workspace.py")),
               "python": platform.python_version(),
               "head_padding": .12, "context_padding": .8, "rounding": "floor_left_top_ceil_right_bottom",
               "seed": 20260907, "torch_threads": 2, "opencv_threads": 1,
               "software": {name: importlib.metadata.version(name) for name in ("torch", "torchvision", "ultralytics", "numpy", "pillow", "opencv-python")}}
    lock = out / "execution_contract.json"
    if lock.exists():
        if json.loads(lock.read_text(encoding="utf-8")) != context:
            raise ValueError("Execution context changed; start a new versioned run")
    else:
        if any(p.name != "runtime_settings" for p in out.iterdir()):
            raise ValueError("Unrecognized existing output directory")
        lock.write_text(json.dumps(context, indent=2) + "\n", encoding="utf-8")
    with sqlite3.connect(source_db.as_uri() + "?mode=ro", uri=True) as source, sqlite3.connect(out / "detection.sqlite", timeout=5) as db:
        db.executescript("""
            CREATE TABLE IF NOT EXISTS jobs (sha256 TEXT PRIMARY KEY, status TEXT, n_predictions INTEGER, n_invalid_rois INTEGER, max_det_flag INTEGER, error_type TEXT, error_message TEXT);
            CREATE TABLE IF NOT EXISTS detections (head_id TEXT PRIMARY KEY, sha256 TEXT, det_index INTEGER, box_json TEXT,
                confidence REAL, class_id INTEGER, roi_status TEXT, reason TEXT, crops_json TEXT);
            CREATE INDEX IF NOT EXISTS detection_image ON detections(sha256);
        """)
        objects = {r[0]: r for r in source.execute("SELECT sha256,relative_path,width,height,decoded_rgb_sha256 FROM objects ORDER BY sha256")}
        db.executemany("INSERT OR IGNORE INTO jobs (sha256,status) VALUES (?, 'pending')", ((s,) for s in objects))
        db.commit()
        completed = 0
        while not limit or completed < limit:
            # The write lock covers selecting, processing and committing a job;
            # a second worker cannot claim the same pending object.
            db.execute("BEGIN IMMEDIATE")
            candidate = db.execute("SELECT sha256 FROM jobs WHERE status='pending' ORDER BY sha256 LIMIT 1").fetchone()
            if candidate is None:
                db.commit()
                break
            sha = candidate[0]
            try:
                image = decode_object(workspace, objects[sha])
                predictions = model.predict(np.asarray(image)[:, :, ::-1].copy(), **PARAMETERS)
                if len(predictions) != 1:
                    raise ValueError("Detector returned an unexpected image count")
                boxes = predictions[0].boxes
                records = [] if boxes is None else [(list(map(float, b)), float(c), int(k)) for b,c,k in
                                                    zip(boxes.xyxy.cpu().tolist(), boxes.conf.cpu().tolist(), boxes.cls.cpu().tolist())]
                records.sort(key=lambda r: (-r[1], r[2], *r[0]))
                invalid = write_detections(db, out, sha, image, records)
                status = "detected_with_invalid_rois" if invalid else "detected" if records else "no_detection"
                db.execute("UPDATE jobs SET status=?,n_predictions=?,n_invalid_rois=?,max_det_flag=? WHERE sha256=?",
                           (status, len(records), invalid, int(len(records) >= PARAMETERS["max_det"]), sha))
                db.commit()
            except Exception as error:
                db.rollback()
                db.execute("UPDATE jobs SET status='error',error_type=?,error_message=? WHERE sha256=?", (type(error).__name__, str(error), sha))
                db.commit()
            completed += 1
            if completed % 25 == 0:
                print(json.dumps({"image_jobs_completed_this_invocation": completed}), flush=True)
        states = dict(db.execute("SELECT status,COUNT(*) FROM jobs GROUP BY status ORDER BY status"))
        counts = {"input_image_objects": len(objects), "completed_this_invocation": completed,
                  "raw_detections": db.execute("SELECT COUNT(*) FROM detections").fetchone()[0],
                  "valid_head_context_pairs": db.execute("SELECT COUNT(*) FROM detections WHERE roi_status='roi_ready'").fetchone()[0],
                  "possible_max_det_truncations": db.execute("SELECT COUNT(*) FROM jobs WHERE max_det_flag=1").fetchone()[0]}
    status = ("CACHED_IMAGE_DETECTION_PARTIAL" if states.get("pending") else
              "CACHED_IMAGE_DETECTION_COMPLETED_WITH_ERRORS" if states.get("error") else "CACHED_IMAGE_DETECTION_COMPLETED")
    result = {"status": status,
              "execution_contract": context, "counts": counts, "job_states": states,
              "source_photo_ids": report["counts"]["source_photo_ids"],
              "source_photo_ids_without_verified_image_in_this_workspace": report["counts"]["photo_ids_without_verified_local_image"],
              "detection_database_sha256": digest(out / "detection.sqlite"), "ecological_models_executed": False,
              "limits": ["Cached-image execution is not full-source processing or detector accuracy.",
                         "The production detector was trained against automatic pseudo-labels; no independent precision/recall is estimated here.",
                         "All detector proposals, invalid ROIs, negative jobs and errors remain distinct; confidence is not calibrated correctness.",
                         "EXIF display orientation is not gravity. PNG crops are lossless derivatives; source objects remain unchanged.",
                         "Raw head/context crops do not establish finalized continuous traits or biological stage."]}
    (out / "detection_report.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0, help="Per-invocation work limit; remaining jobs stay pending. Zero processes all pending cached objects.")
    args = parser.parse_args()
    run(args.workspace, args.weights, args.out_dir, args.limit)


if __name__ == "__main__":
    main()
