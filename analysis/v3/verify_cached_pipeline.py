"""Recount saved cached-image products and verify every linked file identity.

This verifies pipeline integrity and reports diagnostic availability. It does
not evaluate whether detections or masks are biologically correct.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3

import numpy as np
from PIL import Image

from .build_image_workspace import readable
from .detect_cached_images import decode_object
from .image_features import registry
from .measure_cached_heads import load_crop
from .workflow import ROOT, digest, text_digest


def require_equal(actual, expected, label):
    if actual != expected:
        raise ValueError(f"{label} mismatch: {actual!r} != {expected!r}")


def verify_mask(root, record):
    path = root / record["path"]
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Mask path escapes measurement directory")
    require_equal(digest(readable(path)), record["sha256"], "Mask SHA-256")
    with Image.open(readable(path)) as image:
        require_equal(image.mode, "L", "Mask image mode")
        require_equal(image.size, (record["width"], record["height"]), "Mask dimensions")
        if not set(np.unique(np.asarray(image))).issubset({0, 255}):
            raise ValueError("Mask is not binary")


def verify(workspace, detection, measurement, out):
    workspace, detection, measurement, out = [p.resolve() for p in (workspace, detection, measurement, out)]
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Preserve earlier verification receipts")
    reports = [json.loads((p / filename).read_text(encoding="utf-8")) for p, filename in (
        (workspace, "image_workspace_report.json"), (detection, "detection_report.json"), (measurement, "measurement_report.json"))]
    paths = [workspace / "image_workspace.sqlite", detection / "detection.sqlite", measurement / "measurements.sqlite"]
    keys = ["workspace_database_sha256", "detection_database_sha256", "measurement_database_sha256"]
    for path, report, key in zip(paths, reports, keys):
        require_equal(digest(path), report[key], key)
    require_equal(reports[1]["execution_contract"]["source_workspace_sha256"], reports[0][keys[0]], "Detector workspace identity")
    require_equal(reports[2]["execution_contract"]["detector_database_sha256"], reports[1][keys[1]], "Measurement detector identity")
    for report in reports[1:]:
        require_equal(report["job_states"].get("pending", 0), 0, "Pending jobs")
    counts = {"verified_image_objects": 0, "verified_crops": 0, "verified_masks": 0,
              "heads_with_flower_and_green_context_pixels": 0, "heads_with_supported_green_context": 0}
    with sqlite3.connect(paths[0].as_uri()+"?mode=ro", uri=True) as cache, sqlite3.connect(paths[1].as_uri()+"?mode=ro", uri=True) as det, sqlite3.connect(paths[2].as_uri()+"?mode=ro", uri=True) as measured:
        objects = list(cache.execute("SELECT sha256,relative_path,width,height,decoded_rgb_sha256 FROM objects"))
        require_equal({r[0] for r in objects}, {r[0] for r in det.execute("SELECT sha256 FROM jobs")}, "Detector image job set")
        for row in objects:
            decode_object(workspace, row)
            counts["verified_image_objects"] += 1
        heads = list(det.execute("SELECT head_id,crops_json FROM detections"))
        require_equal({r[0] for r in heads}, {r[0] for r in measured.execute("SELECT head_id FROM jobs")}, "Measurement head job set")
        for head_id, crops in heads:
            for crop in json.loads(crops).values():
                load_crop(detection, crop)
                counts["verified_crops"] += 1
            actual = {r[0] for r in measured.execute("SELECT endpoint_id FROM endpoints WHERE head_id=?", (head_id,))}
            require_equal(actual, {r["endpoint_id"] for r in registry()}, "Per-head endpoint set")
        for head_id, payload, masks in measured.execute("SELECT * FROM details"):
            for mask in json.loads(masks).values():
                verify_mask(measurement, mask)
                counts["verified_masks"] += 1
            result = json.loads(payload)
            pairs = result["diagnostics"]["paired_colour"]
            counts["heads_with_flower_and_green_context_pixels"] += int(pairs["floral_union"]["n_pixels"] > 0 and pairs["green_non_head_context"]["n_pixels"] > 0)
            counts["heads_with_supported_green_context"] += int(pairs["green_non_head_context"]["support_status"] == "available")
            stored = {r[0]: r[1:] for r in measured.execute("SELECT endpoint_id,value,original,mirror,mirror_abs_difference,status FROM endpoints WHERE head_id=?", (head_id,))}
            for endpoint in result["endpoints"]:
                require_equal(stored[endpoint["endpoint_id"]], tuple(endpoint[k] for k in ("value", "original", "mirror", "mirror_abs_difference", "status")), "Raw/endpoint table agreement")
        counts["detections_retained"] = len(heads)
        counts["endpoint_rows_retained"] = measured.execute("SELECT COUNT(*) FROM endpoints").fetchone()[0]
        counts["finite_qc_failed_values_retained"] = measured.execute("SELECT COUNT(*) FROM endpoints WHERE value IS NOT NULL AND status!='usable'").fetchone()[0]
        counts["heads_without_any_usable_endpoint_retained"] = measured.execute("SELECT COUNT(*) FROM jobs WHERE head_id NOT IN (SELECT head_id FROM endpoints WHERE status='usable')").fetchone()[0]
        counts["negative_image_jobs_retained"] = det.execute("SELECT COUNT(*) FROM jobs WHERE status='no_detection'").fetchone()[0]
        require_equal(counts["endpoint_rows_retained"], len(heads)*27, "Full endpoint denominator")
        require_equal(measured.execute("SELECT COUNT(*) FROM endpoints WHERE status='usable' AND value IS NULL").fetchone()[0], 0, "Usable missing values")
    result = {"status": "CACHED_PIPELINE_FILE_AND_RELATIONSHIP_INTEGRITY_VERIFIED", "counts": counts,
              "input_database_sha256": {key: report[key] for key, report in zip(keys, reports)},
              "implementation_sha256_text_lf": text_digest(Path(__file__)),
              "measurement_job_states": reports[2]["job_states"],
              "limits": ["Byte identity, complete endpoint slots and linked masks are not biological accuracy or full-source coverage.",
                         "Context pixel availability is not a successful negative-control regression or independent illumination calibration.",
                         "Finite QC-failed values remain inspectable; retaining them does not make them eligible for ecological inference."]}
    out.mkdir(parents=True)
    (out / "cached_pipeline_verification.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workspace", "detection", "measurement", "out-dir"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    verify(args.workspace, args.detection, args.measurement, args.out_dir)


if __name__ == "__main__":
    main()
