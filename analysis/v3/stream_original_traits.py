"""Stream original iNaturalist images for the frozen v3 ecological cohort.

The worker is location/environment blind. It downloads one original image at a
time, runs the pinned detector, computes only the 14 measurement-qualified
endpoints for ecological use, evaluates the predeclared +/-5% bbox uncertainty
for bbox-sensitive endpoints, saves numerical/provenance rows, and releases the
source bytes before the next image. Source images are never written to disk.

This module can run a deterministic pilot subset or a deterministic shard. It is
not itself an ecological model and never reads CHELSA values or v2 results.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import time

from . import image_features as features
from .perturb_cached_heads import transform
from .resolution_stream_gate import LICENSES, _decode, _detector_records, _download, sized_url, text
from .workflow import digest


SELECTION_SALT = "ch1-v3-original-stream-pilot-v1"
BBOX_CONDITIONS = [
    {"id": "bbox_left_5pct", "kind": "bbox_shift", "dx": -0.05, "dy": 0.0},
    {"id": "bbox_right_5pct", "kind": "bbox_shift", "dx": 0.05, "dy": 0.0},
    {"id": "bbox_up_5pct", "kind": "bbox_shift", "dx": 0.0, "dy": -0.05},
    {"id": "bbox_down_5pct", "kind": "bbox_shift", "dx": 0.0, "dy": 0.05},
]
IDENTITY = {"id": "baseline", "kind": "identity"}


def finite(value) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def observation_score(obs_id: str) -> str:
    return hashlib.sha256((SELECTION_SALT + "|" + str(obs_id)).encode("utf-8")).hexdigest()


def select_observation_ids(cohort_csv: Path, n: int = 0, shard_index: int = 0, shard_count: int = 1) -> list[str]:
    if n < 0 or shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("Invalid pilot size or shard specification")
    rows = []
    with cohort_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if "obs_id" not in (reader.fieldnames or []):
            raise ValueError("Cohort lacks obs_id")
        for row in reader:
            obs = text(row.get("obs_id"))
            if not obs:
                raise ValueError("Cohort contains missing obs_id")
            score = observation_score(obs)
            if int(score[:16], 16) % shard_count != shard_index:
                continue
            rows.append((score, obs))
    if len({obs for _, obs in rows}) != len(rows):
        raise ValueError("Cohort obs_id is not unique")
    rows.sort()
    if n:
        rows = rows[:n]
    return [obs for _, obs in rows]


def qualified_endpoints(decision: dict) -> tuple[list[str], list[str]]:
    qualified = sorted(row["endpoint_id"] for row in decision["endpoints"] if row.get("ecological_route") == "stream_original_required")
    bbox = sorted(row["endpoint_id"] for row in decision["endpoints"] if row.get("ecological_route") == "stream_original_required" and row.get("bbox_uncertainty_required") is True)
    if len(qualified) != 14:
        raise ValueError(f"Expected 14 original-stream endpoints, found {len(qualified)}")
    if set(bbox) - set(qualified):
        raise ValueError("bbox uncertainty endpoint is not qualified")
    return qualified, bbox


def photo_schedule(metadata_csv: Path, selected_obs: set[str]) -> tuple[list[dict], dict]:
    by_photo: dict[str, dict] = {}
    links = 0
    states = {"license_unavailable": 0, "url_unavailable": 0, "metadata_conflict": 0}
    selected_with_any_link: set[str] = set()
    with metadata_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"obs_id", "photo_id", "photo_license_code", "medium_image_url", "large_image_url"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Photo metadata lacks required scheduling columns")
        for row in reader:
            obs = text(row.get("obs_id"))
            if obs not in selected_obs:
                continue
            selected_with_any_link.add(obs)
            photo = text(row.get("photo_id"))
            if not photo:
                continue
            links += 1
            license_code = text(row.get("photo_license_code")).casefold()
            base = text(row.get("large_image_url")) or text(row.get("medium_image_url")) or text(row.get("raw_image_url"))
            if license_code not in LICENSES:
                states["license_unavailable"] += 1
                continue
            try:
                original_url = sized_url(base, "original") if base else ""
            except ValueError:
                original_url = ""
            if not original_url:
                states["url_unavailable"] += 1
                continue
            candidate = {"photo_id": photo, "original_url": original_url, "license_code": license_code, "obs_ids": [obs]}
            if photo not in by_photo:
                by_photo[photo] = candidate
            else:
                existing = by_photo[photo]
                if existing["original_url"] != original_url or existing["license_code"] != license_code:
                    existing["conflict"] = True
                if obs not in existing["obs_ids"]:
                    existing["obs_ids"].append(obs)
    schedule = []
    for photo, row in by_photo.items():
        if row.pop("conflict", False):
            states["metadata_conflict"] += 1
            continue
        row["obs_ids"].sort()
        schedule.append(row)
    schedule.sort(key=lambda row: (observation_score(row["obs_ids"][0]), row["photo_id"]))
    report = {
        "selected_observations": len(selected_obs),
        "selected_observations_with_metadata_link": len(selected_with_any_link),
        "metadata_observation_photo_links": links,
        "scheduled_unique_photos": len(schedule),
        "scheduled_observation_photo_links": int(sum(len(row["obs_ids"]) for row in schedule)),
        "unavailable_states": states,
    }
    return schedule, report


def endpoint_map(result: dict) -> dict[str, dict]:
    return {row["endpoint_id"]: row for row in result["endpoints"]}


def measure_head(bgr, detections, focal_index: int, qualified: list[str], bbox_sensitive: list[str]) -> tuple[list[dict], list[dict], dict]:
    box = detections[focal_index][0]
    other_boxes = [record[0] for index, record in enumerate(detections) if index != focal_index]
    head, context, recipe = transform(bgr, box, other_boxes, IDENTITY)
    baseline, _ = features.measure(head, context, recipe["context_box"], recipe["all_excluded_head_boxes"])
    baseline_by = endpoint_map(baseline)
    endpoint_rows = [baseline_by[ep] for ep in qualified]

    changes: dict[str, list[float]] = {ep: [] for ep in bbox_sensitive}
    usable_shift_counts: dict[str, int] = {ep: 0 for ep in bbox_sensitive}
    for condition in BBOX_CONDITIONS:
        shifted_head, shifted_context, shifted_recipe = transform(bgr, box, other_boxes, condition)
        shifted, _ = features.measure(
            shifted_head, shifted_context, shifted_recipe["context_box"], shifted_recipe["all_excluded_head_boxes"]
        )
        shifted_by = endpoint_map(shifted)
        for ep in bbox_sensitive:
            base = baseline_by[ep]
            moved = shifted_by[ep]
            if base["status"] == "usable" and moved["status"] == "usable" and finite(base["value"]) and finite(moved["value"]):
                changes[ep].append(abs(float(moved["value"]) - float(base["value"])))
                usable_shift_counts[ep] += 1

    uncertainty_rows = []
    for ep in bbox_sensitive:
        values = changes[ep]
        uncertainty_rows.append({
            "endpoint_id": ep,
            "bbox_shift_usable_n": usable_shift_counts[ep],
            "bbox_shift_median_abs_change": statistics.median(values) if values else None,
            "bbox_shift_max_abs_change": max(values) if values else None,
        })
    diagnostics = baseline.get("diagnostics", {})
    return endpoint_rows, uncertainty_rows, diagnostics


def run(metadata: Path, cohort: Path, weights: Path, decision_path: Path, out: Path,
        pilot_observations: int = 128, shard_index: int = 0, shard_count: int = 1,
        expected_cohort_sha256: str | None = None) -> dict:
    import requests
    from ultralytics import YOLO
    from .detect_cached_images import MODEL_SHA

    if digest(weights) != MODEL_SHA:
        raise ValueError("Detector weight identity mismatch")
    if expected_cohort_sha256 and digest(cohort) != expected_cohort_sha256:
        raise ValueError("Frozen ecological source cohort identity mismatch")
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    if decision.get("trait_environment_results_inspected") is not False:
        raise ValueError("Measurement decision is not phenotype-blind")
    qualified, bbox_sensitive = qualified_endpoints(decision)
    selected = select_observation_ids(cohort, pilot_observations, shard_index, shard_count)
    selected_set = set(selected)
    schedule, schedule_report = photo_schedule(metadata, selected_set)
    if out.exists():
        raise ValueError("Use a fresh output directory")
    out.mkdir(parents=True)

    with (out / "selected_observations_private.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["obs_id", "selection_sha256"])
        writer.writerows((obs, observation_score(obs)) for obs in selected)

    endpoint_fields = ["photo_id", "head_index", "endpoint_id", "unit", "value", "original", "mirror", "mirror_abs_difference", "status"]
    uncertainty_fields = ["photo_id", "head_index", "endpoint_id", "bbox_shift_usable_n", "bbox_shift_median_abs_change", "bbox_shift_max_abs_change"]
    link_fields = ["photo_id", "obs_id"]
    transfer_fields = ["photo_id", "status", "bytes", "width", "height", "detections", "elapsed_seconds", "error"]

    endpoint_usable = {ep: 0 for ep in qualified}
    bbox_coverage = {ep: {"heads": 0, "four_usable_shifts": 0} for ep in bbox_sensitive}
    context_support = {"non_head_context_available": 0, "green_non_head_context_available": 0, "heads": 0}
    transfer = {"success": 0, "error": 0, "detector_positive": 0, "no_detection": 0, "bytes": 0, "heads": 0}
    max_live_bytes = 0
    started = time.perf_counter()

    model = YOLO(str(weights))
    session = requests.Session()
    session.headers.update({"User-Agent": "azami-ch1-v3-original-stream/1.0"})

    with (out / "endpoint_measurements_private.csv").open("w", encoding="utf-8", newline="") as endpoint_handle, \
         (out / "bbox_uncertainty_private.csv").open("w", encoding="utf-8", newline="") as uncertainty_handle, \
         (out / "photo_observation_links_private.csv").open("w", encoding="utf-8", newline="") as link_handle, \
         (out / "transfer_private.csv").open("w", encoding="utf-8", newline="") as transfer_handle:
        endpoint_writer = csv.DictWriter(endpoint_handle, fieldnames=endpoint_fields)
        uncertainty_writer = csv.DictWriter(uncertainty_handle, fieldnames=uncertainty_fields)
        link_writer = csv.DictWriter(link_handle, fieldnames=link_fields)
        transfer_writer = csv.DictWriter(transfer_handle, fieldnames=transfer_fields)
        endpoint_writer.writeheader(); uncertainty_writer.writeheader(); link_writer.writeheader(); transfer_writer.writeheader()

        for item in schedule:
            for obs in item["obs_ids"]:
                link_writer.writerow({"photo_id": item["photo_id"], "obs_id": obs})
            payload = None
            t0 = time.perf_counter()
            try:
                payload = _download(session, item["original_url"])
                max_live_bytes = max(max_live_bytes, len(payload))
                image, bgr = _decode(payload)
                detections = _detector_records(model, bgr)
                transfer["success"] += 1
                transfer["bytes"] += len(payload)
                if detections:
                    transfer["detector_positive"] += 1
                else:
                    transfer["no_detection"] += 1
                transfer["heads"] += len(detections)
                for head_index in range(len(detections)):
                    endpoints, uncertainty, diagnostics = measure_head(bgr, detections, head_index, qualified, bbox_sensitive)
                    for row in endpoints:
                        endpoint_writer.writerow({"photo_id": item["photo_id"], "head_index": head_index, **{k: row.get(k) for k in endpoint_fields if k not in {"photo_id", "head_index"}}})
                        if row["status"] == "usable" and finite(row["value"]):
                            endpoint_usable[row["endpoint_id"]] += 1
                    for row in uncertainty:
                        uncertainty_writer.writerow({"photo_id": item["photo_id"], "head_index": head_index, **row})
                        bbox_coverage[row["endpoint_id"]]["heads"] += 1
                        if row["bbox_shift_usable_n"] == 4:
                            bbox_coverage[row["endpoint_id"]]["four_usable_shifts"] += 1
                    paired = diagnostics.get("paired_colour", {})
                    context_support["heads"] += 1
                    if paired.get("non_head_context", {}).get("support_status") == "available":
                        context_support["non_head_context_available"] += 1
                    if paired.get("green_non_head_context", {}).get("support_status") == "available":
                        context_support["green_non_head_context_available"] += 1
                transfer_writer.writerow({"photo_id": item["photo_id"], "status": "success", "bytes": len(payload),
                                          "width": image.width, "height": image.height, "detections": len(detections),
                                          "elapsed_seconds": time.perf_counter()-t0, "error": ""})
            except Exception as error:
                transfer["error"] += 1
                transfer_writer.writerow({"photo_id": item["photo_id"], "status": "error", "bytes": 0, "width": "", "height": "",
                                          "detections": "", "elapsed_seconds": time.perf_counter()-t0,
                                          "error": type(error).__name__ + ": " + str(error)})
            finally:
                payload = None

    elapsed = time.perf_counter() - started
    report = {
        "status": "ORIGINAL_STREAM_PILOT_COMPLETE_NO_ECOLOGICAL_MODEL",
        "selection": {"salt": SELECTION_SALT, "pilot_observations": pilot_observations, "shard_index": shard_index, "shard_count": shard_count},
        "cohort_sha256": digest(cohort),
        "metadata_sha256": digest(metadata),
        "measurement_decision_sha256": digest(decision_path),
        "qualified_endpoints": qualified,
        "bbox_uncertainty_endpoints": bbox_sensitive,
        "schedule": schedule_report,
        "transfer": transfer,
        "endpoint_usable_heads": endpoint_usable,
        "bbox_uncertainty_coverage": bbox_coverage,
        "paired_context_support": context_support,
        "elapsed_seconds": elapsed,
        "mean_seconds_per_scheduled_photo": elapsed / len(schedule) if schedule else None,
        "mean_original_bytes_per_successful_photo": transfer["bytes"] / transfer["success"] if transfer["success"] else None,
        "maximum_source_image_bytes_live_at_once": max_live_bytes,
        "source_images_persisted": 0,
        "environment_values_read": 0,
        "ecological_models_executed": 0,
        "technical_uncertainty": {
            "bbox": "per-head +/-5% shift summaries are retained for orientation and gross-shape endpoints",
            "colour": "paired flower/non-head/green-context diagnostics are retained per head; population photometric stress sensitivity is supplied by the frozen 14-condition technical audit rather than re-running six synthetic photometric transforms for every original image"
        },
        "limits": [
            "The pilot measures execution throughput and endpoint yield, not ecological associations.",
            "No source image bytes are written to the output directory.",
            "Detector localization remains pseudo-label-trained; streamed measurement does not create independent detector accuracy validation.",
            "Metadata scheduling does not invent a phenotype for unavailable, unlicensed, failed-download or no-detection photos.",
        ],
    }
    (out / "original_stream_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    parser.add_argument("--decision", type=Path, default=Path("analysis/v3/measurement_qualification_decision_20260908.json"))
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--pilot-observations", type=int, default=128)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--expected-cohort-sha256")
    args = parser.parse_args()
    run(args.metadata, args.cohort, args.weights, args.decision, args.out_dir,
        args.pilot_observations, args.shard_index, args.shard_count, args.expected_cohort_sha256)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
