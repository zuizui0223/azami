"""Stream original iNaturalist images for the frozen v3 ecological cohort.

The worker is location/environment blind. It downloads one original image at a
time, runs the pinned detector, retains all 27 raw endpoint slots with the 14
measurement-qualified routes marked separately, evaluates +/-5% bbox uncertainty
for bbox-sensitive endpoints, saves numerical/provenance rows, and releases the
source bytes before the next image. Source images are never written to disk.

This module can run a deterministic pilot subset or a deterministic shard. It is
not itself an ecological model and never reads CHELSA values or v2 results.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import statistics
import time

from . import image_features as features
from .perturb_cached_heads import transform
from .resolution_stream_gate import LICENSES, _decode, _detector_records, _download, sized_url, text
from .workflow import ROOT, digest, text_digest


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


def photo_schedule(metadata_csv: Path, selected_obs: set[str], link_records: list[dict] | None = None) -> tuple[list[dict], dict]:
    """Retain selected links before admission and check every known photo version.

    The second metadata pass includes conflicting versions attached to an
    unselected observation. A permissive row cannot override a restricted row.
    This is still only a pilot over the supplied (possibly unreconciled) table.
    """
    selected_rows = []
    with metadata_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"obs_id", "photo_id", "photo_license_code", "medium_image_url", "large_image_url"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Photo metadata lacks required scheduling columns")
        for source_row, row in enumerate(reader, start=2):
            obs = text(row.get("obs_id"))
            if obs not in selected_obs:
                continue
            selected_rows.append((source_row, row))
    selected_with_any_link = {text(row.get("obs_id")) for _, row in selected_rows}
    selected_photos = {text(row.get("photo_id")) for _, row in selected_rows} - {""}
    versions: dict[str, set[tuple[str, str]]] = {photo: set() for photo in selected_photos}
    known_links: dict[str, set[str]] = {photo: set() for photo in selected_photos}
    with metadata_csv.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            photo = text(row.get("photo_id"))
            if photo not in selected_photos:
                continue
            license_code = text(row.get("photo_license_code")).casefold()
            base = text(row.get("large_image_url")) or text(row.get("medium_image_url")) or text(row.get("raw_image_url"))
            versions[photo].add((license_code, base))
            if text(row.get("obs_id")):
                known_links[photo].add(text(row["obs_id"]))

    states = {state: 0 for state in ("license_unavailable", "url_unavailable", "metadata_conflict", "photo_id_missing", "metadata_link_missing")}
    by_photo: dict[str, dict] = {}
    retained_links = []
    for source_row, row in selected_rows:
        obs, photo = text(row.get("obs_id")), text(row.get("photo_id"))
        license_code = text(row.get("photo_license_code")).casefold()
        base = text(row.get("large_image_url")) or text(row.get("medium_image_url")) or text(row.get("raw_image_url"))
        original_url = ""
        if not photo:
            status = "photo_id_missing"
        elif len(versions[photo]) > 1:
            status = "metadata_conflict"
        elif license_code not in LICENSES:
            status = "license_unavailable"
        else:
            try:
                original_url = sized_url(base, "original") if base else ""
            except ValueError:
                original_url = ""
            status = "scheduled" if original_url else "url_unavailable"
        retained_links.append({"source_row": source_row, "obs_id": obs, "photo_id": photo,
                               "status": status, "license_code": license_code, "source_url": base,
                               "original_url": original_url,
                               "source_metadata_json": json.dumps(row, ensure_ascii=True),
                               "known_photo_versions_json": json.dumps(sorted(versions.get(photo, set()))),
                               "known_photo_observation_ids_json": json.dumps(sorted(known_links.get(photo, set())))})
        if status != "scheduled":
            states[status] += 1
        else:
            candidate = by_photo.setdefault(photo, {"photo_id": photo, "original_url": original_url,
                "license_code": license_code, "obs_ids": [], "known_metadata_obs_ids": sorted(known_links[photo]),
                "source_versions": sorted(versions[photo])})
            if obs not in candidate["obs_ids"]:
                candidate["obs_ids"].append(obs)
    for obs in sorted(selected_obs - selected_with_any_link):
        retained_links.append({"source_row": "", "obs_id": obs, "photo_id": "", "status": "metadata_link_missing",
                               "license_code": "", "source_url": "", "original_url": "",
                               "source_metadata_json": "{}",
                               "known_photo_versions_json": "[]", "known_photo_observation_ids_json": "[]"})
        states["metadata_link_missing"] += 1
    if link_records is not None:
        link_records.extend(retained_links)
    schedule = list(by_photo.values())
    for row in schedule:
        row["obs_ids"].sort()
    schedule.sort(key=lambda row: (observation_score(row["obs_ids"][0]), row["photo_id"]))
    report = {
        "selected_observations": len(selected_obs),
        "selected_observations_with_metadata_link": len(selected_with_any_link),
        "metadata_observation_photo_links": sum(bool(text(row.get("photo_id"))) for _, row in selected_rows),
        "retained_metadata_rows": len(selected_rows),
        "retained_link_state_rows": len(retained_links),
        "scheduled_unique_photos": len(schedule),
        "scheduled_observation_photo_links": int(sum(len(row["obs_ids"]) for row in schedule)),
        "unavailable_states": states,
    }
    return schedule, report


def endpoint_map(result: dict) -> dict[str, dict]:
    return {row["endpoint_id"]: row for row in result["endpoints"]}


def failure_endpoints(status: str) -> list[dict]:
    return [{"endpoint_id": row["endpoint_id"], "unit": row["unit"], "value": None,
             "original": None, "mirror": None, "mirror_abs_difference": None,
             "status": status} for row in features.registry()]


def measure_head(bgr, detections, focal_index: int, qualified: list[str], bbox_sensitive: list[str]) -> tuple[list[dict], list[dict], dict]:
    box = detections[focal_index][0]
    other_boxes = [record[0] for index, record in enumerate(detections) if index != focal_index]
    head, context, recipe = transform(bgr, box, other_boxes, IDENTITY)
    baseline, _ = features.measure(head, context, recipe["context_box"], recipe["all_excluded_head_boxes"])
    baseline_by = endpoint_map(baseline)
    expected = {row["endpoint_id"] for row in features.registry()}
    if len(baseline["endpoints"]) != 27 or set(baseline_by) != expected:
        raise ValueError("Baseline must retain the complete 27-endpoint registry")
    endpoint_rows = [baseline_by[ep] for ep in sorted(expected)]
    conditions = [{"condition": "baseline", "recipe": recipe,
                   "endpoints": [baseline_by[ep] for ep in bbox_sensitive]}]

    changes: dict[str, list[float]] = {ep: [] for ep in bbox_sensitive}
    usable_shift_counts: dict[str, int] = {ep: 0 for ep in bbox_sensitive}
    for condition in BBOX_CONDITIONS:
        try:
            shifted_head, shifted_context, shifted_recipe = transform(bgr, box, other_boxes, condition)
            shifted, _ = features.measure(
                shifted_head, shifted_context, shifted_recipe["context_box"], shifted_recipe["all_excluded_head_boxes"]
            )
            if set(endpoint_map(shifted)) != expected:
                raise ValueError("Shift endpoint registry is incomplete")
        except Exception as error:
            shifted_recipe = {"condition": condition, "error_type": type(error).__name__, "error": str(error)}
            shifted = {"endpoints": failure_endpoints("bbox_measurement_error")}
        shifted_by = endpoint_map(shifted)
        conditions.append({"condition": condition["id"], "recipe": shifted_recipe,
                           "endpoints": [shifted_by[ep] for ep in bbox_sensitive]})
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
    diagnostics = dict(baseline.get("diagnostics", {}))
    diagnostics["stream_measurement_provenance"] = {
        "baseline_recipe": recipe,
        "detector_box": box,
        "detector_confidence": detections[focal_index][1],
        "detector_class": detections[focal_index][2],
        "raw": baseline.get("raw"), "legacy_qc": baseline.get("legacy_qc"),
        "extended_combined": baseline.get("extended_combined"), "engine_errors": baseline.get("engine_errors"),
        "bbox_conditions": conditions,
    }
    return endpoint_rows, uncertainty_rows, diagnostics


def oriented_pixel_identity(image) -> str:
    """Hash EXIF-oriented RGB uint8 pixels with an explicit dimensional prefix."""
    rgb = image.convert("RGB")
    prefix = b"azami-oriented-RGB-uint8-v1\0" + rgb.width.to_bytes(8, "big") + rgb.height.to_bytes(8, "big")
    return hashlib.sha256(prefix + rgb.tobytes()).hexdigest()


def software_runtime() -> dict:
    """Record this interpreter's runtime, never infer it from requirements text."""
    distributions = {}
    for name in ("torch", "torchvision", "ultralytics", "numpy", "pandas", "pillow",
                 "opencv-python", "opencv-python-headless", "requests"):
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = None
    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "system": platform.system(),
        "machine": platform.machine(),
        "distribution_versions": distributions,
        "loaded_feature_module_versions": {"cv2": features.cv2.__version__, "numpy": features.np.__version__},
        "missing_distribution_metadata": sorted(name for name, version in distributions.items() if version is None),
        "requirements_conformance_verified": False,
    }


def run(metadata: Path, cohort: Path, weights: Path, decision_path: Path, out: Path,
        pilot_observations: int = 128, shard_index: int = 0, shard_count: int = 1,
        expected_cohort_sha256: str | None = None) -> dict:
    # No source bytes may be fetched for an unreconciled, non-durable production pass.
    if pilot_observations == 0:
        raise ValueError("Production stream is blocked: reconciled source/dependence input and durable private numerical destination are not implemented")
    if not 1 <= pilot_observations <= 128:
        raise ValueError("Only a bounded local operational pilot of 1..128 observations is permitted")
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / name) for name in ("local_data", "outputs"))):
        raise ValueError("Private numerical outputs require an external or ignored local_data/outputs directory")
    import requests
    from ultralytics import YOLO
    from .detect_cached_images import MODEL_SHA, PARAMETERS

    if digest(weights) != MODEL_SHA:
        raise ValueError("Detector weight identity mismatch")
    if expected_cohort_sha256 and digest(cohort) != expected_cohort_sha256:
        raise ValueError("Frozen ecological source cohort identity mismatch")
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    if decision.get("trait_environment_results_inspected") is not False:
        raise ValueError("Measurement decision is not phenotype-blind")
    qualified, bbox_sensitive = qualified_endpoints(decision)
    routes = {row["endpoint_id"]: row["ecological_route"] for row in decision["endpoints"]}
    registry_ids = {row["endpoint_id"] for row in features.registry()}
    if len(decision["endpoints"]) != 27 or set(routes) != registry_ids:
        raise ValueError("Decision must retain one route for every registered endpoint")
    selected = select_observation_ids(cohort, pilot_observations, shard_index, shard_count)
    selected_set = set(selected)
    link_records = []
    schedule, schedule_report = photo_schedule(metadata, selected_set, link_records)
    if out.exists():
        raise ValueError("Use a fresh output directory")
    out.mkdir(parents=True)
    execution = {
        "status": "LOCAL_OPERATIONAL_PILOT_NOT_PRODUCTION",
        "model_sha256": MODEL_SHA, "detector_parameters": PARAMETERS,
        "worker_sha256_text_lf": text_digest(Path(__file__)),
        "helper_sha256_text_lf": {name: text_digest(Path(__file__).with_name(name)) for name in
                                   ("resolution_stream_gate.py", "perturb_cached_heads.py", "detect_cached_images.py")},
        "feature_specification": features.specification(),
        "software_runtime": software_runtime(),
        "cohort_sha256": digest(cohort), "metadata_sha256": digest(metadata),
        "measurement_decision_sha256": digest(decision_path),
        "retained_endpoint_count": 27, "ecological_route_endpoint_count": len(qualified),
        "pixel_hash_definition": "SHA256(azami-oriented-RGB-uint8-v1 NUL + width uint64 big-endian + height uint64 big-endian + EXIF-oriented RGB uint8 row-major bytes)",
        "source_reconciliation_verified": False, "durable_private_archive_verified": False,
    }
    (out / "execution_contract.json").write_text(json.dumps(execution, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    with (out / "selected_observations_private.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["obs_id", "selection_sha256"])
        writer.writerows((obs, observation_score(obs)) for obs in selected)

    endpoint_fields = ["photo_id", "head_index", "endpoint_id", "unit", "value", "original", "mirror", "mirror_abs_difference", "status", "ecological_route", "primary_measurement_eligible"]
    uncertainty_fields = ["photo_id", "head_index", "endpoint_id", "bbox_shift_usable_n", "bbox_shift_median_abs_change", "bbox_shift_max_abs_change"]
    bbox_fields = ["photo_id", "head_index", "condition", "endpoint_id", "unit", "value", "original", "mirror", "mirror_abs_difference", "status"]
    link_fields = ["source_row", "photo_id", "obs_id", "status", "license_code", "source_url", "original_url", "source_metadata_json", "known_photo_versions_json", "known_photo_observation_ids_json"]
    transfer_fields = ["photo_id", "status", "bytes", "source_byte_sha256", "oriented_rgb_pixel_sha256", "width", "height", "detections", "elapsed_seconds", "error"]

    endpoint_usable = {ep: 0 for ep in qualified}
    primary_measurement_eligible = {ep: 0 for ep in qualified}
    all_endpoint_usable = {ep: 0 for ep in sorted(registry_ids)}
    bbox_coverage = {ep: {"heads": 0, "four_usable_shifts": 0} for ep in bbox_sensitive}
    context_support = {"non_head_context_available": 0, "green_non_head_context_available": 0, "heads": 0}
    transfer = {"success": 0, "error": 0, "download_success": 0, "detector_positive": 0, "no_detection": 0, "bytes": 0, "heads": 0, "head_measurement_errors": 0}
    max_live_bytes = 0
    started = time.perf_counter()

    model = YOLO(str(weights))
    session = requests.Session()
    session.headers.update({"User-Agent": "azami-ch1-v3-original-stream/1.0"})

    with ExitStack() as stack:
        stack.callback(session.close)
        def csv_writer(name, fields):
            handle = stack.enter_context((out / name).open("w", encoding="utf-8", newline=""))
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            return writer
        endpoint_writer = csv_writer("endpoint_measurements_private.csv", endpoint_fields)
        uncertainty_writer = csv_writer("bbox_uncertainty_private.csv", uncertainty_fields)
        bbox_writer = csv_writer("bbox_measurements_private.csv", bbox_fields)
        link_writer = csv_writer("photo_observation_links_private.csv", link_fields)
        transfer_writer = csv_writer("transfer_private.csv", transfer_fields)
        diagnostic_handle = stack.enter_context((out / "head_diagnostics_private.jsonl").open("w", encoding="utf-8"))
        detection_handle = stack.enter_context((out / "photo_detection_private.jsonl").open("w", encoding="utf-8"))
        link_writer.writerows(link_records)

        for item in schedule:
            payload = None
            t0 = time.perf_counter()
            photo_state = {"photo_id": item["photo_id"], "status": "error", "bytes": 0,
                           "source_byte_sha256": "", "oriented_rgb_pixel_sha256": "", "width": "", "height": "",
                           "detections": "", "elapsed_seconds": 0, "error": ""}
            photo_detection = {**item, "status": "not_detected", "detections": []}
            try:
                payload = _download(session, item["original_url"])
                transfer["download_success"] += 1
                transfer["bytes"] += len(payload)
                photo_state.update(bytes=len(payload), source_byte_sha256=hashlib.sha256(payload).hexdigest())
                max_live_bytes = max(max_live_bytes, len(payload))
                image, bgr = _decode(payload)
                photo_state.update(width=image.width, height=image.height, oriented_rgb_pixel_sha256=oriented_pixel_identity(image))
                detections = _detector_records(model, bgr)
                photo_state["detections"] = len(detections)
                photo_detection.update(status="detected" if detections else "no_detection", detections=[
                    {"head_index": i, "box_xyxy": record[0], "confidence": record[1], "class_id": record[2]}
                    for i, record in enumerate(detections)])
                if detections:
                    transfer["detector_positive"] += 1
                else:
                    transfer["no_detection"] += 1
                transfer["heads"] += len(detections)
                head_errors = 0
                for head_index in range(len(detections)):
                    try:
                        endpoints, uncertainty, diagnostics = measure_head(bgr, detections, head_index, qualified, bbox_sensitive)
                    except Exception as error:
                        head_errors += 1
                        transfer["head_measurement_errors"] += 1
                        endpoints = failure_endpoints("head_measurement_error")
                        uncertainty = [{"endpoint_id": ep, "bbox_shift_usable_n": 0,
                                        "bbox_shift_median_abs_change": None, "bbox_shift_max_abs_change": None} for ep in bbox_sensitive]
                        diagnostics = {"error_type": type(error).__name__, "error": str(error), "stream_measurement_provenance": {
                            "bbox_conditions": [{"condition": c["id"], "recipe": {"condition": c},
                                "endpoints": [row for row in endpoints if row["endpoint_id"] in bbox_sensitive]}
                                for c in [IDENTITY, *BBOX_CONDITIONS]]}}
                    if len(endpoints) != 27 or {row["endpoint_id"] for row in endpoints} != registry_ids:
                        raise ValueError("Per-head endpoint denominator is incomplete")
                    stable_bbox = {row["endpoint_id"] for row in uncertainty if row["bbox_shift_usable_n"] == 4}
                    for row in endpoints:
                        usable = row["status"] == "usable" and finite(row["value"])
                        admitted = row["endpoint_id"] in qualified
                        eligible = bool(admitted and usable and (row["endpoint_id"] not in bbox_sensitive or row["endpoint_id"] in stable_bbox))
                        endpoint_writer.writerow({"photo_id": item["photo_id"], "head_index": head_index,
                            **{k: row.get(k) for k in endpoint_fields if k not in {"photo_id", "head_index", "ecological_route", "primary_measurement_eligible"}},
                            "ecological_route": routes[row["endpoint_id"]], "primary_measurement_eligible": eligible})
                        if eligible:
                            primary_measurement_eligible[row["endpoint_id"]] += 1
                        if usable:
                            all_endpoint_usable[row["endpoint_id"]] += 1
                            if admitted:
                                endpoint_usable[row["endpoint_id"]] += 1
                    for row in uncertainty:
                        uncertainty_writer.writerow({"photo_id": item["photo_id"], "head_index": head_index, **row})
                        bbox_coverage[row["endpoint_id"]]["heads"] += 1
                        if row["bbox_shift_usable_n"] == 4:
                            bbox_coverage[row["endpoint_id"]]["four_usable_shifts"] += 1
                    for condition in diagnostics["stream_measurement_provenance"]["bbox_conditions"]:
                        for row in condition["endpoints"]:
                            bbox_writer.writerow({"photo_id": item["photo_id"], "head_index": head_index, "condition": condition["condition"],
                                                  **{k: row.get(k) for k in bbox_fields if k not in {"photo_id", "head_index", "condition"}}})
                    diagnostic_handle.write(json.dumps(features.clean({"photo_id": item["photo_id"], "head_index": head_index,
                        "detector_box_xyxy": detections[head_index][0], "detector_confidence": detections[head_index][1],
                        "detector_class_id": detections[head_index][2], "diagnostics": diagnostics}), allow_nan=False) + "\n")
                    paired = diagnostics.get("paired_colour", {})
                    context_support["heads"] += 1
                    if paired.get("non_head_context", {}).get("support_status") == "available":
                        context_support["non_head_context_available"] += 1
                    if paired.get("green_non_head_context", {}).get("support_status") == "available":
                        context_support["green_non_head_context_available"] += 1
                if head_errors:
                    photo_state.update(status="measurement_error", error=f"{head_errors} head measurement(s) failed; raw failure slots retained")
                    transfer["error"] += 1
                else:
                    photo_state["status"] = "success"
                    transfer["success"] += 1
            except Exception as error:
                transfer["error"] += 1
                photo_state["error"] = type(error).__name__ + ": " + str(error)
            finally:
                photo_state["elapsed_seconds"] = time.perf_counter() - t0
                transfer_writer.writerow(photo_state)
                photo_detection["transfer"] = photo_state
                detection_handle.write(json.dumps(features.clean(photo_detection), allow_nan=False) + "\n")
                payload = None

    elapsed = time.perf_counter() - started
    report = {
        "status": "ORIGINAL_STREAM_PILOT_COMPLETE_NO_ECOLOGICAL_MODEL",
        "selection": {"salt": SELECTION_SALT, "pilot_observations": pilot_observations, "shard_index": shard_index, "shard_count": shard_count},
        "cohort_sha256": digest(cohort),
        "metadata_sha256": digest(metadata),
        "measurement_decision_sha256": digest(decision_path),
        "qualified_endpoints": qualified,
        "retained_endpoints": sorted(registry_ids),
        "retained_endpoints_per_detected_head": 27,
        "production_execution_authorized": False,
        "operational_pilot_only": True,
        "execution_contract_sha256": digest(out / "execution_contract.json"),
        "model_sha256": MODEL_SHA,
        "worker_sha256_text_lf": execution["worker_sha256_text_lf"],
        "numerical_file_sha256": {path.name: digest(path) for path in sorted(out.iterdir()) if path.is_file()},
        "bbox_uncertainty_endpoints": bbox_sensitive,
        "schedule": schedule_report,
        "transfer": transfer,
        "endpoint_usable_heads": endpoint_usable,
        "primary_measurement_eligible_heads": primary_measurement_eligible,
        "all_endpoint_usable_heads": all_endpoint_usable,
        "bbox_uncertainty_coverage": bbox_coverage,
        "paired_context_support": context_support,
        "elapsed_seconds": elapsed,
        "mean_seconds_per_scheduled_photo": elapsed / len(schedule) if schedule else None,
        "mean_original_bytes_per_downloaded_photo": transfer["bytes"] / transfer["download_success"] if transfer["download_success"] else None,
        "maximum_source_image_bytes_live_at_once": max_live_bytes,
        "source_images_persisted": 0,
        "environment_values_read": 0,
        "ecological_models_executed": 0,
        "technical_uncertainty": {
            "bbox": "per-head baseline and all four +/-5% shift values/statuses, recipes, maximum and median absolute changes are retained for orientation and gross-shape endpoints",
            "colour": "paired flower/non-head/green-context diagnostics are retained per head; population photometric stress sensitivity is supplied by the frozen 14-condition technical audit rather than re-running six synthetic photometric transforms for every original image"
        },
        "limits": [
            "The pilot measures execution throughput and endpoint yield, not ecological associations.",
            "Production is blocked until reconciled source/dependence scheduling and a verified durable private numerical destination are implemented.",
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
