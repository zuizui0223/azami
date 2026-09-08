"""Reopen a completed original-image pilot; verify all 27 slots without ecology.

No image requests, environmental data, v2 results or admission changes. File
identity, row conservation and bbox arithmetic are checked independently of the
worker's counters. This does not validate botanical accuracy or qualify a route.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import statistics

from .workflow import ROOT, digest


FILES = {
    "execution_contract.json", "selected_observations_private.csv",
    "endpoint_measurements_private.csv", "bbox_uncertainty_private.csv",
    "bbox_measurements_private.csv", "photo_observation_links_private.csv",
    "transfer_private.csv", "head_diagnostics_private.jsonl",
    "photo_detection_private.jsonl",
}
CONDITIONS = {"baseline", "bbox_left_5pct", "bbox_right_5pct", "bbox_up_5pct", "bbox_down_5pct"}


def read_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value):
    if value in (None, ""):
        return None
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError("Non-finite numeric serialization")
    return parsed


def require(condition, message):
    if not condition:
        raise ValueError(message)


def unique(rows, fields):
    result = {tuple(str(row[field]) for field in fields): row for row in rows}
    require(len(result) == len(rows), "Duplicate row identity: " + ",".join(fields))
    return result


def verify(directory: Path, decision_path: Path) -> dict:
    report_path = directory / "original_stream_report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    is_unit = report["status"] == "ORIGINAL_STREAM_PHOTO_UNIT_COMPLETE_NO_ECOLOGICAL_MODEL"
    require(is_unit or report["status"] == "ORIGINAL_STREAM_PILOT_COMPLETE_NO_ECOLOGICAL_MODEL", "Not a completed pilot or photo unit")
    require(set(report["numerical_file_sha256"]) == FILES, "Incomplete numerical-file manifest")
    for name in sorted(FILES):
        require(digest(directory / name) == report["numerical_file_sha256"][name], "Saved file hash changed: " + name)
    execution = json.loads((directory / "execution_contract.json").read_text(encoding="utf-8"))
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    require(digest(decision_path) == report["measurement_decision_sha256"] == execution["measurement_decision_sha256"], "Measurement decision identity changed")
    require(digest(directory / "execution_contract.json") == report["execution_contract_sha256"], "Execution identity changed")
    require(report["selection"] == execution["selection"], "Selection identity changed")
    require(report["production_execution_authorized"] is False and report["operational_pilot_only"] is (not is_unit), "Pilot scope was promoted")
    if is_unit:
        require(report.get("raw_measurement_unit_only") is True and execution["input_mode"] == "reconciled_photo_unit"
                and execution["status"] == "BOUNDED_RAW_PHOTO_UNIT_NO_ECOLOGY", "Photo-unit scope differs")
    require(report["source_images_persisted"] == report["environment_values_read"] == report["ecological_models_executed"] == 0, "Pilot scope was exceeded")
    routes = {row["endpoint_id"]: row["ecological_route"] for row in decision["endpoints"]}
    ids = set(routes)
    qualified = {key for key, value in routes.items() if value == "stream_original_required"}
    bbox_ids = set(report["bbox_uncertainty_endpoints"])
    expected_bbox = {row["endpoint_id"] for row in decision["endpoints"] if row.get("bbox_uncertainty_required") and row["endpoint_id"] in qualified}
    require(len(ids) == 27 and len(decision["endpoints"]) == 27 and len(qualified) == 14, "Endpoint denominator changed")
    require(set(report["retained_endpoints"]) == ids and set(report["qualified_endpoints"]) == qualified and bbox_ids == expected_bbox, "Reported endpoint sets changed")

    selected = read_rows(directory / "selected_observations_private.csv")
    selected_ids = {row["obs_id"] for row in selected}
    require(len(selected_ids) == len(selected) == report["schedule"]["selected_observations"], "Selected observation count differs")
    links = read_rows(directory / "photo_observation_links_private.csv")
    require({row["obs_id"] for row in links} == selected_ids, "Observation lost from link ledger")
    transfers = read_rows(directory / "transfer_private.csv")
    transfer_by = unique(transfers, ("photo_id",))
    scheduled_states = {"request_candidate_not_authorized", "scheduled"}
    scheduled = {row["photo_id"] for row in links if row["status"] in scheduled_states}
    require(set(transfer_by) == {(photo,) for photo in scheduled}, "Scheduled transfer missing or extra")
    if execution["input_mode"] in {"reconciled_whole_component_pilot", "reconciled_photo_unit"}:
        require(len(links) == report["schedule"]["selected_photo_links"], "Source links lost")
        require(len(transfers) == report["schedule"]["request_candidates"], "Request denominator differs")
        if is_unit:
            require(len(transfers) == 1, "Photo unit must contain exactly one request")
    detections = [json.loads(line) for line in (directory / "photo_detection_private.jsonl").read_text(encoding="utf-8").splitlines()]
    detection_by = unique(detections, ("photo_id",))
    require(set(detection_by) == set(transfer_by), "Detection/transfer identities differ")
    heads = set()
    for key, photo in detection_by.items():
        transfer = transfer_by[key]
        require(str(photo["transfer"]["status"]) == transfer["status"], "Transfer status differs")
        if transfer["detections"] != "":
            require(len(photo["detections"]) == int(transfer["detections"]), "Head count differs")
        for head in photo["detections"]:
            identity = (key[0], str(head["head_index"]))
            require(identity not in heads, "Duplicate detected head")
            heads.add(identity)
    require(len(heads) == report["transfer"]["heads"], "Reported total heads differs")
    require(sum(row["status"] == "success" for row in transfers) == report["transfer"]["success"], "Reported successes differ")
    require(sum(row["status"] != "success" for row in transfers) == report["transfer"]["error"], "Reported failures differ")

    endpoint_rows = read_rows(directory / "endpoint_measurements_private.csv")
    endpoints = unique(endpoint_rows, ("photo_id", "head_index", "endpoint_id"))
    require(set(endpoints) == {(*head, ep) for head in heads for ep in ids}, "Not exactly 27 raw endpoint slots per detected head")
    bbox = unique(read_rows(directory / "bbox_measurements_private.csv"), ("photo_id", "head_index", "endpoint_id", "condition"))
    uncertainty = unique(read_rows(directory / "bbox_uncertainty_private.csv"), ("photo_id", "head_index", "endpoint_id"))
    expected_uncertainty = {(*head, ep) for head in heads for ep in bbox_ids}
    require(set(uncertainty) == expected_uncertainty, "BBox uncertainty slots differ")
    require(set(bbox) == {(*key, condition) for key in expected_uncertainty for condition in CONDITIONS}, "BBox condition slots differ")
    stable, coverage = set(), {ep: 0 for ep in bbox_ids}
    for key, summary in uncertainty.items():
        base = endpoints[key]
        baseline = bbox[(*key, "baseline")]
        require(all(baseline[k] == base[k] for k in ("value", "original", "mirror", "status")), "BBox baseline differs from raw endpoint")
        changes = []
        for condition in sorted(CONDITIONS - {"baseline"}):
            moved = bbox[(*key, condition)]
            if base["status"] == moved["status"] == "usable":
                require(number(base["value"]) is not None and number(moved["value"]) is not None, "Usable bbox value is missing")
                changes.append(abs(number(base["value"]) - number(moved["value"])))
        require(int(summary["bbox_shift_usable_n"]) == len(changes), "BBox usable-shift count differs")
        for field, expected in (("bbox_shift_median_abs_change", statistics.median(changes) if changes else None), ("bbox_shift_max_abs_change", max(changes) if changes else None)):
            actual = number(summary[field])
            require(actual == expected or (actual is not None and expected is not None and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12)), "BBox arithmetic differs: " + field)
        if len(changes) == 4:
            stable.add(key)
            coverage[key[2]] += 1
    usable, eligible, statuses = Counter(), Counter(), defaultdict(Counter)
    finite_heads, usable_photos, eligible_photos = Counter(), defaultdict(set), defaultdict(set)
    for key, row in endpoints.items():
        ep = key[2]
        value = number(row["value"])
        require(row["ecological_route"] == routes[ep], "Saved route changed")
        is_usable = row["status"] == "usable"
        require(not is_usable or value is not None, "Usable endpoint is missing")
        admitted = is_usable and ep in qualified and (ep not in bbox_ids or key in stable)
        require(row["primary_measurement_eligible"] == str(admitted), "Measurement eligibility differs from saved decision")
        statuses[ep][row["status"]] += 1
        finite_heads[ep] += int(value is not None)
        usable[ep] += int(is_usable)
        eligible[ep] += int(admitted)
        if is_usable:
            usable_photos[ep].add(key[0])
        if admitted:
            eligible_photos[ep].add(key[0])
    require({ep: usable[ep] for ep in ids} == report["all_endpoint_usable_heads"], "All-27 usable counts differ")
    require({ep: usable[ep] for ep in qualified} == report["endpoint_usable_heads"], "Route usable counts differ")
    require({ep: eligible[ep] for ep in qualified} == report["primary_measurement_eligible_heads"], "Eligible counts differ")
    require({ep: {"heads": len(heads), "four_usable_shifts": coverage[ep]} for ep in bbox_ids} == report["bbox_uncertainty_coverage"], "BBox coverage differs")
    diagnostics = [json.loads(line) for line in (directory / "head_diagnostics_private.jsonl").read_text(encoding="utf-8").splitlines()]
    require(set(unique(diagnostics, ("photo_id", "head_index"))) == heads, "Head diagnostics lost or duplicated")
    photo_observations = defaultdict(set)
    for row in links:
        photo_observations[row["photo_id"]].add(row["obs_id"])
    def supported_observations(photos):
        return len({obs for photo in photos for obs in photo_observations[photo]})
    result = {
        "status": "ORIGINAL_STREAM_ALL27_NUMERICAL_INTEGRITY_VERIFIED",
        "input_report_sha256": digest(report_path),
        "verified_file_sha256": report["numerical_file_sha256"],
        "selected_observations": len(selected), "scheduled_photos": len(transfers),
        "transfer_status_counts": dict(Counter(row["status"] for row in transfers)),
        "detected_heads": len(heads), "raw_endpoint_slots": len(endpoints),
        "bbox_slots": len(bbox), "endpoint_status_counts": {ep: dict(statuses[ep]) for ep in sorted(ids)},
        "endpoint_coverage": {ep: {
            "route": routes[ep], "finite_raw_heads_including_qc_failure": finite_heads[ep],
            "raw_usable_heads": usable[ep], "raw_usable_photos": len(usable_photos[ep]),
            "raw_usable_observations": supported_observations(usable_photos[ep]),
            "measurement_eligible_photos": len(eligible_photos[ep]),
            "measurement_eligible_observations": supported_observations(eligible_photos[ep]),
        } for ep in sorted(ids)},
        "measurement_eligible_heads": {ep: eligible[ep] for ep in sorted(ids)},
        "technical_hold_endpoint_count": len(ids - qualified),
        "environment_values_read": 0, "v2_result_files_read": 0, "ecological_models_executed": 0,
        "ecological_fitting_authorized": False, "production_image_execution_authorized": False,
        "limits": ["Numerical integrity and operational coverage, not botanical accuracy or inference admission.",
                   "Missing or held endpoints retain their slots; this audit changes no thresholds or routes.",
                   "A completed pilot can contain failed transfers and unmeasurable heads.",
                   "Image hashes identify consumed bytes; source images were not saved or replayed here.",
                   "Observation coverage counts links; it does not aggregate traits or make shared photos independent."],
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--decision", type=Path, default=ROOT / "analysis/v3/measurement_qualification_decision_20260908.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), "Preserve previous verification output")
    report = verify(args.input_dir, args.decision)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({key: report[key] for key in ("status", "selected_observations", "scheduled_photos", "detected_heads", "raw_endpoint_slots", "bbox_slots", "transfer_status_counts")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
