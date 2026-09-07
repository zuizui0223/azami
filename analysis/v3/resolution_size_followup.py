"""Post-inspection streamed image-size follow-up for Chapter 1 v3.

The completed 32-photo medium-vs-original pilot showed much higher endpoint
eligibility at original resolution but did not reach its predeclared minimum of
15 common-usable head pairs for any endpoint. This follow-up therefore asks a
narrower technical question before any full-source image execution: does the
standard iNaturalist ``large`` derivative retain the usable information seen in
``original``?

Selection is fixed by a new hash salt and never uses taxon, geography,
environment or trait values. Source image bytes are held for one pair only and
are never written to the output directory. This is post-inspection technical
follow-up, not preregistration or independent accuracy validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
from pathlib import Path

from .resolution_stream_gate import (
    LICENSES,
    _crop_measure,
    _decode,
    _detector_records,
    _download,
    _sha_bytes,
    greedy_matches,
    percentile,
    rank_correlation,
    sized_url,
    text,
)

FOLLOWUP_SALT = "ch1-v3-large-original-followup-v1"


def eligible(row: dict[str, str]) -> bool:
    if not text(row.get("obs_id")) or not text(row.get("photo_id")):
        return False
    if text(row.get("photo_license_code")).lower() not in LICENSES:
        return False
    base = text(row.get("large_image_url")) or text(row.get("medium_image_url")) or text(row.get("raw_image_url"))
    if not base:
        return False
    try:
        return bool(sized_url(base, "large") and sized_url(base, "original"))
    except ValueError:
        return False


def selection_score(row: dict[str, str]) -> int:
    payload = f"{FOLLOWUP_SALT}|{text(row.get('obs_id'))}|{text(row.get('photo_id'))}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest(), "big")


def select_rows(rows, n: int) -> list[dict[str, str]]:
    if n <= 0:
        raise ValueError("Follow-up size must be positive")
    seen_obs: set[str] = set()
    heap: list[tuple[int, str, dict[str, str]]] = []
    for row in rows:
        if not eligible(row):
            continue
        obs = text(row["obs_id"])
        if obs in seen_obs:
            continue
        seen_obs.add(obs)
        score = selection_score(row)
        item = (-score, text(row.get("photo_id")), dict(row))
        if len(heap) < n:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
    chosen = sorted([(-neg, row) for neg, _, row in heap], key=lambda x: (x[0], text(x[1].get("photo_id"))))
    result = []
    for score, row in chosen:
        base = text(row.get("large_image_url")) or text(row.get("medium_image_url")) or text(row.get("raw_image_url"))
        row = dict(row)
        row["selection_score_sha256"] = f"{score:064x}"
        row["baseline_url"] = sized_url(base, "large")
        row["comparison_url"] = sized_url(base, "original")
        result.append(row)
    return result


def exact_two_sided_binomial_half(successes: int, failures: int) -> float | None:
    """Exact two-sided p-value for discordant paired states under p=0.5."""
    n = successes + failures
    if n <= 0:
        return None
    tail = min(successes, failures)
    probability = sum(math.comb(n, k) for k in range(tail + 1)) / (2 ** n)
    return min(1.0, 2.0 * probability)


def summarize_endpoint(rows: list[dict], min_common: int = 30,
                       min_discordant: int = 15, comparison_gain_share: float = .80,
                       eligibility_p_max: float = .01, rank_threshold: float = .90,
                       normalized_delta_threshold: float = .15) -> dict:
    baseline_usable = sum(r["baseline_status"] == "usable" for r in rows)
    comparison_usable = sum(r["comparison_status"] == "usable" for r in rows)
    both = [r for r in rows if r["baseline_status"] == r["comparison_status"] == "usable"
            and r["baseline_value"] is not None and r["comparison_value"] is not None]
    gains = sum(r["comparison_status"] == "usable" and r["baseline_status"] != "usable" for r in rows)
    losses = sum(r["baseline_status"] == "usable" and r["comparison_status"] != "usable" for r in rows)
    discordant = gains + losses
    gain_share = gains / discordant if discordant else None
    eligibility_p = exact_two_sided_binomial_half(gains, losses)

    b = [float(r["baseline_value"]) for r in both]
    c = [float(r["comparison_value"]) for r in both]
    deltas = [abs(a-d) for a, d in zip(b, c)]
    rho = rank_correlation(b, c)
    iqr = normalized = None
    if len(c) >= 4:
        q25, q75 = percentile(c, .25), percentile(c, .75)
        iqr = q75 - q25 if q25 is not None and q75 is not None else None
        if iqr and iqr > 0:
            normalized = percentile(deltas, .5) / iqr

    flags = []
    eligibility_decidable = discordant >= min_discordant
    numeric_decidable = len(both) >= min_common
    if (eligibility_decidable and gain_share is not None and gain_share >= comparison_gain_share
            and eligibility_p is not None and eligibility_p <= eligibility_p_max):
        flags.append("original_eligibility_gain_over_large")
    if numeric_decidable and rho is not None and rho < rank_threshold:
        flags.append("large_original_rank_disagreement")
    if numeric_decidable and normalized is not None and normalized >= normalized_delta_threshold:
        flags.append("large_original_native_unit_shift")

    if flags:
        status = "ORIGINAL_MAY_ADD_INFORMATION_BEYOND_LARGE"
    elif numeric_decidable or eligibility_decidable:
        status = "NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE"
    else:
        status = "INSUFFICIENT_FOLLOWUP_INFORMATION"
    return {
        "rows": len(rows),
        "large_usable": baseline_usable,
        "original_usable": comparison_usable,
        "both_usable": len(both),
        "original_gain": gains,
        "original_loss": losses,
        "discordant_eligibility_pairs": discordant,
        "original_gain_share_among_discordant": gain_share,
        "exact_two_sided_binomial_p": eligibility_p,
        "rank_correlation": rho,
        "median_abs_difference": percentile(deltas, .5),
        "p95_abs_difference": percentile(deltas, .95),
        "original_iqr": iqr,
        "median_abs_difference_over_original_iqr": normalized,
        "numeric_decidable": numeric_decidable,
        "eligibility_decidable": eligibility_decidable,
        "status": status,
        "flags": flags,
    }


def run(metadata: Path, weights: Path, out: Path, photos: int = 128, min_iou: float = .50):
    import requests
    from ultralytics import YOLO
    from .detect_cached_images import MODEL_SHA
    from .workflow import digest

    if digest(weights) != MODEL_SHA:
        raise ValueError("Detector weight identity mismatch")
    if photos <= 0 or not 0 <= min_iou <= 1:
        raise ValueError("Invalid follow-up size or IoU")
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Use a fresh output directory")

    with metadata.open(encoding="utf-8-sig", newline="") as handle:
        selected = select_rows(csv.DictReader(handle), photos)
    if len(selected) != photos:
        raise ValueError(f"Requested {photos} photos but only {len(selected)} eligible observations were selected")

    private_fields = ["pair_id", "obs_id", "photo_id", "photo_license_code", "baseline_url", "comparison_url", "selection_score_sha256"]
    private_rows = []
    for i, row in enumerate(selected):
        private_rows.append({"pair_id": f"followup_{i:04d}", **{k: text(row.get(k)) for k in private_fields if k != "pair_id"}})
    with (out / "pilot_manifest_private.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=private_fields)
        writer.writeheader()
        writer.writerows(private_rows)

    model = YOLO(str(weights))
    session = requests.Session()
    session.headers.update({"User-Agent": "azami-ch1-v3-large-original-followup/1.0"})
    transfer_rows, endpoint_rows = [], []
    max_live_bytes = 0
    for pair, source in zip(private_rows, selected):
        pid = pair["pair_id"]
        baseline_data = comparison_data = None
        try:
            baseline_data = _download(session, source["baseline_url"])
            comparison_data = _download(session, source["comparison_url"])
            max_live_bytes = max(max_live_bytes, len(baseline_data) + len(comparison_data))
            baseline_pil, baseline_bgr = _decode(baseline_data)
            comparison_pil, comparison_bgr = _decode(comparison_data)
            br = _detector_records(model, baseline_bgr)
            cr = _detector_records(model, comparison_bgr)
            matches = greedy_matches([r[0] for r in br], [r[0] for r in cr], baseline_pil.size, comparison_pil.size, min_iou)
            transfer_rows.append({
                "pair_id": pid, "status": "success",
                "large_bytes": len(baseline_data), "original_bytes": len(comparison_data),
                "large_sha256": _sha_bytes(baseline_data), "original_sha256": _sha_bytes(comparison_data),
                "large_width": baseline_pil.width, "large_height": baseline_pil.height,
                "original_width": comparison_pil.width, "original_height": comparison_pil.height,
                "large_detections": len(br), "original_detections": len(cr),
                "matched_heads": len(matches), "error": ""
            })
            for match_index, (bi, ci, iou) in enumerate(matches):
                bm = _crop_measure(baseline_bgr, br, bi)
                cm = _crop_measure(comparison_bgr, cr, ci)
                b_by = {e["endpoint_id"]: e for e in bm["endpoints"]}
                c_by = {e["endpoint_id"]: e for e in cm["endpoints"]}
                for endpoint_id in sorted(b_by):
                    be, ce = b_by[endpoint_id], c_by[endpoint_id]
                    endpoint_rows.append({
                        "pair_id": pid, "match_index": match_index, "normalized_box_iou": iou,
                        "endpoint_id": endpoint_id, "unit": be["unit"],
                        "baseline_size": "large", "comparison_size": "original",
                        "baseline_value": be["value"], "baseline_status": be["status"],
                        "comparison_value": ce["value"], "comparison_status": ce["status"]
                    })
        except Exception as error:
            transfer_rows.append({
                "pair_id": pid, "status": "error", "large_bytes": 0, "original_bytes": 0,
                "large_sha256": "", "original_sha256": "", "large_width": "", "large_height": "",
                "original_width": "", "original_height": "", "large_detections": "", "original_detections": "",
                "matched_heads": 0, "error": type(error).__name__ + ": " + str(error)
            })
        finally:
            baseline_data = comparison_data = None

    transfer_fields = list(transfer_rows[0]) if transfer_rows else ["pair_id", "status"]
    with (out / "pair_transfer_and_detection.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=transfer_fields)
        writer.writeheader()
        writer.writerows(transfer_rows)
    endpoint_fields = ["pair_id", "match_index", "normalized_box_iou", "endpoint_id", "unit",
                       "baseline_size", "comparison_size", "baseline_value", "baseline_status",
                       "comparison_value", "comparison_status"]
    with (out / "matched_endpoint_measurements.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=endpoint_fields)
        writer.writeheader()
        writer.writerows(endpoint_rows)

    by_endpoint: dict[str, list[dict]] = {}
    for row in endpoint_rows:
        by_endpoint.setdefault(row["endpoint_id"], []).append(row)
    endpoint_summary = {key: summarize_endpoint(rows) for key, rows in sorted(by_endpoint.items())}
    total_transferred = sum(int(r.get("large_bytes") or 0) + int(r.get("original_bytes") or 0) for r in transfer_rows)
    status_counts: dict[str, int] = {}
    for value in endpoint_summary.values():
        status_counts[value["status"]] = status_counts.get(value["status"], 0) + 1
    report = {
        "status": "LARGE_ORIGINAL_STREAM_FOLLOWUP_COMPLETE",
        "study_status": "post_inspection_technical_followup_after_32_photo_gate",
        "photos_requested": photos,
        "photos_selected": len(selected),
        "successful_pairs": sum(r["status"] == "success" for r in transfer_rows),
        "matched_head_pairs": sum(int(r.get("matched_heads") or 0) for r in transfer_rows),
        "bytes_transferred_total": total_transferred,
        "maximum_source_image_bytes_live_at_once": max_live_bytes,
        "source_images_persisted": 0,
        "full_original_archive_required_by_design": False,
        "decision_thresholds": {
            "minimum_common_usable_head_pairs": 30,
            "minimum_discordant_eligibility_pairs": 15,
            "original_gain_share_among_discordant": .80,
            "eligibility_exact_p_max": .01,
            "rank_correlation_min": .90,
            "median_abs_delta_over_original_iqr_max": .15
        },
        "endpoint_status_counts": status_counts,
        "endpoint_summary": endpoint_summary,
        "ecological_models_executed": False,
        "limits": [
            "This follow-up was designed after inspecting the completed 32-photo medium/original technical gate.",
            "It compares pipeline information retention between iNaturalist large and original server image versions, not physical botanical accuracy.",
            "Eligibility changes can reflect detector localization, crop pixel dimensions, JPEG resizing and endpoint QC thresholds.",
            "No source image bytes are retained, and no outcome can require a permanent full-source original-image archive.",
            "A resolution-sensitive endpoint can motivate streamed original processing for that endpoint; it does not validate the endpoint biologically."
        ]
    }
    (out / "resolution_followup_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--photos", type=int, default=128)
    parser.add_argument("--min-iou", type=float, default=.50)
    args = parser.parse_args()
    print(json.dumps(run(args.metadata, args.weights, args.out_dir, args.photos, args.min_iou), indent=2))


if __name__ == "__main__":
    main()
