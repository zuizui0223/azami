"""Bounded medium-vs-original image-resolution gate for Chapter 1 v3.

This pilot never retains source image bytes after a photo pair has been measured.
It samples rows by a deterministic hash before image retrieval, downloads one
medium/original pair at a time, runs the same pinned detector and 27-endpoint
feature engine on both versions, matches detections in normalized image space,
and retains only manifests, hashes, dimensions and numerical measurements.

The gate asks whether native-resolution processing is worth a larger streamed
pass. It does not establish physical trait accuracy and does not require or
justify a full original-image archive.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
from pathlib import Path
import re
from typing import Iterable

LICENSES = {
    "cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa",
    "cc-by-nd", "cc-by-nc-nd",
}
SELECTION_SALT = "ch1-v3-resolution-stream-gate-v1"
SIZE_RE = re.compile(r"/(square|small|medium|large|original)(\.[^/?#]+)(?=\?|#|$)", re.I)


def text(value) -> str:
    return "" if value is None else str(value).strip()


def sized_url(url: str, size: str) -> str:
    if size not in {"square", "small", "medium", "large", "original"}:
        raise ValueError("Unsupported image size")
    url = text(url)
    if not url:
        return ""
    replaced, n = SIZE_RE.subn(lambda m: f"/{size}{m.group(2)}", url, count=1)
    if n != 1:
        raise ValueError("Image URL has no recognized iNaturalist size component")
    return replaced


def score_row(row: dict[str, str]) -> int:
    payload = f"{SELECTION_SALT}|{text(row.get('obs_id'))}|{text(row.get('photo_id'))}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest(), "big")


def eligible_source_row(row: dict[str, str]) -> bool:
    if not text(row.get("obs_id")) or not text(row.get("photo_id")):
        return False
    if text(row.get("photo_license_code")).lower() not in LICENSES:
        return False
    base = text(row.get("medium_image_url")) or text(row.get("large_image_url")) or text(row.get("raw_image_url"))
    try:
        return bool(base and sized_url(base, "medium") and sized_url(base, "original"))
    except ValueError:
        return False


def select_rows(rows: Iterable[dict[str, str]], n: int) -> list[dict[str, str]]:
    """Return the n lowest hash-scored first eligible photos, one per observation."""
    if n <= 0:
        raise ValueError("Pilot size must be positive")
    seen_obs: set[str] = set()
    heap: list[tuple[int, str, dict[str, str]]] = []
    for row in rows:
        if not eligible_source_row(row):
            continue
        obs = text(row["obs_id"])
        if obs in seen_obs:
            continue
        seen_obs.add(obs)
        score = score_row(row)
        item = (-score, text(row.get("photo_id")), dict(row))
        if len(heap) < n:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
    chosen = sorted([(-neg, row) for neg, _, row in heap], key=lambda item: (item[0], text(item[1].get("photo_id"))))
    result = []
    for score, row in chosen:
        base = text(row.get("medium_image_url")) or text(row.get("large_image_url")) or text(row.get("raw_image_url"))
        row = dict(row)
        row["selection_score_sha256"] = f"{score:064x}"
        row["medium_url"] = sized_url(base, "medium")
        row["original_url"] = sized_url(base, "original")
        result.append(row)
    return result


def norm_box(box, width: int, height: int) -> tuple[float, float, float, float]:
    if width <= 0 or height <= 0 or len(box) != 4:
        raise ValueError("Invalid box or image extent")
    x1, y1, x2, y2 = map(float, box)
    if not all(math.isfinite(v) for v in (x1, y1, x2, y2)) or x2 <= x1 or y2 <= y1:
        raise ValueError("Invalid detector box")
    return x1 / width, y1 / height, x2 / width, y2 / height


def box_iou(a, b) -> float:
    ax1, ay1, ax2, ay2 = a
    bx1, by1, bx2, by2 = b
    ix1, iy1, ix2, iy2 = max(ax1, bx1), max(ay1, by1), min(ax2, bx2), min(ay2, by2)
    intersection = max(0.0, ix2 - ix1) * max(0.0, iy2 - iy1)
    union = (ax2-ax1)*(ay2-ay1) + (bx2-bx1)*(by2-by1) - intersection
    return intersection / union if union > 0 else 0.0


def greedy_matches(medium_boxes, original_boxes, medium_size, original_size, min_iou: float = 0.50):
    if not 0 <= min_iou <= 1:
        raise ValueError("min_iou must be in [0,1]")
    mw, mh = medium_size
    ow, oh = original_size
    candidates = []
    for mi, m in enumerate(medium_boxes):
        mn = norm_box(m, mw, mh)
        for oi, o in enumerate(original_boxes):
            score = box_iou(mn, norm_box(o, ow, oh))
            if score >= min_iou:
                candidates.append((score, mi, oi))
    candidates.sort(key=lambda x: (-x[0], x[1], x[2]))
    used_m, used_o, matches = set(), set(), []
    for score, mi, oi in candidates:
        if mi in used_m or oi in used_o:
            continue
        used_m.add(mi)
        used_o.add(oi)
        matches.append((mi, oi, score))
    return matches


def percentile(values: list[float], q: float):
    finite = sorted(v for v in values if math.isfinite(v))
    if not finite:
        return None
    if len(finite) == 1:
        return finite[0]
    pos = (len(finite)-1) * q
    lo, hi = math.floor(pos), math.ceil(pos)
    if lo == hi:
        return finite[lo]
    return finite[lo] * (hi-pos) + finite[hi] * (pos-lo)


def rank_correlation(x: list[float], y: list[float]):
    if len(x) != len(y) or len(x) < 3:
        return None
    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: (values[i], i))
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i + 1
            while j < len(order) and values[order[j]] == values[order[i]]:
                j += 1
            rank = (i + 1 + j) / 2.0
            for k in range(i, j):
                out[order[k]] = rank
            i = j
        return out
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx)/len(rx), sum(ry)/len(ry)
    num = sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    dx = math.sqrt(sum((a-mx)**2 for a in rx))
    dy = math.sqrt(sum((b-my)**2 for b in ry))
    return num/(dx*dy) if dx and dy else None


def summarize_endpoint(rows: list[dict], min_common: int = 15,
                       gain_threshold: float = 0.10,
                       rank_threshold: float = 0.90,
                       normalized_delta_threshold: float = 0.15) -> dict:
    medium_usable = sum(r["medium_status"] == "usable" for r in rows)
    original_usable = sum(r["original_status"] == "usable" for r in rows)
    both_usable_rows = [r for r in rows if r["medium_status"] == r["original_status"] == "usable"
                        and r["medium_value"] is not None and r["original_value"] is not None]
    gains = sum(r["original_status"] == "usable" and r["medium_status"] != "usable" for r in rows)
    losses = sum(r["medium_status"] == "usable" and r["original_status"] != "usable" for r in rows)
    m = [float(r["medium_value"]) for r in both_usable_rows]
    o = [float(r["original_value"]) for r in both_usable_rows]
    deltas = [abs(a-b) for a,b in zip(m,o)]
    iqr = None
    normalized_median = None
    if len(o) >= 4:
        q25, q75 = percentile(o, .25), percentile(o, .75)
        iqr = q75-q25 if q25 is not None and q75 is not None else None
        if iqr and iqr > 0:
            normalized_median = percentile(deltas, .5) / iqr
    gain_fraction = gains / max(1, original_usable)
    rho = rank_correlation(m, o)
    enough = len(both_usable_rows) >= min_common
    flags = []
    if enough and gain_fraction >= gain_threshold:
        flags.append("original_eligibility_gain")
    if enough and rho is not None and rho < rank_threshold:
        flags.append("rank_disagreement")
    if enough and normalized_median is not None and normalized_median >= normalized_delta_threshold:
        flags.append("large_native_unit_shift_relative_to_original_iqr")
    status = "INSUFFICIENT_COMMON_USABLE" if not enough else "ORIGINAL_MAY_ADD_INFORMATION" if flags else "NO_PILOT_EVIDENCE_ORIGINAL_REQUIRED"
    return {
        "rows": len(rows), "medium_usable": medium_usable, "original_usable": original_usable,
        "both_usable": len(both_usable_rows), "original_gain": gains, "original_loss": losses,
        "original_gain_fraction_of_original_usable": gain_fraction,
        "spearman_like_rank_correlation": rho,
        "median_abs_difference": percentile(deltas, .5), "p95_abs_difference": percentile(deltas, .95),
        "original_iqr": iqr, "median_abs_difference_over_original_iqr": normalized_median,
        "status": status, "flags": flags,
    }


def _sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _download(session, url: str, timeout: int = 90):
    response = session.get(url, timeout=(15, timeout))
    response.raise_for_status()
    data = response.content
    if len(data) < 1000:
        raise ValueError("Downloaded image payload is too small")
    return data


def _decode(data: bytes):
    import io
    import numpy as np
    from PIL import Image, ImageOps
    with Image.open(io.BytesIO(data)) as opened:
        image = ImageOps.exif_transpose(opened).convert("RGB")
        image.load()
    return image, np.asarray(image)[:, :, ::-1].copy()


def _detector_records(model, bgr):
    from .detect_cached_images import PARAMETERS
    predictions = model.predict(bgr, **PARAMETERS)
    if len(predictions) != 1:
        raise ValueError("Detector returned an unexpected image count")
    boxes = predictions[0].boxes
    records = [] if boxes is None else [(list(map(float, b)), float(c), int(k)) for b,c,k in
        zip(boxes.xyxy.cpu().tolist(), boxes.conf.cpu().tolist(), boxes.cls.cpu().tolist())]
    records = [r for r in records if r[2] == 0]
    records.sort(key=lambda r: (-r[1], r[2], *r[0]))
    return records


def _crop_measure(bgr, records, index):
    from .detect_cached_images import crop_box
    from . import image_features as features
    h, w = bgr.shape[:2]
    padded_heads = [crop_box(r[0], w, h, .12)[0] for r in records]
    head_box = padded_heads[index]
    context_box = crop_box(records[index][0], w, h, .8)[0]
    x1,y1,x2,y2 = head_box
    cx1,cy1,cx2,cy2 = context_box
    head = bgr[y1:y2, x1:x2].copy()
    context = bgr[cy1:cy2, cx1:cx2].copy()
    result, _ = features.measure(head, context, context_box, padded_heads)
    return result


def run(metadata: Path, weights: Path, out: Path, pilot_photos: int = 32, min_iou: float = .50):
    import requests
    from ultralytics import YOLO
    from .detect_cached_images import MODEL_SHA
    from .workflow import digest
    if digest(weights) != MODEL_SHA:
        raise ValueError("Detector weight identity mismatch")
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Use a fresh output directory")
    with metadata.open(encoding="utf-8-sig", newline="") as handle:
        selected = select_rows(csv.DictReader(handle), pilot_photos)
    private_fields = ["pair_id","obs_id","photo_id","photo_license_code","medium_url","original_url","selection_score_sha256"]
    private_rows = []
    for i,row in enumerate(selected):
        private_rows.append({"pair_id": f"pair_{i:04d}", **{k:text(row.get(k)) for k in private_fields if k != "pair_id"}})
    with (out/"pilot_manifest_private.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=private_fields)
        writer.writeheader()
        writer.writerows(private_rows)
    model = YOLO(str(weights))
    session = requests.Session()
    session.headers.update({"User-Agent": "azami-ch1-v3-resolution-stream-gate/1.0"})
    transfer_rows, endpoint_rows = [], []
    max_live_bytes = 0
    for pair, source in zip(private_rows, selected):
        pid = pair["pair_id"]
        medium_data = original_data = None
        try:
            medium_data = _download(session, source["medium_url"])
            original_data = _download(session, source["original_url"])
            max_live_bytes = max(max_live_bytes, len(medium_data)+len(original_data))
            medium_pil, medium_bgr = _decode(medium_data)
            original_pil, original_bgr = _decode(original_data)
            mr, orr = _detector_records(model, medium_bgr), _detector_records(model, original_bgr)
            matches = greedy_matches([r[0] for r in mr], [r[0] for r in orr], medium_pil.size, original_pil.size, min_iou)
            transfer_rows.append({"pair_id":pid,"status":"success","medium_bytes":len(medium_data),"original_bytes":len(original_data),
                "medium_sha256":_sha_bytes(medium_data),"original_sha256":_sha_bytes(original_data),
                "medium_width":medium_pil.width,"medium_height":medium_pil.height,"original_width":original_pil.width,"original_height":original_pil.height,
                "medium_detections":len(mr),"original_detections":len(orr),"matched_heads":len(matches),"error":""})
            for match_index,(mi,oi,iou) in enumerate(matches):
                mm, oo = _crop_measure(medium_bgr,mr,mi), _crop_measure(original_bgr,orr,oi)
                m_by = {e["endpoint_id"]:e for e in mm["endpoints"]}
                o_by = {e["endpoint_id"]:e for e in oo["endpoints"]}
                for endpoint_id in sorted(m_by):
                    me, oe = m_by[endpoint_id], o_by[endpoint_id]
                    endpoint_rows.append({"pair_id":pid,"match_index":match_index,"normalized_box_iou":iou,"endpoint_id":endpoint_id,"unit":me["unit"],
                        "medium_value":me["value"],"medium_status":me["status"],"original_value":oe["value"],"original_status":oe["status"]})
        except Exception as error:
            transfer_rows.append({"pair_id":pid,"status":"error","medium_bytes":0,"original_bytes":0,"medium_sha256":"","original_sha256":"",
                "medium_width":"","medium_height":"","original_width":"","original_height":"","medium_detections":"","original_detections":"","matched_heads":0,
                "error":type(error).__name__+": "+str(error)})
        finally:
            medium_data = original_data = None
    tfields = list(transfer_rows[0]) if transfer_rows else ["pair_id","status"]
    with (out/"pair_transfer_and_detection.csv").open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=tfields)
        writer.writeheader()
        writer.writerows(transfer_rows)
    efields=["pair_id","match_index","normalized_box_iou","endpoint_id","unit","medium_value","medium_status","original_value","original_status"]
    with (out/"matched_endpoint_measurements.csv").open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=efields)
        writer.writeheader()
        writer.writerows(endpoint_rows)
    by_endpoint={}
    for row in endpoint_rows:
        by_endpoint.setdefault(row["endpoint_id"],[]).append(row)
    summary={ep:summarize_endpoint(rows) for ep,rows in sorted(by_endpoint.items())}
    total_transferred=sum(int(r.get("medium_bytes") or 0)+int(r.get("original_bytes") or 0) for r in transfer_rows)
    report={"status":"RESOLUTION_STREAM_GATE_COMPLETE","pilot_photos_requested":pilot_photos,"pilot_photos_selected":len(selected),
        "successful_pairs":sum(r["status"]=="success" for r in transfer_rows),"matched_head_pairs":sum(int(r.get("matched_heads") or 0) for r in transfer_rows),
        "bytes_transferred_total":total_transferred,"maximum_source_image_bytes_live_at_once":max_live_bytes,
        "source_images_persisted":0,"full_original_archive_required_by_this_gate":False,
        "endpoint_summary":summary,"ecological_models_executed":False,
        "limits":["This is a deterministic technical pilot, not a representative ecological cohort.",
                  "A medium/original difference combines server-side resizing, JPEG encoding, detector localization and feature extraction.",
                  "No source image bytes are retained; hashes and numerical outputs do not establish physical accuracy.",
                  "A flagged endpoint justifies a larger streamed native-resolution test, not a 1.5-TB archive."]}
    (out/"resolution_gate_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata",type=Path,required=True)
    parser.add_argument("--weights",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--pilot-photos",type=int,default=32)
    parser.add_argument("--min-iou",type=float,default=.50)
    args=parser.parse_args()
    print(json.dumps(run(args.metadata,args.weights,args.out_dir,args.pilot_photos,args.min_iou),indent=2))

if __name__ == "__main__":
    main()
