"""Blinded, native-only Cirsium head anatomy review packets.

Input data are private, not committed. Original bytes and WCVP native membership
are checked independently of observed flower colour or detector success. A/B
packets have different opaque tokens and order; exact source IDs and attribution
are in the PRIVATE trusted codebook, never in the public repository.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps

SOURCE_SHA = "5600dbe9db39dddbf0f74a5ba9534d2abfd865be0c3071f781c1c0a4efc7aece"
REVIEW_STATES = {"unreviewed", "completed", "unscorable"}
FAILURE_REASONS = {
    "no_head_visible", "unresolved_head", "head_obscured", "severe_blur",
    "too_small", "not_flowering", "insufficient_colour", "unresolved_tissue", "other",
}
WHORLS = {"outer", "middle", "inner", "unknown"}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def native_roster(path):
    data = pd.read_csv(
        path,
        usecols=lambda column: column in {"obs_id", "native_primary", "accepted_key"},
        dtype={"obs_id": str, "accepted_key": str},
        low_memory=False,
    )
    if not {"obs_id", "native_primary", "accepted_key"} <= set(data):
        raise ValueError("Native eligibility must be explicitly supplied")
    if data.obs_id.isna().any() or data.obs_id.duplicated().any():
        raise ValueError("Source observation identifiers must be unique and nonmissing")
    if data.native_primary.isna().any():
        raise ValueError("Unknown source-native status")
    if data.native_primary.dtype != bool:
        status = data.native_primary.astype(str).str.lower()
        if not status.isin({"true", "false"}).all():
            raise ValueError("Unrecognized source-native status")
        data["native_primary"] = status.eq("true")
    return data.loc[data.native_primary, ["obs_id", "accepted_key"]].dropna().set_index("obs_id")


def source_audit(archive):
    records = json.loads(archive.read("acquisition_ledger.json"))
    if not 1 <= len(records) <= 100:
        raise ValueError("Process bounded batches of 1–100 original photographs")
    if len({str(x["obs_id"]) for x in records}) != len(records):
        raise ValueError("Duplicate source observation")
    if len({str(x["photo_id"]) for x in records}) != len(records):
        raise ValueError("Duplicate source photograph")
    return records


def blank_review(token, filename):
    return {
        "image_token": token, "image_file": "images/" + filename,
        "review_status": "unreviewed", "failure_reason": None,
        "flowering_stage": None, "view_pose": None,
        "head_base": None, "floral_end": None, "peduncle_point": None,
        "involucre_rim": None, "floral_polygons": None,
        "involucre_polygon": None, "anatomical_phyllaries": None,
        "visibility_and_colour_confidence": None, "reviewer_comments": None,
    }


def build(source_zip, native_csv, outdir, source_sha256=SOURCE_SHA):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    raw = Path(source_zip).read_bytes()
    if sha(raw) != source_sha256:
        raise ValueError("Original-source ZIP checksum mismatch")
    native = native_roster(native_csv)
    codebook = []
    entries = {"A": [], "B": []}
    pictures = {"A": {}, "B": {}}
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        records = source_audit(archive)
        native_records = [x for x in records if str(x["obs_id"]) in native.index]
        for x in native_records:
            if x["status"] != "downloaded_original":
                continue
            blob = archive.read(x["original_path"])
            if sha(blob) != x["sha256"]:
                raise ValueError("Original image byte hash mismatch")
            image = ImageOps.exif_transpose(Image.open(io.BytesIO(blob))).convert("RGB")
            image.thumbnail((1400, 1400))
            encoded = io.BytesIO()
            image.save(encoded, "JPEG", quality=88)
            obs_id = str(x["obs_id"])
            record = {
                "original_photo_id": str(x["photo_id"]),
                "original_observation_id": obs_id,
                "original_sha256": x["sha256"],
                "accepted_key": str(native.loc[obs_id, "accepted_key"]),
                "source_license": x.get("license"),
                "source_attribution": x.get("attribution"),
                "source_width": x.get("width"),
                "source_height": x.get("height"),
                "preview_width": image.width,
                "preview_height": image.height,
            }
            for person in ("A", "B"):
                token = person + "-" + sha(
                    (person + "|" + sha(x["sha256"].encode())).encode()
                )[:12]
                filename = token + ".jpg"
                entries[person].append(blank_review(token, filename))
                pictures[person][filename] = encoded.getvalue()
                record["token_" + person] = token
            codebook.append(record)
    readme = (
        "INDEPENDENT HEAD ANATOMY REVIEW: do not share A/B annotations.\n"
        "Each photo can be nonflowering, blurred, without heads, or unscorable.\n"
        "Use normalized [x,y] coordinates, 0-1, against the supplied preview.\n"
        "Set review_status=completed or unscorable; give failure_reason if unscorable.\n"
        "Do not equate missing phyllaries (null) with truly absent phyllaries ([]).\n"
        "Annotate floral tissue ROI independent of pixel hue; retain white parts.\n"
        "For visible phyllaries, use base/tip landmarks and whorl (unknown allowed).\n"
        "Do not assume photographed image vertical is gravity, or shine is sticky.\n"
        "Same-head mask and landmark provenance must be retained.\n"
        "Colour pixels will be remeasured from verified ORIGINALS, not JPEG previews.\n"
    )
    for person in ("A", "B"):
        entries[person].sort(key=lambda row: sha(
            (person + "|review-order|" + row["image_token"]).encode()
        ))
        with zipfile.ZipFile(
            out / ("blinded_review_" + person + ".zip"), "w", zipfile.ZIP_DEFLATED
        ) as archive:
            for filename, blob in pictures[person].items():
                archive.writestr("images/" + filename, blob)
            archive.writestr(
                "annotations_to_complete.json",
                json.dumps(entries[person], ensure_ascii=False, indent=2),
            )
            archive.writestr("READ_ME.txt", readme)
    codebook_path = out / "PRIVATE_codebook.json"
    codebook_path.write_text(json.dumps(codebook, ensure_ascii=False, indent=2))
    receipt = {
        "source_photos": len(records),
        "source_native_observations": len(native_records),
        "native_original_photos_recovered": len(codebook),
        "native_download_failures": len(native_records) - len(codebook),
        "native_accepted_taxa": len({x["accepted_key"] for x in codebook}),
        "packet_entries_per_rater": {person: len(entries[person]) for person in entries},
        "independent_labels_completed": 0,
        "independent_reference_validated": False,
        "all_native_originals_included_before_trait_filter": True,
        "native_join_key": "observation_id_not_photo_id",
        "original_source_sha256": source_sha256,
        "private_codebook_sha256": sha(codebook_path.read_bytes()),
        "private_image_and_source_ids_published": False,
        "ecological_inference_ready": False,
    }
    (out / "AGGREGATE_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def polygon_mask(polygons, shape):
    if polygons is None:
        return None
    height, width = shape
    canvas = Image.new("1", (width, height), 0)
    artist = ImageDraw.Draw(canvas)
    for polygon in polygons:
        if len(polygon) < 3:
            raise ValueError("Polygon needs at least three landmarks")
        for xx, yy in polygon:
            if not all(np.isfinite([xx, yy])) or not (0 <= xx <= 1 and 0 <= yy <= 1):
                raise ValueError("Polygon coordinates out of bounds")
        artist.polygon([(round(xx * (width-1)), round(yy * (height-1)))
                        for xx, yy in polygon], fill=1)
    return np.asarray(canvas, dtype=bool)


def mask_iou(first, second):
    if first is None or second is None:
        return None
    union = int(np.logical_or(first, second).sum())
    return float(np.logical_and(first, second).sum() / union) if union else None


def validate_review(rows, expected):
    if len(rows) != len(expected) or {r.get("image_token") for r in rows} != set(expected):
        raise ValueError("Incomplete or mismatched reviewer token set")
    if len(set(r.get("image_token") for r in rows)) != len(rows):
        raise ValueError("Duplicate rater token")
    for row in rows:
        state = row.get("review_status")
        if state not in REVIEW_STATES:
            raise ValueError("Invalid review state")
        if state == "unscorable" and row.get("failure_reason") not in FAILURE_REASONS:
            raise ValueError("Unscorable image requires explicit reason")
        if state == "completed" and row.get("failure_reason") is not None:
            raise ValueError("Completed review cannot carry failure reason")
        for field in ("head_base", "floral_end", "peduncle_point", "involucre_rim"):
            point = row.get(field)
            if point is not None and (
                len(point) != 2 or not np.isfinite(point).all()
                or any(v < 0 or v > 1 for v in point)
            ):
                raise ValueError("Out-of-range landmark " + field)
        members = row.get("anatomical_phyllaries")
        if members is not None:
            for item in members:
                if item.get("whorl", "unknown") not in WHORLS:
                    raise ValueError("Incorrect phyllary whorl")
                for field in ("base", "tip"):
                    point = item.get(field)
                    if point is None or len(point) != 2 or not np.isfinite(point).all() or not all(
                        0 <= v <= 1 for v in point
                    ):
                        raise ValueError("Invalid phyllary landmark")


def evaluate(codebook, review_a, review_b):
    book = json.loads(Path(codebook).read_text())
    a = json.loads(Path(review_a).read_text())
    b = json.loads(Path(review_b).read_text())
    validate_review(a, [item["token_A"] for item in book])
    validate_review(b, [item["token_B"] for item in book])
    a = {item["image_token"]: item for item in a}
    b = {item["image_token"]: item for item in b}
    both = 0
    stage_matches, contour_scores, point_errors, count_diffs = [], [], [], []
    for identity in book:
        left = a[identity["token_A"]]
        right = b[identity["token_B"]]
        if left["review_status"] != "unreviewed" and right["review_status"] != "unreviewed":
            both += 1
        if left["review_status"] != "completed" or right["review_status"] != "completed":
            continue
        if left.get("flowering_stage") is not None and right.get("flowering_stage") is not None:
            stage_matches.append(left["flowering_stage"] == right["flowering_stage"])
        for field in ("floral_polygons", "involucre_polygon"):
            p = left.get(field)
            q = right.get(field)
            if field == "involucre_polygon":
                p = [p] if p else None
                q = [q] if q else None
            score = mask_iou(polygon_mask(p, (256, 256)), polygon_mask(q, (256, 256)))
            if score is not None:
                contour_scores.append((field, score))
        for field in ("head_base", "floral_end", "peduncle_point", "involucre_rim"):
            x, y = left.get(field), right.get(field)
            if x is not None and y is not None:
                point_errors.append((field, float(np.linalg.norm(np.asarray(x)-y))))
        x = left.get("anatomical_phyllaries")
        y = right.get("anatomical_phyllaries")
        if x is not None and y is not None:
            count_diffs.append(abs(len(x) - len(y)))
    result = {
        "source_photo_records": len(book),
        "both_raters_attempted": both,
        "both_raters_attempted_fraction": both/len(book) if book else 0,
        "stage_agreement_fraction": float(np.mean(stage_matches)) if stage_matches else None,
        "n_paired_stages": len(stage_matches),
        "median_polygon_iou": {
            field: float(np.median([v for k,v in contour_scores if k==field]))
            if any(k==field for k,_ in contour_scores) else None
            for field in ("floral_polygons", "involucre_polygon")
        },
        "median_normalized_landmark_disagreement": {
            field: float(np.median([v for k,v in point_errors if k==field]))
            if any(k==field for k,_ in point_errors) else None
            for field in ("head_base", "floral_end", "peduncle_point", "involucre_rim")
        },
        "median_disagreement_in_annotated_visible_phyllary_count":
            float(np.median(count_diffs)) if count_diffs else None,
        "independent_anatomical_validation_complete": False,
        "assessment": "PENDING_INDEPENDENT_REVIEW" if both==0
                      else "RATERS_COMPLETE_PARTLY_OR_FULLY_ADJUDICATION_PENDING",
    }
    return result


def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="command",required=True)
    a=sub.add_parser("prepare")
    a.add_argument("--source-zip",required=True)
    a.add_argument("--source-sha256",default=SOURCE_SHA)
    a.add_argument("--native-ledger",required=True)
    a.add_argument("--out-dir",required=True)
    a=sub.add_parser("evaluate")
    a.add_argument("--private-codebook",required=True)
    a.add_argument("--review-a",required=True)
    a.add_argument("--review-b",required=True)
    a.add_argument("--out-file",required=True)
    args=p.parse_args()
    if args.command=="prepare":
        result=build(args.source_zip,args.native_ledger,args.out_dir,args.source_sha256)
    else:
        result=evaluate(args.private_codebook,args.review_a,args.review_b)
        Path(args.out_file).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
