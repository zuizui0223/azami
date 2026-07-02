r'''
Global Cirsium flowering-annotation + flower-only pipeline

Goal:
- Fetch global Cirsium observations from iNaturalist using id_above pagination.
- First filter observations by iNaturalist plant phenology annotation: Flowers and Fruits = Flowers.
- Inspect photos attached to those flowering observations.
- Download only a small temporary image first.
- Run YOLO immediately.
- Keep only photos where YOLO detects a Cirsium flower/head.
- Save flower-only images, YOLO crops, crop predictions, observation summaries, and maps.

Run:
    cd C:\Users\zuizui\cirsium_inat
    python global_cirsium_flowering_annotation_then_yolo_pipeline.py

Notes:
- This version does NOT stop at MAX_OBS by default. It is designed for all available records.
- To test first, set MAX_OBS = 5000 or MAX_PHOTOS = 10000.
- If interrupted, rerun the same command. It resumes from checkpoint files.
'''

from pathlib import Path
import time
import re
import json
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt

from ultralytics import YOLO
from tensorflow import keras


# ============================================================
# USER SETTINGS
# ============================================================
BASE = Path(r"C:\Users\zuizui\cirsium_inat")

# Full run: None. Test run: e.g. 5000 or 10000.
MAX_OBS = None
MAX_PHOTOS = None

PER_PAGE = 200
REQUEST_SLEEP = 0.8
DOWNLOAD_SLEEP = 0.10

GENUS_NAME = "Cirsium"

# iNaturalist controlled annotation filter:
# Flowers and Fruits = 12; Flowers = 13
# This prefilters observations to records annotated as flowering before downloading images.
USE_INAT_FLOWERING_ANNOTATION = True
FLOWERS_FRUITS_TERM_ID = 12
FLOWERS_TERM_VALUE_ID = 13

ONLY_GEOTAGGED = True
ONLY_PHOTOS = True
EXCLUDE_OBSCURED = False

# For screening, small is much faster. The kept image is also small by default.
# If you want nicer kept images later, set KEEP_IMAGE_SIZE = "medium".
SCREEN_IMAGE_SIZE = "small"
KEEP_IMAGE_SIZE = "small"

# If True, non-flower downloaded temp images are removed immediately.
DELETE_NON_FLOWER_TEMP = True

# YOLO flower/head detector
YOLO_CANDIDATES = [
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\weights\best.pt"),
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\best.pt"),
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\train\weights\best.pt"),
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\detect\train\weights\best.pt"),
]

# Trained upward/nodding crop classifier
HEAD_MODEL = BASE / "cirsium_head_direction_model_yolo_crops_fixed" / "cirsium_head_direction_yolo_crop_model.keras"

YOLO_CONF_THRES = 0.25
PAD_RATIO = 0.12

IMG_SIZE = 224
BATCH_SIZE = 32
NODDING_HIGH = 0.70
UPWARD_HIGH = 0.30

OUTDIR = BASE / "global_cirsium_flowering_annotation_then_yolo"
TMP_DIR = OUTDIR / "tmp_screen_images"
FLOWER_IMG_DIR = OUTDIR / "flower_only_images"
CROP_DIR = OUTDIR / "head_crops"
MAP_DIR = OUTDIR / "maps"

for d in [OUTDIR, TMP_DIR, FLOWER_IMG_DIR, CROP_DIR, MAP_DIR]:
    d.mkdir(parents=True, exist_ok=True)

STATE_JSON = OUTDIR / "download_state.json"
ALL_PHOTO_META_CSV = OUTDIR / "flowering_annotation_photo_metadata_seen.csv"
FLOWER_PHOTO_META_CSV = OUTDIR / "flowering_annotation_yolo_flower_photo_metadata.csv"
SCREENING_CSV = OUTDIR / "flowering_annotation_yolo_screening_results.csv"
CROP_META_CSV = OUTDIR / "flowering_annotation_yolo_crop_metadata.csv"
CROP_PRED_CSV = OUTDIR / "flowering_annotation_crop_predictions.csv"
OBS_SUMMARY_CSV = OUTDIR / "flowering_annotation_observation_head_direction_summary.csv"

MAP_ALL = MAP_DIR / "global_cirsium_flowering_annotation_head_direction_high_confidence_map.png"
MAP_NODDING = MAP_DIR / "global_cirsium_flowering_annotation_nodding_high_confidence_map.png"
MAP_UPWARD = MAP_DIR / "global_cirsium_flowering_annotation_upward_high_confidence_map.png"


# ============================================================
# HTTP helpers
# ============================================================
def make_session():
    session = requests.Session()
    retry = Retry(
        total=2,
        connect=2,
        read=2,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=30, pool_maxsize=30)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({
        "User-Agent": "global-cirsium-flower-only-pipeline/0.3 "
                      "(research use; contact: rachelzhang0223@gmail.com)"
    })
    return session


SESSION = make_session()
DNS_FAIL_STREAK = 0
DNS_FAIL_STOP = 80


def is_dns_error(err):
    msg = str(err)
    return (
        "NameResolutionError" in msg
        or "getaddrinfo failed" in msg
        or "Failed to resolve" in msg
        or "Temporary failure in name resolution" in msg
    )


def safe_get_json(url, params=None, max_retry=5):
    last_err = None
    for attempt in range(1, max_retry + 1):
        try:
            r = SESSION.get(url, params=params, timeout=(20, 60))
            if r.status_code == 429:
                wait = 15 * attempt
                print(f"[WARN] API 429 rate limited. sleeping {wait}s")
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_err = e
            wait = 3 * attempt
            print(f"[WARN] API request failed attempt {attempt}/{max_retry}: {e}")
            time.sleep(wait)
    raise RuntimeError(f"API request failed: {url} params={params} err={last_err}")


def prefer_image_size_url(url, preferred="small"):
    if url is None:
        return None
    if preferred in ["square", "small", "medium", "large", "original"]:
        return re.sub(r"/(square|small|medium|large|original)\.", f"/{preferred}.", str(url))
    return str(url)


def safe_download(url, out_path, max_retry=2):
    global DNS_FAIL_STREAK
    out_path = Path(out_path)

    if out_path.exists() and out_path.stat().st_size > 1000:
        return True
    if out_path.exists() and out_path.stat().st_size <= 1000:
        out_path.unlink(missing_ok=True)

    last_err = None
    for attempt in range(1, max_retry + 1):
        try:
            r = SESSION.get(url, timeout=(15, 35), stream=True)

            if r.status_code == 404:
                print(f"[SKIP] 404 image not found: {url}")
                DNS_FAIL_STREAK = 0
                return False
            if r.status_code == 429:
                wait = 15 * attempt
                print(f"[WARN] image 429 rate limited. sleeping {wait}s")
                time.sleep(wait)
                continue

            r.raise_for_status()
            content_type = r.headers.get("Content-Type", "")
            if "image" not in content_type.lower():
                print(f"[SKIP] not image: {url} content-type={content_type}")
                DNS_FAIL_STREAK = 0
                return False

            out_path.parent.mkdir(parents=True, exist_ok=True)
            n_bytes = 0
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 128):
                    if chunk:
                        f.write(chunk)
                        n_bytes += len(chunk)

            if n_bytes < 1000 or out_path.stat().st_size < 1000:
                print(f"[SKIP] tiny/broken image: {url}")
                out_path.unlink(missing_ok=True)
                DNS_FAIL_STREAK = 0
                return False

            DNS_FAIL_STREAK = 0
            time.sleep(DOWNLOAD_SLEEP)
            return True

        except requests.exceptions.ConnectionError as e:
            last_err = e
            if is_dns_error(e):
                DNS_FAIL_STREAK += 1
                print(f"[DNS SKIP] {DNS_FAIL_STREAK}/{DNS_FAIL_STOP}: {url}")
                time.sleep(1.5)
                if DNS_FAIL_STREAK >= DNS_FAIL_STOP:
                    raise RuntimeError(
                        "Too many consecutive DNS failures. Try ipconfig /flushdns, change Wi-Fi/VPN, "
                        "or use smartphone tethering, then rerun."
                    )
                return False
            print(f"[WARN] connection error attempt {attempt}/{max_retry}: {e}")
            time.sleep(2 * attempt)

        except requests.exceptions.ReadTimeout as e:
            last_err = e
            print(f"[TIMEOUT SKIP] {url}")
            time.sleep(1)
            return False

        except Exception as e:
            last_err = e
            print(f"[WARN] image download failed attempt {attempt}/{max_retry}: {e}")
            time.sleep(2 * attempt)

    print(f"[WARN] failed image download: {url} err={last_err}")
    return False


# ============================================================
# Generic helpers
# ============================================================
def append_csv_rows(rows, path):
    if not rows:
        return
    df = pd.DataFrame(rows)
    header = not Path(path).exists()
    df.to_csv(path, mode="a", header=header, index=False, encoding="utf-8-sig")


def load_csv_or_empty(path):
    if Path(path).exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def find_existing(candidates, name):
    for p in candidates:
        if p.exists():
            print(f"[INFO] found {name}: {p}")
            return p
    raise FileNotFoundError(f"{name} not found. Checked:\n" + "\n".join(str(p) for p in candidates))


def ext_from_url(url):
    path = urlparse(str(url)).path
    ext = Path(path).suffix.lower()
    if ext in [".jpg", ".jpeg", ".png", ".webp"]:
        return ext
    return ".jpg"


def safe_crop_box(x1, y1, x2, y2, w, h, pad_ratio=0.12):
    bw = x2 - x1
    bh = y2 - y1
    pad_x = bw * pad_ratio
    pad_y = bh * pad_ratio
    x1 = max(0, int(x1 - pad_x))
    y1 = max(0, int(y1 - pad_y))
    x2 = min(w, int(x2 + pad_x))
    y2 = min(h, int(y2 + pad_y))
    return x1, y1, x2, y2


# ============================================================
# iNaturalist metadata
# ============================================================
def resolve_inat_taxon_id(name="Cirsium"):
    data = safe_get_json(
        "https://api.inaturalist.org/v1/taxa",
        {"q": name, "rank": "genus", "is_active": "true", "per_page": 10},
    )
    results = data.get("results", [])
    for t in results:
        if str(t.get("name", "")).lower() == name.lower() and t.get("rank") == "genus":
            print(f"[INFO] resolved {name} taxon_id: {t.get('id')}")
            return int(t["id"])
    if not results:
        raise RuntimeError(f"No taxon found for {name}")
    print("[WARN] exact genus match not found; using first result")
    print(json.dumps(results[0], ensure_ascii=False, indent=2)[:1000])
    return int(results[0]["id"])


def make_photo_rows(obs):
    rows = []
    obs_id = obs.get("id")
    taxon = obs.get("taxon") or {}
    user = obs.get("user") or {}
    geojson = obs.get("geojson") or {}
    coords = geojson.get("coordinates", [None, None])
    lon = coords[0] if len(coords) > 0 else None
    lat = coords[1] if len(coords) > 1 else None

    if EXCLUDE_OBSCURED:
        if obs.get("geoprivacy") in ["obscured", "private"] or obs.get("obscured"):
            return []

    annotations = obs.get("annotations") or []
    has_flowering_annotation = any(
        (a.get("controlled_attribute", {}).get("id") == FLOWERS_FRUITS_TERM_ID
         and a.get("controlled_value", {}).get("id") == FLOWERS_TERM_VALUE_ID)
        for a in annotations
    )
    # Safe annotation summary. Avoid nested quotes inside f-strings.
    annotation_parts = []
    for a in annotations:
        attr = a.get("controlled_attribute") or {}
        val = a.get("controlled_value") or {}
        attr_label = attr.get("label") or attr.get("id")
        val_label = val.get("label") or val.get("id")
        annotation_parts.append(f"{attr_label}={val_label}")
    annotation_summary = "; ".join(map(str, annotation_parts))

    for k, photo in enumerate(obs.get("photos") or []):
        photo_id = photo.get("id")
        raw_url = photo.get("url") or photo.get("medium_url") or photo.get("large_url")
        if not raw_url:
            continue
        screen_url = prefer_image_size_url(raw_url, SCREEN_IMAGE_SIZE)
        keep_url = prefer_image_size_url(raw_url, KEEP_IMAGE_SIZE)
        ext = ext_from_url(screen_url)
        source_name = f"obs_{obs_id}_photo_{photo_id}{ext}"
        temp_path = TMP_DIR / source_name
        flower_path = FLOWER_IMG_DIR / source_name

        rows.append({
            "obs_id": obs_id,
            "photo_id": photo_id,
            "photo_index": k,
            "taxon_id": taxon.get("id"),
            "taxon_name": taxon.get("name"),
            "preferred_common_name": taxon.get("preferred_common_name"),
            "rank": taxon.get("rank"),
            "observed_on": obs.get("observed_on"),
            "created_at": obs.get("created_at"),
            "latitude": lat,
            "longitude": lon,
            "positional_accuracy": obs.get("positional_accuracy"),
            "quality_grade": obs.get("quality_grade"),
            "geoprivacy": obs.get("geoprivacy"),
            "obscured": obs.get("obscured"),
            "place_guess": obs.get("place_guess"),
            "user_login": user.get("login"),
            "license_code": obs.get("license_code"),
            "photo_license_code": photo.get("license_code"),
            "inat_flowers_annotation": has_flowering_annotation,
            "inat_annotation_summary": annotation_summary,
            "raw_image_url": raw_url,
            "screen_image_url": screen_url,
            "keep_image_url": keep_url,
            "source_name": source_name,
            "temp_path": str(temp_path),
            "flower_local_path": str(flower_path),
        })
    return rows


def load_state():
    if STATE_JSON.exists():
        with open(STATE_JSON, "r", encoding="utf-8") as f:
            return json.load(f)
    return {"last_obs_id": 0, "n_obs_seen": 0, "n_photos_seen": 0}


def save_state(state):
    tmp = STATE_JSON.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)
    tmp.replace(STATE_JSON)


# ============================================================
# YOLO screen and crop on the fly
# ============================================================
def screen_one_photo_with_yolo(row, yolo_model):
    """Download one photo, run YOLO, keep only flower/head-positive images and crops."""
    source_name = row["source_name"]
    temp_path = Path(row["temp_path"])
    flower_path = Path(row["flower_local_path"])
    screen_url = row["screen_image_url"]
    keep_url = row["keep_image_url"]

    # If already has flower image and crop metadata, skip not handled here; caller checks processed set.
    ok = safe_download(screen_url, temp_path)
    if not ok:
        return {
            "source_name": source_name,
            "obs_id": row["obs_id"],
            "photo_id": row["photo_id"],
            "download_ok": False,
            "has_flower_yolo": False,
            "n_detections": 0,
            "screen_status": "download_failed",
        }, []

    crop_rows = []
    try:
        im = Image.open(temp_path).convert("RGB")
        W, H = im.size
        results = yolo_model.predict(source=str(temp_path), conf=YOLO_CONF_THRES, verbose=False)

        if len(results) == 0 or results[0].boxes is None or len(results[0].boxes) == 0:
            if DELETE_NON_FLOWER_TEMP:
                temp_path.unlink(missing_ok=True)
            return {
                "source_name": source_name,
                "obs_id": row["obs_id"],
                "photo_id": row["photo_id"],
                "download_ok": True,
                "has_flower_yolo": False,
                "n_detections": 0,
                "screen_status": "no_flower_detected",
            }, []

        # Flower/head detected. Keep image.
        # If KEEP_IMAGE_SIZE == SCREEN_IMAGE_SIZE, just move/copy temp to flower path.
        flower_path.parent.mkdir(parents=True, exist_ok=True)
        if KEEP_IMAGE_SIZE == SCREEN_IMAGE_SIZE:
            if not flower_path.exists():
                temp_path.replace(flower_path)
            else:
                temp_path.unlink(missing_ok=True)
            kept_image_path = flower_path
        else:
            # Redownload better image, but fall back to screen image if medium/original fails.
            if safe_download(keep_url, flower_path):
                temp_path.unlink(missing_ok=True)
                kept_image_path = flower_path
                im = Image.open(kept_image_path).convert("RGB")
                W, H = im.size
                results = yolo_model.predict(source=str(kept_image_path), conf=YOLO_CONF_THRES, verbose=False)
            else:
                if not flower_path.exists():
                    temp_path.replace(flower_path)
                kept_image_path = flower_path

        # Crop from the kept image and current YOLO result.
        im = Image.open(kept_image_path).convert("RGB")
        W, H = im.size
        results = yolo_model.predict(source=str(kept_image_path), conf=YOLO_CONF_THRES, verbose=False)
        boxes = results[0].boxes if len(results) > 0 else []

        for j, box in enumerate(boxes):
            xyxy = box.xyxy[0].detach().cpu().numpy()
            conf = float(box.conf[0].detach().cpu().numpy()) if box.conf is not None else np.nan
            cls = int(box.cls[0].detach().cpu().numpy()) if box.cls is not None else -1
            x1, y1, x2, y2 = map(float, xyxy)
            cx1, cy1, cx2, cy2 = safe_crop_box(x1, y1, x2, y2, W, H, PAD_RATIO)
            crop = im.crop((cx1, cy1, cx2, cy2))
            crop_name = f"{Path(source_name).stem}_det{j:02d}_conf{conf:.3f}.jpg"
            crop_path = CROP_DIR / crop_name
            crop.save(crop_path, quality=95)

            cr = dict(row)
            cr.update({
                "source_image": str(kept_image_path),
                "crop_path": str(crop_path),
                "crop_name": crop_name,
                "det_index": j,
                "yolo_conf": conf,
                "yolo_class": cls,
                "image_width": W,
                "image_height": H,
                "x1": x1,
                "y1": y1,
                "x2": x2,
                "y2": y2,
                "crop_x1": cx1,
                "crop_y1": cy1,
                "crop_x2": cx2,
                "crop_y2": cy2,
            })
            crop_rows.append(cr)

        return {
            "source_name": source_name,
            "obs_id": row["obs_id"],
            "photo_id": row["photo_id"],
            "download_ok": True,
            "has_flower_yolo": True,
            "n_detections": len(crop_rows),
            "screen_status": "flower_detected",
            "flower_local_path": str(kept_image_path),
        }, crop_rows

    except Exception as e:
        print(f"[WARN] YOLO screening failed: {source_name} | {e}")
        temp_path.unlink(missing_ok=True)
        return {
            "source_name": source_name,
            "obs_id": row["obs_id"],
            "photo_id": row["photo_id"],
            "download_ok": True,
            "has_flower_yolo": False,
            "n_detections": 0,
            "screen_status": "yolo_error",
            "error": str(e),
        }, []


def collect_flower_only_global_photos():
    taxon_id = resolve_inat_taxon_id(GENUS_NAME)
    yolo_path = find_existing(YOLO_CANDIDATES, "YOLO model")
    yolo_model = YOLO(str(yolo_path))

    state = load_state()
    last_obs_id = int(state.get("last_obs_id", 0) or 0)
    n_obs_seen = int(state.get("n_obs_seen", 0) or 0)
    n_photos_seen = int(state.get("n_photos_seen", 0) or 0)

    screening_old = load_csv_or_empty(SCREENING_CSV)
    processed_sources = set(screening_old["source_name"].astype(str)) if "source_name" in screening_old.columns else set()
    print(f"[INFO] resume last_obs_id={last_obs_id} processed_photos={len(processed_sources)}")

    batch_index = 0
    total_flower_photos = 0
    total_crops = 0

    while True:
        params = {
            "taxon_id": taxon_id,
            "per_page": PER_PAGE,
            "order_by": "id",
            "order": "asc",
            "photos": "true" if ONLY_PHOTOS else "false",
        }
        if last_obs_id > 0:
            params["id_above"] = last_obs_id
        if ONLY_GEOTAGGED:
            params["geo"] = "true"
        if USE_INAT_FLOWERING_ANNOTATION:
            # iNat controlled terms: Flowers and Fruits = 12, Flowers = 13
            params["term_id"] = FLOWERS_FRUITS_TERM_ID
            params["term_value_id"] = FLOWERS_TERM_VALUE_ID

        data = safe_get_json("https://api.inaturalist.org/v1/observations", params)
        results = data.get("results", [])
        if not results:
            print("[INFO] no more observations returned")
            break

        batch_index += 1
        page_photo_rows = []
        for obs in results:
            page_photo_rows.extend(make_photo_rows(obs))

        n_obs_seen += len(results)
        n_photos_seen += len(page_photo_rows)
        last_obs_id = max(int(o.get("id")) for o in results if o.get("id") is not None)

        print(
            f"[INFO] batch={batch_index} obs={len(results)} photo_rows={len(page_photo_rows)} "
            f"last_obs_id={last_obs_id} total_obs_seen={n_obs_seen} total_photos_seen={n_photos_seen}"
        )

        # Save all seen photo metadata, including non-flower photos.
        append_csv_rows(page_photo_rows, ALL_PHOTO_META_CSV)

        screen_rows = []
        flower_rows = []
        crop_rows = []

        for idx, row in enumerate(page_photo_rows, start=1):
            source_name = str(row["source_name"])
            if source_name in processed_sources:
                continue

            screen_result, crops = screen_one_photo_with_yolo(row, yolo_model)
            screen_rows.append(screen_result)
            processed_sources.add(source_name)

            if screen_result.get("has_flower_yolo"):
                fr = dict(row)
                fr.update(screen_result)
                flower_rows.append(fr)
                crop_rows.extend(crops)

            if idx % 50 == 0:
                print(
                    f"[INFO] screened {idx}/{len(page_photo_rows)} in batch | "
                    f"flower_photos_batch={len(flower_rows)} crops_batch={len(crop_rows)}"
                )

            if MAX_PHOTOS is not None and len(processed_sources) >= MAX_PHOTOS:
                print(f"[INFO] reached MAX_PHOTOS={MAX_PHOTOS}")
                break

        append_csv_rows(screen_rows, SCREENING_CSV)
        append_csv_rows(flower_rows, FLOWER_PHOTO_META_CSV)
        append_csv_rows(crop_rows, CROP_META_CSV)

        total_flower_photos += len(flower_rows)
        total_crops += len(crop_rows)

        state = {
            "last_obs_id": last_obs_id,
            "n_obs_seen": n_obs_seen,
            "n_photos_seen": n_photos_seen,
            "processed_photos": len(processed_sources),
            "updated_at_unix": time.time(),
        }
        save_state(state)

        print(
            f"[INFO] checkpoint saved | flower_photos_this_batch={len(flower_rows)} "
            f"crops_this_batch={len(crop_rows)} processed_photos={len(processed_sources)}"
        )

        if MAX_OBS is not None and n_obs_seen >= MAX_OBS:
            print(f"[INFO] reached MAX_OBS={MAX_OBS}")
            break
        if MAX_PHOTOS is not None and len(processed_sources) >= MAX_PHOTOS:
            print(f"[INFO] reached MAX_PHOTOS={MAX_PHOTOS}")
            break

        time.sleep(REQUEST_SLEEP)

    print("\n=== flower-only collection done ===")
    print(f"all photo metadata : {ALL_PHOTO_META_CSV}")
    print(f"screening results  : {SCREENING_CSV}")
    print(f"flower photo meta  : {FLOWER_PHOTO_META_CSV}")
    print(f"crop metadata      : {CROP_META_CSV}")

    return pd.read_csv(FLOWER_PHOTO_META_CSV) if FLOWER_PHOTO_META_CSV.exists() else pd.DataFrame()


# ============================================================
# Direction prediction
# ============================================================
def load_images_numpy(paths):
    imgs = []
    for p in paths:
        im = Image.open(p).convert("RGB").resize((IMG_SIZE, IMG_SIZE))
        arr = np.asarray(im).astype("float32")
        imgs.append(arr)
    return np.stack(imgs, axis=0)


def predict_in_batches(model, paths):
    probs = []
    for i in range(0, len(paths), BATCH_SIZE):
        batch_paths = paths[i:i + BATCH_SIZE]
        X = load_images_numpy(batch_paths)
        p = model.predict(X, verbose=0).reshape(-1)
        probs.append(p)
        if i == 0 or (i // BATCH_SIZE) % 10 == 0:
            print(f"[INFO] predicted crops {min(i + BATCH_SIZE, len(paths))}/{len(paths)}")
    return np.concatenate(probs) if probs else np.array([])


def decision_from_prob(p):
    if p >= NODDING_HIGH:
        return "nodding_high"
    if p <= UPWARD_HIGH:
        return "upward_high"
    return "uncertain"


def confidence_from_prob(p):
    return max(float(p), 1.0 - float(p))


def predict_crop_direction():
    if not CROP_META_CSV.exists():
        print(f"[WARN] no crop metadata found: {CROP_META_CSV}")
        return pd.DataFrame()

    crop_df = pd.read_csv(CROP_META_CSV)
    if len(crop_df) == 0:
        print("[WARN] crop metadata is empty")
        return crop_df

    if not HEAD_MODEL.exists():
        raise FileNotFoundError(f"trained crop classifier not found: {HEAD_MODEL}")

    old_pred = load_csv_or_empty(CROP_PRED_CSV)
    done_crops = set(old_pred["crop_name"].astype(str)) if "crop_name" in old_pred.columns else set()

    todo = crop_df[~crop_df["crop_name"].astype(str).isin(done_crops)].copy()
    todo["crop_exists"] = todo["crop_path"].astype(str).apply(lambda x: Path(x).exists())
    todo = todo[todo["crop_exists"]].copy()

    print(f"[INFO] total crops={len(crop_df)} already_predicted={len(done_crops)} todo={len(todo)}")
    if len(todo) == 0:
        return old_pred

    model = keras.models.load_model(HEAD_MODEL)
    paths = [Path(x) for x in todo["crop_path"].astype(str).tolist()]
    prob = predict_in_batches(model, paths)

    todo["pred_prob_nodding"] = prob
    todo["pred_label_050"] = np.where(todo["pred_prob_nodding"] >= 0.5, "nodding", "upward")
    todo["confidence"] = todo["pred_prob_nodding"].apply(confidence_from_prob)
    todo["decision"] = todo["pred_prob_nodding"].apply(decision_from_prob)

    append_csv_rows(todo.to_dict("records"), CROP_PRED_CSV)
    out = pd.read_csv(CROP_PRED_CSV)

    print("\n=== crop prediction summary ===")
    print(out["decision"].value_counts(dropna=False))
    print(f"[INFO] saved crop predictions: {CROP_PRED_CSV}")
    return out


# ============================================================
# Aggregate and maps
# ============================================================
def aggregate_to_observations():
    if not FLOWER_PHOTO_META_CSV.exists() or not CROP_PRED_CSV.exists():
        print("[WARN] missing flower metadata or crop predictions")
        return pd.DataFrame()

    flower_df = pd.read_csv(FLOWER_PHOTO_META_CSV)
    pred_df = pd.read_csv(CROP_PRED_CSV)

    rows = []
    for source_name, g in pred_df.groupby("source_name", dropna=False):
        n_nod = int((g["decision"] == "nodding_high").sum())
        n_up = int((g["decision"] == "upward_high").sum())
        n_unc = int((g["decision"] == "uncertain").sum())
        best_idx = g["confidence"].idxmax()
        best = g.loc[best_idx]

        if n_nod > 0 and n_up > 0:
            source_decision = "mixed_high"
        elif n_nod > 0:
            source_decision = "nodding_high"
        elif n_up > 0 and n_unc == 0:
            source_decision = "upward_high"
        elif n_up > 0 and n_unc > 0:
            source_decision = "upward_high_with_uncertain"
        else:
            source_decision = "uncertain"

        rows.append({
            "source_name": source_name,
            "n_crops": len(g),
            "source_decision": source_decision,
            "best_crop_path": best["crop_path"],
            "best_prob_nodding": best["pred_prob_nodding"],
            "best_confidence": best["confidence"],
            "best_label_050": best["pred_label_050"],
            "max_prob_nodding": g["pred_prob_nodding"].max(),
            "min_prob_nodding": g["pred_prob_nodding"].min(),
            "mean_prob_nodding": g["pred_prob_nodding"].mean(),
            "n_nodding_high": n_nod,
            "n_upward_high": n_up,
            "n_uncertain": n_unc,
        })

    summary = pd.DataFrame(rows)
    out = flower_df.drop_duplicates("source_name").merge(summary, on="source_name", how="left")
    out = out.dropna(subset=["source_decision"]).copy()

    out["head_direction_conservative"] = np.nan
    out.loc[out["source_decision"] == "nodding_high", "head_direction_conservative"] = "nodding"
    out.loc[out["source_decision"] == "upward_high", "head_direction_conservative"] = "upward"

    out["head_direction_lenient"] = np.nan
    out.loc[out["source_decision"].isin(["nodding_high"]), "head_direction_lenient"] = "nodding"
    out.loc[out["source_decision"].isin(["upward_high", "upward_high_with_uncertain"]), "head_direction_lenient"] = "upward"

    out.to_csv(OBS_SUMMARY_CSV, index=False, encoding="utf-8-sig")
    print("\n=== observation summary ===")
    print(out["source_decision"].value_counts(dropna=False))
    print("\n=== conservative labels ===")
    print(out["head_direction_conservative"].value_counts(dropna=False))
    print(f"[INFO] saved observation summary: {OBS_SUMMARY_CSV}")
    return out


def plot_world_points(df, out_path, title, label_col="head_direction_conservative"):
    d = df.dropna(subset=["latitude", "longitude", label_col]).copy()
    if len(d) == 0:
        print(f"[WARN] no data for map: {title}")
        return
    d["latitude"] = pd.to_numeric(d["latitude"], errors="coerce")
    d["longitude"] = pd.to_numeric(d["longitude"], errors="coerce")
    d = d.dropna(subset=["latitude", "longitude"])
    d = d[(d["latitude"].between(-90, 90)) & (d["longitude"].between(-180, 180))]

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    ax.grid(True, linewidth=0.4, alpha=0.4)

    up = d[d[label_col] == "upward"]
    nod = d[d[label_col] == "nodding"]

    if len(up) > 0:
        ax.scatter(up["longitude"], up["latitude"], s=8, alpha=0.55, label=f"upward (n={len(up)})")
    if len(nod) > 0:
        ax.scatter(nod["longitude"], nod["latitude"], s=10, alpha=0.70, marker="^", label=f"nodding (n={len(nod)})")

    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"[INFO] saved map: {out_path}")


def make_maps(obs_summary):
    if obs_summary is None or len(obs_summary) == 0:
        print("[WARN] no observation summary for maps")
        return
    plot_world_points(obs_summary, MAP_ALL, "Global Cirsium flower-only head orientation predictions")
    plot_world_points(
        obs_summary[obs_summary["head_direction_conservative"] == "nodding"],
        MAP_NODDING,
        "Global Cirsium flower-only predicted nodding records",
    )
    plot_world_points(
        obs_summary[obs_summary["head_direction_conservative"] == "upward"],
        MAP_UPWARD,
        "Global Cirsium flower-only predicted upward records",
    )


# ============================================================
# Main
# ============================================================
def main():
    print("\n==============================================")
    print("GLOBAL CIRSIUM FLOWERING-ANNOTATION THEN YOLO PIPELINE")
    print("==============================================\n")
    print(f"[INFO] OUTDIR: {OUTDIR}")
    print(f"[INFO] MAX_OBS={MAX_OBS} MAX_PHOTOS={MAX_PHOTOS}")
    print(f"[INFO] iNat flowering annotation filter={USE_INAT_FLOWERING_ANNOTATION} term_id={FLOWERS_FRUITS_TERM_ID} term_value_id={FLOWERS_TERM_VALUE_ID}")

    collect_flower_only_global_photos()
    pred_df = predict_crop_direction()
    obs_summary = aggregate_to_observations()
    make_maps(obs_summary)

    print("\n==============================")
    print("DONE")
    print("==============================")
    print(f"All photo metadata seen: {ALL_PHOTO_META_CSV}")
    print(f"YOLO screening results : {SCREENING_CSV}")
    print(f"Flower-only metadata   : {FLOWER_PHOTO_META_CSV}")
    print(f"Flower-only images     : {FLOWER_IMG_DIR}")
    print(f"Head crops             : {CROP_DIR}")
    print(f"Crop predictions       : {CROP_PRED_CSV}")
    print(f"Observation summary    : {OBS_SUMMARY_CSV}")
    print(f"Maps                   : {MAP_DIR}")


if __name__ == "__main__":
    main()
