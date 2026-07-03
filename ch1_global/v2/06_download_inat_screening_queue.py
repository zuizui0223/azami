#!/usr/bin/env python3
"""Download a pre-built Ch.1 image-screening queue, resumably and without taxonomy in filenames.

Run only after inspecting `queue_species_qc.csv` and `queue_spatial_qc.csv` from
05_build_image_screening_queue.py. This downloader does not select images, run
YOLO, or create trait labels.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


DOWNLOAD_VERSION = "1.0.0"
RESULT_FIELDS = ["queue_id", "photo_id", "download_status", "http_status", "bytes_written", "image_local_path", "image_url", "message", "finished_at_utc"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download a Ch.1 iNaturalist screening queue.")
    parser.add_argument("--queue", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--image-size", choices=["small", "medium", "large"], default="small")
    parser.add_argument("--max-images", type=int, default=0, help="0 means all queue rows that are not already successful.")
    parser.add_argument("--sleep-sec", type=float, default=0.12)
    parser.add_argument("--timeout-sec", type=int, default=60)
    parser.add_argument("--user-agent", default="azami-ch1-cirsium-image-download/1.0 (research use; contact: rachelzhang0223@gmail.com)")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def make_session(user_agent: str) -> requests.Session:
    session = requests.Session()
    retry = Retry(total=4, connect=4, read=4, backoff_factor=1.2, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"], raise_on_status=False)
    adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": user_agent})
    return session


def ext_from_url(url: str) -> str:
    match = re.search(r"\.(jpg|jpeg|png|webp)(?:$|\?)", url.lower())
    return f".{match.group(1)}" if match else ".jpg"


def selected_url(row: pd.Series, size: str) -> str:
    preferred = text(row.get(f"{size}_image_url", ""))
    if preferred:
        return preferred
    for candidate in ["small_image_url", "medium_image_url", "large_image_url"]:
        value = text(row.get(candidate, ""))
        if value:
            return value
    return ""


def load_successes(results_path: Path) -> set[str]:
    if not results_path.exists():
        return set()
    try:
        table = pd.read_csv(results_path, dtype=str)
    except pd.errors.EmptyDataError:
        return set()
    if not {"queue_id", "download_status"}.issubset(table.columns):
        raise ValueError(f"Existing results missing queue_id/download_status: {results_path}")
    return set(table.loc[table["download_status"].eq("success"), "queue_id"].dropna().astype(str))


def append_result(results_path: Path, result: dict[str, Any]) -> None:
    write_header = not results_path.exists() or results_path.stat().st_size == 0
    with results_path.open("a", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(result)


def download_one(session: requests.Session, url: str, destination: Path, timeout_sec: int) -> tuple[str, int | None, int, str]:
    if destination.exists() and destination.stat().st_size > 1000:
        return "success", None, int(destination.stat().st_size), "already present"
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with session.get(url, stream=True, timeout=(15, timeout_sec)) as response:
            status = response.status_code
            if status == 429:
                return "rate_limited", status, 0, "HTTP 429"
            if status >= 400:
                return "failed", status, 0, f"HTTP {status}"
            content_type = text(response.headers.get("Content-Type", "")).lower()
            if "image" not in content_type:
                return "failed", status, 0, f"non-image content type: {content_type}"
            destination.parent.mkdir(parents=True, exist_ok=True)
            n_bytes = 0
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 128):
                    if chunk:
                        handle.write(chunk)
                        n_bytes += len(chunk)
        if n_bytes < 1000:
            temporary.unlink(missing_ok=True)
            return "failed", status, n_bytes, "downloaded file is too small"
        os.replace(temporary, destination)
        return "success", status, n_bytes, ""
    except requests.RequestException as error:
        temporary.unlink(missing_ok=True)
        return "failed", None, 0, str(error)


def main() -> None:
    args = parse_args()
    if args.max_images < 0 or args.sleep_sec < 0:
        raise ValueError("--max-images and --sleep-sec must be >= 0")
    queue_path = Path(args.queue)
    queue = pd.read_csv(queue_path, dtype=str, keep_default_na=False)
    required = {"queue_id", "photo_id", "screen_download_filename"}
    missing = required.difference(queue.columns)
    if missing:
        raise ValueError(f"Queue lacks required columns: {sorted(missing)}")
    if queue["queue_id"].duplicated().any():
        raise ValueError("Queue has duplicate queue_id values")

    out_dir = Path(args.out_dir)
    images_dir = out_dir / "images"
    results_path = out_dir / "download_results.csv"
    manifest_path = out_dir / "download_manifest.json"
    out_dir.mkdir(parents=True, exist_ok=True)
    images_dir.mkdir(parents=True, exist_ok=True)
    session = make_session(args.user_agent)
    successes = load_successes(results_path)
    candidates = queue.loc[~queue["queue_id"].astype(str).isin(successes)].copy()
    if args.max_images:
        candidates = candidates.head(args.max_images)
    print(f"[INFO] queued={len(queue)} existing_success={len(successes)} to_attempt={len(candidates)}")

    attempted = 0
    n_success = 0
    for _, row in candidates.iterrows():
        url = selected_url(row, args.image_size)
        filename = text(row["screen_download_filename"])
        if not filename:
            filename = f"inat_photo_{text(row['photo_id'])}{ext_from_url(url)}"
        destination = images_dir / filename
        if not url:
            result = {"queue_id": row["queue_id"], "photo_id": row["photo_id"], "download_status": "failed", "http_status": "", "bytes_written": 0, "image_local_path": "", "image_url": "", "message": "no image URL", "finished_at_utc": utc_now()}
        else:
            status, http_status, n_bytes, message = download_one(session, url, destination, args.timeout_sec)
            result = {"queue_id": row["queue_id"], "photo_id": row["photo_id"], "download_status": status, "http_status": http_status or "", "bytes_written": n_bytes, "image_local_path": str(destination.resolve()) if status == "success" else "", "image_url": url, "message": message, "finished_at_utc": utc_now()}
        append_result(results_path, result)
        attempted += 1
        n_success += int(result["download_status"] == "success")
        print(f"[INFO] {attempted}/{len(candidates)} {row['queue_id']} -> {result['download_status']}")
        if result["download_status"] == "rate_limited":
            raise RuntimeError("Stopped after HTTP 429. Wait before resuming; successful rows will be skipped.")
        time.sleep(args.sleep_sec)

    manifest_path.write_text(json.dumps({"download_version": DOWNLOAD_VERSION, "finished_at_utc": utc_now(), "queue": str(queue_path.resolve()), "out_dir": str(out_dir.resolve()), "image_size": args.image_size, "attempted_this_run": attempted, "successful_this_run": n_success}, ensure_ascii=False, indent=2), encoding="utf-8")
    print("[OK] queue download finished", out_dir.resolve())


if __name__ == "__main__":
    main()
