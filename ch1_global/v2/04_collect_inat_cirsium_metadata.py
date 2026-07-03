#!/usr/bin/env python3
"""Collect reproducible global iNaturalist metadata for Chapter 1 v2.

This stage deliberately collects **metadata only**. It does not download images,
run YOLO, or infer traits. That separation prevents an old classifier from
silently deciding which observations exist in the Chapter 1 sampling frame.

Outputs
-------
- observation_raw.ndjson: one complete API observation object per line
- photo_metadata.csv: one row per photo, with observation/taxon/licence/geo fields
- checkpoint.json: restart state using stable observation ID pagination
- collection_provenance.json: exact endpoint, taxon resolution, parameters, and time
- collection_summary.csv: a compact accounting table

Examples
--------
# Pilot: stop after 5,000 observations.
python ch1_global/v2/04_collect_inat_cirsium_metadata.py \
  --out-dir "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_raw" \
  --max-observations 5000

# Full metadata inventory: 0 means no imposed observation limit.
python ch1_global/v2/04_collect_inat_cirsium_metadata.py \
  --out-dir "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_raw" \
  --max-observations 0
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


COLLECTOR_VERSION = "1.0.0"
TAXA_ENDPOINT = "https://api.inaturalist.org/v1/taxa"
OBSERVATIONS_ENDPOINT = "https://api.inaturalist.org/v1/observations"
PHOTO_FIELDS = [
    "obs_id", "photo_id", "photo_index", "photo_key",
    "taxon_id", "taxon_name", "taxon_rank", "species_guess",
    "preferred_common_name", "observed_on", "created_at", "updated_at",
    "quality_grade", "captive", "user_id", "user_login",
    "latitude", "longitude", "positional_accuracy", "geoprivacy", "obscured",
    "coordinate_usable_for_environment", "inat_flowers_annotation", "inat_annotation_summary",
    "observation_license_code", "photo_license_code", "photo_attribution",
    "raw_image_url", "small_image_url", "medium_image_url", "large_image_url",
    "source_file_stem",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect global Cirsium iNaturalist photo metadata only.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--taxon-name", default="Cirsium")
    parser.add_argument("--taxon-id", type=int, default=None, help="Skip taxon lookup only when an exact current iNaturalist genus ID is known.")
    parser.add_argument("--per-page", type=int, default=200)
    parser.add_argument("--max-observations", type=int, default=5000, help="0 means no imposed limit; start with a pilot first.")
    parser.add_argument("--quality-grade", choices=["any", "research", "needs_id", "casual"], default="any")
    parser.add_argument("--photo-license", default="any", help="Optional iNaturalist photo_license filter, e.g. cc-by or cc0.")
    parser.add_argument("--sleep-sec", type=float, default=0.8)
    parser.add_argument("--timeout-sec", type=int, default=90)
    parser.add_argument("--restart", action="store_true", help="Start a new collection only in a new output directory; existing files are never overwritten.")
    parser.add_argument("--user-agent", default="azami-ch1-cirsium-metadata/1.0 (research use; contact: rachelzhang0223@gmail.com)")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def as_bool(value: Any) -> bool:
    return bool(value is True or text(value).lower() in {"true", "1", "yes"})


def finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def make_session(user_agent: str) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": user_agent, "Accept": "application/json"})
    return session


def get_json(session: requests.Session, url: str, params: dict[str, Any], timeout_sec: int) -> dict[str, Any]:
    response = session.get(url, params=params, timeout=(20, timeout_sec))
    if response.status_code == 429:
        retry_after = response.headers.get("Retry-After", "")
        raise RuntimeError(f"iNaturalist returned HTTP 429. Pause and rerun later. Retry-After={retry_after}")
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, dict):
        raise RuntimeError(f"Expected a JSON object from {url}, received {type(payload).__name__}")
    return payload


def resolve_genus_taxon(session: requests.Session, name: str, timeout_sec: int) -> dict[str, Any]:
    payload = get_json(
        session,
        TAXA_ENDPOINT,
        {"q": name, "rank": "genus", "is_active": "true", "per_page": 30},
        timeout_sec,
    )
    for taxon in payload.get("results", []):
        if text(taxon.get("name")).lower() == name.lower() and text(taxon.get("rank")).lower() == "genus":
            return taxon
    raise RuntimeError(f"Could not resolve an active genus-level exact match for {name!r} in iNaturalist taxa API")


def preferred_photo_url(photo: dict[str, Any], size: str) -> str:
    url = text(photo.get("url") or photo.get("medium_url") or photo.get("large_url"))
    if not url:
        return ""
    import re
    return re.sub(r"/(square|small|medium|large|original)\.", f"/{size}.", url)


def parse_coordinates(observation: dict[str, Any]) -> tuple[float | None, float | None]:
    geojson = observation.get("geojson") or {}
    coordinates = geojson.get("coordinates") if isinstance(geojson, dict) else None
    if isinstance(coordinates, list) and len(coordinates) >= 2:
        lon, lat = finite_float(coordinates[0]), finite_float(coordinates[1])
        if lat is not None and lon is not None:
            return lat, lon
    location = text(observation.get("location"))
    if "," in location:
        lat_text, lon_text = location.split(",", 1)
        return finite_float(lat_text), finite_float(lon_text)
    return None, None


def annotation_summary(observation: dict[str, Any]) -> tuple[bool, str]:
    flowers = False
    parts: list[str] = []
    for annotation in observation.get("annotations") or []:
        attribute = annotation.get("controlled_attribute") or {}
        value = annotation.get("controlled_value") or {}
        attr_id = attribute.get("id")
        value_id = value.get("id")
        attr_label = text(attribute.get("label") or attr_id)
        value_label = text(value.get("label") or value_id)
        parts.append(f"{attr_label}={value_label}")
        if attr_id == 12 and value_id == 13:
            flowers = True
    return flowers, "; ".join(parts)


def observation_to_photo_rows(observation: dict[str, Any]) -> list[dict[str, Any]]:
    taxon = observation.get("taxon") or {}
    user = observation.get("user") or {}
    obs_id = observation.get("id")
    if obs_id is None:
        return []
    lat, lon = parse_coordinates(observation)
    geoprivacy = text(observation.get("geoprivacy")).lower()
    obscured = as_bool(observation.get("obscured"))
    coordinate_usable = lat is not None and lon is not None and geoprivacy not in {"obscured", "private"} and not obscured
    has_flowers, summary = annotation_summary(observation)
    rows: list[dict[str, Any]] = []
    for photo_index, photo in enumerate(observation.get("photos") or []):
        photo_id = photo.get("id")
        if photo_id is None:
            continue
        rows.append({
            "obs_id": obs_id,
            "photo_id": photo_id,
            "photo_index": photo_index,
            "photo_key": f"obs_{obs_id}_photo_{photo_id}",
            "taxon_id": taxon.get("id"),
            "taxon_name": text(taxon.get("name")),
            "taxon_rank": text(taxon.get("rank")),
            "species_guess": text(observation.get("species_guess")),
            "preferred_common_name": text(taxon.get("preferred_common_name")),
            "observed_on": text(observation.get("observed_on")),
            "created_at": text(observation.get("created_at")),
            "updated_at": text(observation.get("updated_at")),
            "quality_grade": text(observation.get("quality_grade")),
            "captive": as_bool(observation.get("captive")),
            "user_id": (user.get("id") if isinstance(user, dict) else None),
            "user_login": text(user.get("login") if isinstance(user, dict) else ""),
            "latitude": lat,
            "longitude": lon,
            "positional_accuracy": observation.get("positional_accuracy"),
            "geoprivacy": geoprivacy,
            "obscured": obscured,
            "coordinate_usable_for_environment": coordinate_usable,
            "inat_flowers_annotation": has_flowers,
            "inat_annotation_summary": summary,
            "observation_license_code": text(observation.get("license_code")),
            "photo_license_code": text(photo.get("license_code")),
            "photo_attribution": text(photo.get("attribution")),
            "raw_image_url": text(photo.get("url") or photo.get("medium_url") or photo.get("large_url")),
            "small_image_url": preferred_photo_url(photo, "small"),
            "medium_image_url": preferred_photo_url(photo, "medium"),
            "large_image_url": preferred_photo_url(photo, "large"),
            "source_file_stem": f"inat_photo_{photo_id}",
        })
    return rows


def atomic_json_dump(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    os.replace(temporary, path)


def append_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    should_write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=PHOTO_FIELDS, extrasaction="ignore")
        if should_write_header:
            writer.writeheader()
        writer.writerows(rows)


def load_seen_photo_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        return set(pd.read_csv(path, usecols=["photo_id"], dtype=str)["photo_id"].dropna().tolist())
    except ValueError:
        raise RuntimeError(f"Existing metadata file has no photo_id column: {path}")


def load_checkpoint(path: Path, restart: bool) -> dict[str, Any]:
    if restart:
        if path.exists():
            raise RuntimeError("--restart requires a fresh --out-dir. Refusing to overwrite an existing checkpoint.")
        return {"last_obs_id": 0, "observations_seen": 0, "photos_written": 0, "pages_completed": 0}
    if not path.exists():
        return {"last_obs_id": 0, "observations_seen": 0, "photos_written": 0, "pages_completed": 0}
    with path.open("r", encoding="utf-8") as handle:
        state = json.load(handle)
    return {
        "last_obs_id": int(state.get("last_obs_id", 0)),
        "observations_seen": int(state.get("observations_seen", 0)),
        "photos_written": int(state.get("photos_written", 0)),
        "pages_completed": int(state.get("pages_completed", 0)),
    }


def write_summary(out_dir: Path, state: dict[str, Any], total_results: int | None) -> None:
    summary = pd.DataFrame([{
        "collector_version": COLLECTOR_VERSION,
        "completed_at_utc": utc_now(),
        "observations_seen": state["observations_seen"],
        "photos_written": state["photos_written"],
        "last_obs_id": state["last_obs_id"],
        "pages_completed": state["pages_completed"],
        "api_total_results_on_first_page": total_results,
    }])
    summary.to_csv(out_dir / "collection_summary.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    args = parse_args()
    if not 1 <= args.per_page <= 200:
        raise ValueError("--per-page must be between 1 and 200")
    if args.max_observations < 0:
        raise ValueError("--max-observations must be >= 0")
    if args.sleep_sec < 0:
        raise ValueError("--sleep-sec must be >= 0")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = out_dir / "observation_raw.ndjson"
    metadata_path = out_dir / "photo_metadata.csv"
    checkpoint_path = out_dir / "checkpoint.json"
    provenance_path = out_dir / "collection_provenance.json"
    session = make_session(args.user_agent)

    if args.taxon_id is None:
        resolved_taxon = resolve_genus_taxon(session, args.taxon_name, args.timeout_sec)
        taxon_id = int(resolved_taxon["id"])
    else:
        resolved_taxon = {"id": args.taxon_id, "name": args.taxon_name, "rank": "genus", "resolution": "user_supplied"}
        taxon_id = int(args.taxon_id)

    state = load_checkpoint(checkpoint_path, args.restart)
    seen_photo_ids = load_seen_photo_ids(metadata_path)
    provenance = {
        "collector_version": COLLECTOR_VERSION,
        "started_at_utc": utc_now(),
        "taxon_resolution": resolved_taxon,
        "endpoints": {"taxa": TAXA_ENDPOINT, "observations": OBSERVATIONS_ENDPOINT},
        "parameters": {
            "taxon_id": taxon_id,
            "taxon_name": args.taxon_name,
            "per_page": args.per_page,
            "max_observations": args.max_observations,
            "quality_grade": args.quality_grade,
            "photo_license": args.photo_license,
            "pagination": "order_by=id; order=asc; id_above=checkpoint.last_obs_id",
            "has": "photos",
        },
        "resume_state_at_start": state,
        "preexisting_photo_ids": len(seen_photo_ids),
        "note": "Metadata collection retains obscured/private records but flags coordinate_usable_for_environment. It does not filter on flower annotations.",
    }
    atomic_json_dump(provenance, provenance_path)

    base_params: dict[str, Any] = {
        "taxon_id": taxon_id,
        "has[]": "photos",
        "per_page": args.per_page,
        "order_by": "id",
        "order": "asc",
    }
    if args.quality_grade != "any":
        base_params["quality_grade"] = args.quality_grade
    if args.photo_license != "any":
        base_params["photo_license"] = args.photo_license

    total_results: int | None = None
    print("[INFO] collecting iNaturalist metadata only")
    print(json.dumps({"out_dir": str(out_dir.resolve()), "taxon_id": taxon_id, "resume_from_obs_id": state["last_obs_id"]}, ensure_ascii=False))

    while True:
        if args.max_observations and state["observations_seen"] >= args.max_observations:
            print("[INFO] reached --max-observations")
            break
        params = dict(base_params)
        params["id_above"] = state["last_obs_id"]
        payload = get_json(session, OBSERVATIONS_ENDPOINT, params, args.timeout_sec)
        if total_results is None:
            raw_total = payload.get("total_results")
            total_results = int(raw_total) if raw_total is not None else None
        observations = payload.get("results") or []
        if not observations:
            print("[INFO] API returned no further observations")
            break

        allowed = observations
        if args.max_observations:
            remaining = args.max_observations - state["observations_seen"]
            allowed = observations[:remaining]
        rows: list[dict[str, Any]] = []
        with raw_path.open("a", encoding="utf-8") as raw_handle:
            for observation in allowed:
                obs_id = observation.get("id")
                if obs_id is None:
                    continue
                raw_handle.write(json.dumps(observation, ensure_ascii=False) + "\n")
                for row in observation_to_photo_rows(observation):
                    if str(row["photo_id"]) not in seen_photo_ids:
                        rows.append(row)
                        seen_photo_ids.add(str(row["photo_id"]))
        append_rows(metadata_path, rows)

        ids = [int(obs["id"]) for obs in allowed if obs.get("id") is not None]
        if not ids:
            raise RuntimeError("Page contained no usable observation ID; cannot safely continue id_above pagination")
        state["last_obs_id"] = max(ids)
        state["observations_seen"] += len(allowed)
        state["photos_written"] += len(rows)
        state["pages_completed"] += 1
        atomic_json_dump(state, checkpoint_path)
        write_summary(out_dir, state, total_results)
        print(f"[INFO] page={state['pages_completed']} observations={state['observations_seen']} new_photos={len(rows)} total_photos={state['photos_written']} last_obs_id={state['last_obs_id']}")

        if len(allowed) < len(observations):
            break
        time.sleep(args.sleep_sec)

    provenance["completed_at_utc"] = utc_now()
    provenance["resume_state_at_end"] = state
    provenance["api_total_results_on_first_page"] = total_results
    atomic_json_dump(provenance, provenance_path)
    write_summary(out_dir, state, total_results)
    print("[OK] metadata collection finished", out_dir.resolve())


if __name__ == "__main__":
    main()
