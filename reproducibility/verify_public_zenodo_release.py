#!/usr/bin/env python3
"""Verify the published Chapter 1 v2 Zenodo bundle without author credentials.

This script intentionally uses only anonymous HTTPS requests. It retrieves the
public Zenodo record metadata, downloads the published bundle from the file URL
reported by Zenodo, verifies the bundle SHA-256 and size, and then verifies the
four embedded frozen analysis inputs byte-for-byte.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

RECORD_ID = "22295791"
DOI = "10.5281/zenodo.22295791"
RECORD_API = f"https://zenodo.org/api/records/{RECORD_ID}"
RECORD_URL = f"https://zenodo.org/records/{RECORD_ID}"
BUNDLE_NAME = "azami_ch1_v2_reproduction_inputs_2026-09-04.zip"
BUNDLE_SIZE = 56_942_044
BUNDLE_SHA256 = "50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677"
EXPECTED_INPUTS = {
    "artifact-9612943217-continuous.zip": "101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e",
    "artifact-9632715852-multilevel.zip": "51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6",
    "artifact-8227254443-historical.zip": "499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc",
    "artifact-8983877726-spatial.zip": "151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca",
}
USER_AGENT = "azami-third-party-repro-verifier/1.0 (+https://github.com/zuizui0223/azami)"


def _request(url: str) -> urllib.request.addinfourl:
    # Deliberately no Authorization header, cookies, tokens, or secret links.
    req = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json, application/octet-stream;q=0.9, */*;q=0.8"},
    )
    return urllib.request.urlopen(req, timeout=120)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_zip_member(zf: zipfile.ZipFile, member: str) -> str:
    digest = hashlib.sha256()
    with zf.open(member) as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalise_files(payload: dict[str, Any]) -> list[dict[str, Any]]:
    files = payload.get("files", [])
    if isinstance(files, list):
        return [row for row in files if isinstance(row, dict)]
    if isinstance(files, dict):
        entries = files.get("entries", files)
        if isinstance(entries, list):
            return [row for row in entries if isinstance(row, dict)]
        if isinstance(entries, dict):
            out: list[dict[str, Any]] = []
            for key, value in entries.items():
                if isinstance(value, dict):
                    row = dict(value)
                    row.setdefault("key", key)
                    out.append(row)
            return out
    return []


def _file_name(row: dict[str, Any]) -> str | None:
    for key in ("key", "filename", "name"):
        value = row.get(key)
        if isinstance(value, str) and value:
            return value
    return None


def _file_size(row: dict[str, Any]) -> int | None:
    for key in ("size", "filesize"):
        value = row.get(key)
        if isinstance(value, int):
            return value
    return None


def _download_url(row: dict[str, Any]) -> str:
    links = row.get("links", {})
    if isinstance(links, dict):
        for key in ("content", "download", "self"):
            value = links.get(key)
            if isinstance(value, str) and value.startswith("https://"):
                return value
    return f"{RECORD_URL}/files/{BUNDLE_NAME}?download=1"


def _assert_public_access(payload: dict[str, Any]) -> dict[str, Any]:
    access = payload.get("access", {})
    if not isinstance(access, dict):
        return {}
    observed: dict[str, Any] = {}
    for key in ("record", "files", "status"):
        if key in access:
            observed[key] = access[key]
    # New Zenodo commonly reports record/files as 'public'. Older records may
    # expose an access_right/access.status of 'open'. Reject only explicit
    # restricted/closed/embargoed values; the successful anonymous file GET is
    # the decisive access test.
    forbidden = {"restricted", "closed", "embargoed", "private"}
    for value in observed.values():
        if isinstance(value, str) and value.lower() in forbidden:
            raise RuntimeError(f"Zenodo record is not credential-free public: access={observed}")
    return observed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("_zenodo_public_verification"))
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    report_path = (args.report or (out / "zenodo_public_verification.json")).resolve()
    bundle_path = out / BUNDLE_NAME
    extracted = out / "reproduction_inputs"
    extracted.mkdir(parents=True, exist_ok=True)

    with _request(RECORD_API) as response:
        if getattr(response, "status", 200) != 200:
            raise RuntimeError(f"Anonymous Zenodo metadata request failed: HTTP {response.status}")
        payload = json.load(response)

    observed_doi = payload.get("doi")
    if observed_doi != DOI:
        # Some API representations put identifiers under pids.
        pids = payload.get("pids", {})
        doi_value = None
        if isinstance(pids, dict):
            doi_entry = pids.get("doi", {})
            if isinstance(doi_entry, dict):
                doi_value = doi_entry.get("identifier")
        if doi_value != DOI:
            raise RuntimeError(f"Unexpected DOI: payload doi={observed_doi!r}, pids doi={doi_value!r}")

    access = _assert_public_access(payload)
    files = _normalise_files(payload)
    matching = [row for row in files if _file_name(row) == BUNDLE_NAME]
    if len(matching) != 1:
        raise RuntimeError(f"Expected exactly one public bundle named {BUNDLE_NAME!r}; found {len(matching)}")
    entry = matching[0]
    reported_size = _file_size(entry)
    if reported_size is not None and reported_size != BUNDLE_SIZE:
        raise RuntimeError(f"Zenodo metadata size mismatch: {reported_size} != {BUNDLE_SIZE}")

    download_url = _download_url(entry)
    with _request(download_url) as response, bundle_path.open("wb") as fh:
        if getattr(response, "status", 200) != 200:
            raise RuntimeError(f"Anonymous Zenodo file request failed: HTTP {response.status}")
        shutil.copyfileobj(response, fh, length=8 * 1024 * 1024)

    observed_size = bundle_path.stat().st_size
    observed_bundle_sha = _sha256_file(bundle_path)
    if observed_size != BUNDLE_SIZE:
        raise RuntimeError(f"Downloaded bundle size mismatch: {observed_size} != {BUNDLE_SIZE}")
    if observed_bundle_sha != BUNDLE_SHA256:
        raise RuntimeError(f"Downloaded bundle SHA-256 mismatch: {observed_bundle_sha} != {BUNDLE_SHA256}")

    embedded: dict[str, dict[str, Any]] = {}
    with zipfile.ZipFile(bundle_path) as zf:
        names = set(zf.namelist())
        for filename, expected_sha in EXPECTED_INPUTS.items():
            if filename not in names:
                raise RuntimeError(f"Missing embedded input: {filename}")
            observed_sha = _sha256_zip_member(zf, filename)
            if observed_sha != expected_sha:
                raise RuntimeError(f"Embedded SHA-256 mismatch for {filename}: {observed_sha} != {expected_sha}")
            target = extracted / filename
            with zf.open(filename) as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
            embedded[filename] = {
                "sha256": observed_sha,
                "size_bytes": target.stat().st_size,
                "status": "PASS",
            }

    report = {
        "schema_version": 1,
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "verification_mode": "anonymous_credential_free_https",
        "authorization_header_used": False,
        "record_id": RECORD_ID,
        "doi": DOI,
        "record_url": RECORD_URL,
        "record_api": RECORD_API,
        "access_metadata": access,
        "public_file_download_url": download_url,
        "bundle": {
            "filename": BUNDLE_NAME,
            "size_bytes": observed_size,
            "sha256": observed_bundle_sha,
            "status": "PASS",
        },
        "embedded_inputs": embedded,
        "status": "PASS",
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"PUBLIC ZENODO VERIFICATION FAILED: {exc}", file=sys.stderr)
        raise
