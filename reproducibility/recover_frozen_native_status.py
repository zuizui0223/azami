#!/usr/bin/env python3
"""Recover the exact frozen native-status CSV bytes from public Git history.

Git stores the historical CSV with LF line endings, while the frozen sampling
report hashed the original artifact written with CRLF line endings. The table
content is identical. This helper validates the historical Git blob identity,
then restores CRLF deterministically and requires the frozen artifact SHA-256.
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
from pathlib import Path

TAG = "azami-ch1-v2-2026-08-27"
PATH = "analysis_outputs/native_range_sensitivity_v1/observation_native_status.csv"
EXPECTED_GIT_BLOB_SHA = "b98af47482fd86b1353546573492519659cda848"
EXPECTED_LF_SHA256 = "9686b8f515deef3b3aa9311b5137317a094af195e0f2414f2b1b70d7c72b5021"
EXPECTED_CRLF_SHA256 = "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a"
EXPECTED_ROWS_WITH_HEADER = 46277


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    spec = f"{TAG}:{PATH}"
    blob_sha = subprocess.check_output(["git", "rev-parse", spec], text=True).strip()
    if blob_sha != EXPECTED_GIT_BLOB_SHA:
        raise SystemExit(f"Historical native-status Git blob changed: {blob_sha}; expected {EXPECTED_GIT_BLOB_SHA}")

    raw = subprocess.check_output(["git", "cat-file", "blob", blob_sha])
    raw_sha = sha256(raw)
    if raw_sha != EXPECTED_LF_SHA256:
        raise SystemExit(f"Historical LF-native-status SHA-256 mismatch: {raw_sha}; expected {EXPECTED_LF_SHA256}")
    if b"\r\n" in raw:
        raise SystemExit("Historical Git blob unexpectedly already contains CRLF bytes")
    if raw.count(b"\n") != EXPECTED_ROWS_WITH_HEADER:
        raise SystemExit(
            f"Historical native-status line-count mismatch: {raw.count(b'\n')}; expected {EXPECTED_ROWS_WITH_HEADER}"
        )

    exact = raw.replace(b"\n", b"\r\n")
    exact_sha = sha256(exact)
    if exact_sha != EXPECTED_CRLF_SHA256:
        raise SystemExit(f"Recovered CRLF-native-status SHA-256 mismatch: {exact_sha}; expected {EXPECTED_CRLF_SHA256}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(exact)
    print(f"PASS {args.out} sha256={exact_sha} bytes={len(exact)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
