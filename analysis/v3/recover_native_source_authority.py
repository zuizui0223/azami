"""Replayable source-only recovery of WCVP/TDWG authority, without v2 traits.

Uses an isolated copy of the unchanged historical classification helpers. Raw
HTTP responses are retained privately, hash-checked on reuse, and bound into a
manifest for offline replay. Reacquisition is not historical-response recovery;
only a later exact cohort hash comparison can establish the old cohort's bytes.
"""
from __future__ import annotations

import argparse
import base64
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import threading
import time

import pandas as pd
import requests

from .build_native_range_join import classify_frame, load_contract
from .build_native_range_join_from_source import read_source, normalized_state
from .enriched_source_cohort import pinned_hash
from .workflow import ROOT, digest, text_digest


def private_directory(path: Path) -> Path:
    path = path.resolve()
    if path.is_relative_to(ROOT) and (not path.relative_to(ROOT).parts or path.relative_to(ROOT).parts[0] != "local_data"):
        raise ValueError("Private source output inside the repository must use local_data")
    return path


def write_new_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def recover_lf_source(source: Path, destination: Path, expected_sha256: str) -> dict:
    """Make a new LF byte view only if it exactly reproduces the historical pin.

    Original Windows output is never overwritten. This changes serialization,
    not rows or data rules; mismatch fails before a candidate is written.
    """
    private_directory(destination.parent)
    original = source.read_bytes()
    candidate = original.replace(b"\r\n", b"\n")
    actual = hashlib.sha256(candidate).hexdigest()
    if actual != expected_sha256:
        raise ValueError("LF source candidate does not reproduce the historical source hash")
    with destination.open("xb") as handle:
        handle.write(candidate)
    report = {"input_sha256": hashlib.sha256(original).hexdigest(),
              "output_sha256": actual, "crlf_pairs_replaced": original.count(b"\r\n"),
              "historical_source_bytes_verified": True, "original_preserved": True}
    write_new_json(destination.with_suffix(".recovery.json"), report)
    return report


class ResponseCache:
    def __init__(self, directory: Path, *, offline: bool = False,
                 expected_manifest_sha256: str | None = None):
        self.directory = private_directory(directory)
        self.offline = offline
        self.used: dict[str, str] = {}
        self.lock = threading.Lock()
        self.frozen: dict[str, str] | None = None
        manifest = self.directory / "authority_cache_manifest.json"
        if offline:
            pinned_hash(manifest, expected_manifest_sha256 or "", "authority cache manifest")
        if manifest.exists():
            frozen = json.loads(manifest.read_text(encoding="utf-8"))
            if frozen.get("schema_version") != 1:
                raise ValueError("Unknown authority cache manifest schema")
            self.frozen = frozen["entries"]
        elif not offline:
            self.directory.mkdir(parents=True, exist_ok=True)

    def get(self, url: str, timeout: float = 30, retries: int = 4, *, slot: str = "") -> bytes:
        key = hashlib.sha256((slot + "\n" + url).encode()).hexdigest() + ".json"
        path = self.directory / key
        if self.frozen is not None:
            if key not in self.frozen or not path.exists() or digest(path) != self.frozen[key]:
                raise ValueError("Frozen authority response is missing or changed")
        if not path.exists():
            if self.offline:
                raise ValueError("Offline replay lacks a required authority response")
            if retries < 1 or timeout <= 0:
                raise ValueError("Positive request timeout and retries required")
            for attempt in range(retries):
                try:
                    response = requests.get(url, timeout=timeout, headers={
                        "Accept": "application/json", "User-Agent": "azami-v3-source-recovery/1.0"})
                    response.raise_for_status()
                    content = response.content
                    if not isinstance(json.loads(content), dict):
                        raise ValueError("Authority response must be a JSON object")
                    break
                except (requests.RequestException, ValueError):
                    if attempt + 1 == retries:
                        raise
                    time.sleep(min(2 ** attempt, 8))
            write_new_json(path, {"url": url, "slot": slot,
                "acquired_at_utc": datetime.now(timezone.utc).isoformat(),
                "http_status": response.status_code,
                "response_sha256": hashlib.sha256(content).hexdigest(),
                "content_base64": base64.b64encode(content).decode("ascii")})
        envelope = json.loads(path.read_text(encoding="utf-8"))
        content = base64.b64decode(envelope["content_base64"], validate=True)
        if (envelope["url"] != url or envelope["slot"] != slot
                or envelope["http_status"] != 200
                or hashlib.sha256(content).hexdigest() != envelope["response_sha256"]):
            raise ValueError("Authority response envelope identity mismatch")
        if not isinstance(json.loads(content), dict):
            raise ValueError("Cached authority response must be a JSON object")
        with self.lock:
            self.used[key] = digest(path)
        return content

    def json(self, url: str, timeout: float = 30, retries: int = 4) -> dict:
        return json.loads(self.get(url, timeout, retries))

    def seal(self) -> str:
        if self.frozen is not None and self.used != self.frozen:
            raise ValueError("Authority replay did not consume the exact frozen request set")
        manifest = self.directory / "authority_cache_manifest.json"
        if not manifest.exists():
            write_new_json(manifest, {"schema_version": 1, "entries": self.used,
                "evidence_kind": "reacquired_responses_not_historical_response_identity"})
        return digest(manifest)


def historical_helpers(cache: ResponseCache):
    """Inject only HTTP acquisition into an isolated unchanged helper module."""
    path = ROOT / "analysis" / "rebuild_frozen_native_status.py"
    spec = importlib.util.spec_from_file_location("azami_native_recovery_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load historical native classification helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.request_json = cache.json
    return module


def verify_metadata(metadata: dict, source: dict) -> None:
    if str(metadata.get("doi", "")).casefold() != source["dataset_doi"].casefold():
        raise ValueError("WCVP dataset DOI differs from the contract")
    if metadata.get("modified") != source["dataset_modified"]:
        raise ValueError("WCVP dataset modification timestamp differs from the contract")


def run(source_csv: Path, out_dir: Path, cache_dir: Path, *, expected_source_sha256: str,
        contract_path: Path, offline: bool = False, expected_cache_manifest_sha256: str | None = None,
        workers: int = 4, timeout: float = 30, retries: int = 4) -> dict:
    if not 1 <= workers <= 4:
        raise ValueError("Use one to four workers for bounded authority acquisition")
    out_dir = private_directory(out_dir)
    if out_dir.exists():
        raise ValueError("Output exists; preserve previous recovery or replay")
    source_hash = pinned_hash(source_csv, expected_source_sha256, "source observations")
    contract = load_contract(contract_path)
    frame = read_source(source_csv)
    source = contract["source_taxonomy_and_distribution"]
    cache = ResponseCache(cache_dir, offline=offline,
                          expected_manifest_sha256=expected_cache_manifest_sha256)
    helpers = historical_helpers(cache)
    metadata_url = f"https://api.gbif.org/v1/dataset/{source['dataset_key']}"
    before = json.loads(cache.get(metadata_url, timeout, retries, slot="before"))
    verify_metadata(before, source)
    names = sorted({helpers.normalized(n) for n in frame["source_taxon_name"].dropna()
                    if helpers.normalized(n)})

    def resolve(name):
        return helpers.fetch_name_resolution([name], source["dataset_key"], timeout, retries, 0)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        resolution = pd.concat(list(pool.map(resolve, names)), ignore_index=True)
    print(json.dumps({"stage": "names_recovered", "source_names": len(names)}), flush=True)
    keys = sorted({int(k) for k in resolution["accepted_key"] if helpers.normalized(k)})

    def distribution(key):
        return helpers.fetch_distributions(pd.DataFrame({"accepted_key": [key]}), timeout, retries, 0)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(distribution, keys))
    distributions = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    print(json.dumps({"stage": "distributions_recovered", "accepted_keys": len(keys)}), flush=True)
    geo = contract["tdwg_level3"]
    geometry_url = f"https://raw.githubusercontent.com/tdwg/wgsrpd/{geo['commit']}/{geo['path']}"
    geometry = cache.get(geometry_url, timeout, retries)
    if hashlib.sha256(geometry).hexdigest() != geo["expected_sha256"]:
        raise ValueError("Pinned TDWG geometry differs from the contract")
    after = json.loads(cache.get(metadata_url, timeout, retries, slot="after"))
    verify_metadata(after, source)
    manifest_sha = cache.seal()
    joined = classify_frame(frame, resolution, distributions, json.loads(geometry))
    primary = (joined["native_range_status"].eq("native")
        & normalized_state(joined["captive_state"]).eq("false")
        & joined["date_status"].eq("exact_day")
        & joined["coordinate_status"].eq("public_location_present_precision_not_gated")
        & joined["taxon_resolution_status"].eq("resolved_unique_accepted_key"))
    out_dir.mkdir(parents=True)
    tables = {"native_range_join.csv": joined, "taxon_resolution.csv": resolution,
              "wcvp_distribution_records.csv": distributions}
    for name, table in tables.items():
        table.to_csv(out_dir / name, index=False, lineterminator="\n")
    implementations = [Path(__file__), ROOT / "analysis/rebuild_frozen_native_status.py",
                       ROOT / "analysis/v3/build_native_range_join.py",
                       ROOT / "analysis/v3/build_native_range_join_from_source.py"]
    report = {"status": "NATIVE_AUTHORITY_REPLAYED" if offline else "NATIVE_AUTHORITY_REACQUIRED",
        "input_source_sha256": source_hash, "contract_sha256": digest(contract_path),
        "authority_cache_manifest_sha256": manifest_sha, "authority_responses": len(cache.used),
        "n_observations": len(joined), "n_source_taxa": len(names),
        "status_counts": {str(k): int(v) for k, v in joined["native_range_status"].value_counts().items()},
        "primary_native_wild_public_exact_date_resolved_rows": int(primary.sum()),
        "outputs_sha256": {name: digest(out_dir / name) for name in tables},
        "implementations_sha256_text_lf": {str(p.relative_to(ROOT)).replace('\\', '/'): text_digest(p) for p in implementations},
        "runtime": {"python": platform.python_version(), "pandas": pd.__version__},
        "source_rows_deleted": 0, "trait_files_read": 0, "ecological_models_executed": 0,
        "ecological_fitting_authorized": False, "historical_cohort_identity_verified": False,
        "limits": ["Historical API response bytes were not archived; these are reacquired responses.",
                   "The pinned dataset metadata and geometry do not alone prove historical cohort identity.",
                   "The unchanged historical name helper uses its first 100 search records; this is a recovery route, not a new comprehensive taxonomic audit.",
                   "Native status is regional, and source identification is not independently validated."]}
    write_new_json(out_dir / "native_authority_recovery_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--contract-path", type=Path, default=ROOT / "analysis/v3/native_range_join_contract.json")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--expected-cache-manifest-sha256")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    print(json.dumps(run(**vars(args)), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
