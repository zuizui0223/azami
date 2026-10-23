"""Fail-closed, small-range audit of Ren 2024 published Senecioneae nuclear ASTRAL tree.

The Dryad ZIP is 1.58 GB; this code ONLY requests small HTTP 206 byte ranges.
HTTP 200 (server ignoring Range), source SHA/size mismatch, or excess transfer
causes a STOP without full archive download. A matched label is not a verified
taxon concept, specimen, inferred evolutionary event, or rain adaptation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

import requests
from Bio import Phylo

DOI = "10.5061/dryad.573n5tbf7"
VERSION_ID = 340361
ARCHIVE_ID = 3337673
ARCHIVE_SHA256 = "962495cf8131fadd3a1cdf1ea88204432f2aaa81c9e6c0a97435756f0ad88139"
ARCHIVE_SIZE = 1583779312
TREE_SUFFIX = "1_G727_hyb/2_astral_analysis/astral.tre"
MAX_TOTAL_RANGE_BYTES = 32 * 1024 * 1024
MAX_REMOTE_REQUEST_BYTES = 2 * 1024 * 1024
MAX_TREE_MEMBER_UNCOMPRESSED = 3 * 1024 * 1024
MAX_RANGE_CALLS = 256
META_URL = f"https://datadryad.org/api/v2/versions/{VERSION_ID}/files"
RANGE_URLS = [
    f"https://datadryad.org/api/v2/files/{ARCHIVE_ID}/download",
    f"https://datadryad.org/downloads/file_stream/{ARCHIVE_ID}",
]


class StrictRangeError(RuntimeError):
    """Fail-closed refusal of a non-ranged or unverified remote archive."""


def check_official_manifest(session, *, timeout=25):
    response = session.get(META_URL, timeout=timeout, headers={"Accept": "application/json"})
    response.raise_for_status()
    manifest = response.json()
    records = manifest.get("_embedded", {}).get("stash:files", [])
    matches = [x for x in records
               if str(x.get("_links", {}).get("self", {}).get("href", "")).endswith("/" + str(ARCHIVE_ID))]
    if len(matches) != 1:
        raise ValueError("Dryad published version has no unique pinned source ZIP")
    record = matches[0]
    if (int(record.get("size", -1)) != ARCHIVE_SIZE
            or record.get("digest") != ARCHIVE_SHA256
            or record.get("digestType") != "sha-256"):
        raise ValueError("Dryad source ZIP size/SHA256 differs from source-pinned publication")
    return {"archive_id": ARCHIVE_ID, "archive_bytes": ARCHIVE_SIZE,
            "published_zip_sha256": ARCHIVE_SHA256,
            "published_sha_verified_against_live_dryad_metadata": True,
            "metadata_endpoint": META_URL,
            "manifest_member_count": len(records)}


class StrictHTTPRangeReader(io.RawIOBase):
    """Bounded seekable HTTP reader requiring exact Content-Range for every byte."""
    def __init__(self, session, url, size, *, budget=MAX_TOTAL_RANGE_BYTES,
                 max_request=MAX_REMOTE_REQUEST_BYTES):
        self._session = session
        self.url = url
        self._size = size
        self._pos = 0
        self._budget = budget
        self._max_request = max_request
        self.bytes_received = 0
        self.range_requests = 0

    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self._pos

    def seek(self, offset, whence=io.SEEK_SET):
        if whence == io.SEEK_SET: position = offset
        elif whence == io.SEEK_CUR: position = self._pos + offset
        elif whence == io.SEEK_END: position = self._size + offset
        else: raise ValueError("Bad seek mode")
        if not 0 <= position <= self._size:
            raise StrictRangeError("Seek outside published ZIP size")
        self._pos = int(position)
        return self._pos

    def _fetch_exact(self, start, end):
        n = end - start + 1
        if (n <= 0 or n > self._max_request
                or self.bytes_received + n > self._budget
                or self.range_requests >= MAX_RANGE_CALLS):
            raise StrictRangeError("Bounded partial-ZIP request budget exceeded")
        headers = {"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity",
                   "User-Agent": "azami-T3-published-tree-audit/1.0"}
        response = self._session.get(self.url, headers=headers,
                                     stream=True, timeout=(10, 30), allow_redirects=True)
        try:
            if response.status_code != 206:
                raise StrictRangeError(
                    f"Range required: source returned HTTP {response.status_code}; no full download permitted")
            actual = response.headers.get("Content-Range", "")
            expected = f"bytes {start}-{end}/{self._size}"
            if actual != expected:
                raise StrictRangeError("Content-Range does not match pinned byte offsets/size: " + actual[:120])
            if response.headers.get("Content-Encoding", "identity").lower() != "identity":
                raise StrictRangeError("HTTP entity encoded; ZIP byte ranges not trustworthy")
            output = bytearray()
            for chunk in response.iter_content(chunk_size=262144):
                output.extend(chunk)
                if len(output) > n:
                    raise StrictRangeError("Returned more than authorized partial-ZIP bytes")
            if len(output) != n:
                raise StrictRangeError("Returned partial-ZIP body truncated")
            self.bytes_received += n
            self.range_requests += 1
            return bytes(output)
        finally:
            response.close()

    def read(self, size=-1):
        if self.closed: raise ValueError("Reader closed")
        if size is None or size < 0: size = self._size - self._pos
        size = min(size, self._size - self._pos)
        if size == 0: return b""
        chunks = []
        while size:
            chunk_size = min(size, self._max_request)
            chunks.append(self._fetch_exact(self._pos, self._pos + chunk_size - 1))
            self._pos += chunk_size
            size -= chunk_size
        return b"".join(chunks)


def zip_member_tree(session, source_url, *, expected_bytes=ARCHIVE_SIZE):
    stream = StrictHTTPRangeReader(session, source_url, expected_bytes)
    with zipfile.ZipFile(stream, "r") as archive:
        names = [n for n in archive.namelist()
                 if n.replace("\\", "/").endswith(TREE_SUFFIX)]
        if len(names) != 1:
            raise ValueError(f"Missing/ambiguous genuine ASTRAL tree ZIP member: {len(names)} matches")
        info = archive.getinfo(names[0])
        if info.file_size > MAX_TREE_MEMBER_UNCOMPRESSED:
            raise StrictRangeError("Unbounded ASTRAL member size; refuse extraction")
        newick = archive.read(info)
        if len(newick) != info.file_size:
            raise ValueError("Decompressed Newick member has invalid length")
        n_members = len(archive.namelist())
    if not newick.strip().endswith(b";"):
        raise ValueError("Original published ASTRAL tree not Newick-like")
    return {"zip_internal_path": names[0], "zip_total_listed_members": n_members,
            "newick": newick.decode("utf-8-sig"),
            "newick_SHA256": hashlib.sha256(newick).hexdigest(),
            "zip_bytes_downloaded_via_206": stream.bytes_received,
            "zip_range_requests": stream.range_requests,
            "remote_archive_bytes": expected_bytes,
            "downloaded_entire_original_zip": False,
            "http_source": source_url}


def tip_binomial_matches(tree_text, candidates):
    try:
        tree = Phylo.read(io.StringIO(tree_text), "newick")
    except Exception as exc:
        raise ValueError("Published ASTRAL Newick could not be parsed") from exc
    tips = [x.name for x in tree.get_terminals()]
    if any(not t for t in tips) or len(set(tips)) != len(tips):
        raise ValueError("Unnamed or duplicate original nuclear tree tips")
    matches = []
    for item in candidates:
        name = item["taxon"]
        parts = name.split()[:2]
        if len(parts) != 2:
            raise ValueError("Non-binomial T3 source name")
        pattern = re.compile(r"(?<![A-Za-z])" + re.escape(parts[0])
                             + r"(?:[_\s.-]+)" + re.escape(parts[1]) + r"(?![A-Za-z])", re.I)
        labels = [tip for tip in tips if pattern.search(tip)]
        matches.append({"source_name": name, "source_state": item["source_posture"],
                        "taxonomic_status": item["taxonomic_status"],
                        "matched_original_tip_labels": labels,
                        "binomial_string_match_count": len(labels),
                        "biological_taxon_concept_or_voucher_verified": False})
    return {"published_astral_tree_tip_count": len(tips),
            "published_original_tree_tips": tips,
            "candidate_binomial_matches": matches,
            "n_candidate_source_names_with_exact_binomial_string_match":
                sum(x["binomial_string_match_count"] > 0 for x in matches),
            "sufficient_tree_and_traits_to_infer_independent_evolutionary_gains": False,
            "head_orientation_true_gravity_and_anthesis_state_tied_to_voucher": False,
            "rooted_tribal_species_tree_provenance_verified": False,
            "biological_rain_protection_causality_inferred": False}


def audit(session, candidates, urls=RANGE_URLS):
    metadata = check_official_manifest(session)
    failures = []
    for url in urls:
        try:
            details = zip_member_tree(session, url)
            report = tip_binomial_matches(details.pop("newick"), candidates)
            return {"status": "EMPIRICAL_SOURCE_NUCLEAR_TIPS_VERIFIED_STRICT_RANGE",
                    "source": metadata, "tree": details,
                    "taxon_comparison": report,
                    "failure_attempts_prior_to_success": failures}
        except (requests.RequestException, OSError, ValueError, StrictRangeError,
                zipfile.BadZipFile, UnicodeError) as exc:
            failures.append({"url": url, "error_class": type(exc).__name__,
                             "message": str(exc)[:220]})
    return {"status": "NO_AUTHENTIC_RANGE_ASTRAL_TREE_EXTRACTED_STOP",
            "source": metadata, "retrieval_attempts": failures,
            "no_tree_tip_overlap_claim": True,
            "no_independent_nodding_origin_claim": True}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-table", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    opts = p.parse_args()
    from analysis.v4.t3_source_evidence_gate import read_evidence
    rows = read_evidence(opts.candidate_table)
    try:
        result = audit(requests.Session(), rows)
    except Exception as exc:
        result = {"status": "PUBLIC_SOURCE_METADATA_UNAVAILABLE_FAIL_CLOSED",
                  "exception": type(exc).__name__, "reason": str(exc)[:300],
                  "unverified_original_tree_tips_are_not_biological_evidence": True}
    opts.output.parent.mkdir(parents=True, exist_ok=True)
    opts.output.write_text(json.dumps(result, indent=2, sort_keys=True,
                                      ensure_ascii=False, allow_nan=False) + "\n",
                           encoding="utf-8")
    print("T3_DRYAD727_AST_TREE_SOURCE_AUDIT=" +
          json.dumps(result, ensure_ascii=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
