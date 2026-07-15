#!/usr/bin/env python3
"""Audit public molecular-data coverage for a frozen Chapter 1 taxon list.

The script queries NCBI E-utilities and writes taxon-level counts for all nucleotide
records, ITS-like records, chloroplast/plastid records and SRA studies. It does not
infer phylogeny or ploidy. Results are a time-stamped database audit and must be
reviewed for synonyms, vouchers and misidentifications before biological use.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxa", required=True, help="CSV containing a taxon_name column")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL"))
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
    parser.add_argument("--delay", type=float, default=None)
    return parser.parse_args()


def ncbi_count(db: str, term: str, email: str | None, api_key: str | None) -> int:
    params = {"db": db, "term": term, "retmode": "json", "retmax": 0, "tool": "azami_phylo_audit"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = ESEARCH + "?" + urllib.parse.urlencode(params)
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "azami-phylo-audit/1.0"})
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.load(response)
            return int(payload["esearchresult"]["count"])
        except (urllib.error.URLError, TimeoutError, KeyError, ValueError) as error:
            last_error = error
            time.sleep(2**attempt)
    raise RuntimeError(f"NCBI query failed after retries: {db} {term}") from last_error


def read_taxa(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.DictReader(handle)
        if not rows.fieldnames or "taxon_name" not in rows.fieldnames:
            raise ValueError("Input CSV must contain a taxon_name column")
        taxa = sorted({row["taxon_name"].strip() for row in rows if row.get("taxon_name", "").strip()})
    if not taxa:
        raise ValueError("No taxon names found")
    return taxa


def main() -> None:
    args = parse_args()
    taxa_path = Path(args.taxa)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    taxa = read_taxa(taxa_path)
    delay = args.delay if args.delay is not None else (0.11 if args.api_key else 0.34)

    queries = {
        "nuccore_all": ("nuccore", "{taxon}[Organism]"),
        "nuccore_its": ("nuccore", "{taxon}[Organism] AND (ITS[All Fields] OR internal transcribed spacer[All Fields])"),
        "nuccore_plastid": ("nuccore", "{taxon}[Organism] AND (chloroplast[All Fields] OR plastid[All Fields])"),
        "nuccore_complete_plastome": ("nuccore", "{taxon}[Organism] AND (chloroplast[Title] OR plastid[Title]) AND complete genome[Title]"),
        "sra_all": ("sra", "{taxon}[Organism]"),
    }

    output_rows: list[dict[str, object]] = []
    for index, taxon in enumerate(taxa, start=1):
        row: dict[str, object] = {"taxon_name": taxon}
        for column, (db, template) in queries.items():
            row[column] = ncbi_count(db, template.format(taxon=taxon), args.email, args.api_key)
            time.sleep(delay)
        output_rows.append(row)
        print(f"[{index}/{len(taxa)}] {taxon}: {row}", flush=True)

    csv_path = out_dir / "cirsium_molecular_database_coverage.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["taxon_name", *queries])
        writer.writeheader()
        writer.writerows(output_rows)

    metadata = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(taxa_path),
        "input_taxa": len(taxa),
        "service": "NCBI E-utilities esearch",
        "queries": {key: {"db": db, "template": template} for key, (db, template) in queries.items()},
        "interpretation": "Record counts are discovery metadata, not validated orthology, voucher identity, ploidy or species-tree coverage.",
    }
    (out_dir / "cirsium_molecular_database_coverage_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
