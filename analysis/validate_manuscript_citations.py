#!/usr/bin/env python3
"""Fail when active manuscript citations and references diverge."""
from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BODY = [
    ROOT / "manuscript" / "00_title_abstract.md",
    ROOT / "manuscript" / "01_introduction.md",
    ROOT / "manuscript" / "02_methods.md",
    ROOT / "manuscript" / "03_results.md",
    ROOT / "manuscript" / "04_discussion.md",
    ROOT / "manuscript" / "05_conclusion_and_declarations.md",
]
REFERENCES = ROOT / "manuscript" / "06_references.md"


def normalize(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower())


def first_author(fragment: str) -> str:
    fragment = re.sub(r"\bet\s+al\.?\s*$", "", fragment, flags=re.IGNORECASE).strip()
    fragment = fragment.split("&", 1)[0].strip()
    return normalize(fragment.rstrip(".,"))


def expanded_years(fragment: str) -> set[str]:
    expanded = re.sub(r"(\d{4})([a-z]),([a-z])", r"\1\2, \1\3", fragment)
    return set(re.findall(r"(?:19|20)\d{2}[a-z]?", expanded))


def main() -> None:
    reference_keys: set[tuple[str, str]] = set()
    reference_rows = 0
    for line in REFERENCES.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^([^,]+),.*\(((?:19|20)\d{2}[a-z]?)\)\.", line)
        if match:
            reference_keys.add((normalize(match.group(1)), match.group(2)))
            reference_rows += 1

    body = "\n".join(path.read_text(encoding="utf-8") for path in BODY)
    citation_keys: set[tuple[str, str]] = set()
    for group in re.findall(r"\(([^()]*(?:19|20)\d{2}[^()]*)\)", body):
        for item in group.split(";"):
            year_match = re.search(r"(?:19|20)\d{2}", item)
            if not year_match:
                continue
            author = first_author(item[: year_match.start()])
            for year in expanded_years(item[year_match.start() :]):
                citation_keys.add((author, year))

    missing = sorted(citation_keys - reference_keys)
    unused = sorted(reference_keys - citation_keys)
    if missing or unused:
        raise AssertionError(json.dumps({"missing_references": missing, "unused_references": unused}, indent=2))
    print(json.dumps({
        "status": "PASS",
        "reference_rows": reference_rows,
        "unique_author_year_citations": len(citation_keys),
        "missing_references": 0,
        "unused_references": 0,
    }, indent=2))


if __name__ == "__main__":
    main()
