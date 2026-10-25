"""Validation contract for the Chapter 1 public-release metadata.

This module deliberately validates only information that must be settled before
`build_current_release_bundle --final` may claim a release-ready package. It
does not choose author order, affiliations, licensing strategy or archive
strategy for the authors.
"""
from __future__ import annotations

import json
from pathlib import Path

V2_RECORD_DOI = "10.5281/zenodo.22295791"
V2_CONCEPT_DOI = "10.5281/zenodo.22295790"
SCOPE_ID = "azami_ch1_v3_current"
ARCHIVE_STRATEGIES = {"new_version_existing_concept", "linked_code_data_records"}
PLACEHOLDER_TOKENS = {"", "tbd", "todo", "unknown", "placeholder", "none", "null"}


def _text(value, field: str) -> str:
    if not isinstance(value, str) or value.strip().lower() in PLACEHOLDER_TOKENS:
        raise ValueError(f"release metadata field {field!r} must be final non-placeholder text")
    return value.strip()


def validate_release_metadata(path: Path, *, expected_head: str) -> dict:
    """Return validated metadata or fail closed.

    The JSON is an internal release contract, not a direct Zenodo REST payload.
    It keeps author-owned decisions explicit while binding the release to the
    exact code commit and the preserved v2 record.
    """
    obj = json.loads(path.read_text(encoding="utf-8"))
    if obj.get("schema_version") != 1:
        raise ValueError("release metadata schema_version must equal 1")
    if obj.get("release_approved") is not True:
        raise ValueError("release metadata must set release_approved=true")
    if obj.get("scientific_scope") != SCOPE_ID:
        raise ValueError(f"scientific_scope must equal {SCOPE_ID}")

    head = _text(obj.get("final_code_commit"), "final_code_commit")
    if head != expected_head:
        raise ValueError(f"release metadata final_code_commit {head} does not match HEAD {expected_head}")

    _text(obj.get("title"), "title")
    _text(obj.get("description"), "description")

    creators = obj.get("creators")
    if not isinstance(creators, list) or not creators:
        raise ValueError("release metadata creators must be a non-empty ordered list")
    for i, creator in enumerate(creators):
        if not isinstance(creator, dict):
            raise ValueError(f"creator {i} must be an object")
        _text(creator.get("name"), f"creators[{i}].name")
        if "orcid" in creator and creator["orcid"] not in (None, ""):
            _text(creator["orcid"], f"creators[{i}].orcid")
        if "affiliation" in creator and creator["affiliation"] not in (None, ""):
            _text(creator["affiliation"], f"creators[{i}].affiliation")

    strategy = obj.get("archive_strategy")
    if strategy not in ARCHIVE_STRATEGIES:
        raise ValueError(f"archive_strategy must be one of {sorted(ARCHIVE_STRATEGIES)}")

    preserved = obj.get("preserved_v2_record")
    if not isinstance(preserved, dict):
        raise ValueError("preserved_v2_record must be an object")
    if preserved.get("record_doi") != V2_RECORD_DOI or preserved.get("concept_doi") != V2_CONCEPT_DOI:
        raise ValueError("release metadata must preserve the exact published v2 record/concept DOI")
    if preserved.get("modify_existing_v2") is not False:
        raise ValueError("modify_existing_v2 must be false")

    licensing = obj.get("licensing")
    if not isinstance(licensing, dict):
        raise ValueError("licensing must be an object")
    _text(licensing.get("software"), "licensing.software")
    _text(licensing.get("data_and_third_party_strategy"), "licensing.data_and_third_party_strategy")

    if obj.get("anonymous_redownload_required") is not True:
        raise ValueError("anonymous_redownload_required must be true")
    if obj.get("clean_replay_required") is not True:
        raise ValueError("clean_replay_required must be true")

    return obj
