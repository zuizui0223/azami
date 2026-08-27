#!/usr/bin/env python3
"""Validate the submission-facing v2 hypothesis recovery ledger."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS_PATH = ROOT / "analysis" / "ch1" / "hypothesis_recovery_status.json"

EXPECTED_IDS = {
    "H1_CONTINUOUS_EXTRACTION",
    "H2_TAXON_MEAN_INFORMATION_LOSS",
    "H3_ENDPOINT_GEOGRAPHY",
    "H4_SCALE_INVARIANCE",
    "H5_SINGLE_CAPITULUM_SYNDROME",
    "H6_ADAPTIVE_PATTERN_CANDIDATES",
}


def main() -> int:
    data = json.loads(STATUS_PATH.read_text(encoding="utf-8"))
    assert data["schema_version"] == 2
    assert data["release_label"] == "Azami Chapter 1 v2"
    assert data["status"] == "complete_for_present_geography_submission_lane"

    rows = data["questions"]
    assert {row["id"] for row in rows} == EXPECTED_IDS
    assert len(rows) == len(EXPECTED_IDS)
    by_id = {row["id"]: row for row in rows}
    assert by_id["H1_CONTINUOUS_EXTRACTION"]["status"] == "supported"
    assert by_id["H2_TAXON_MEAN_INFORMATION_LOSS"]["status"] == "supported"
    assert by_id["H3_ENDPOINT_GEOGRAPHY"]["status"] == "supported"
    assert by_id["H4_SCALE_INVARIANCE"]["status"] == "rejected"
    assert by_id["H5_SINGLE_CAPITULUM_SYNDROME"]["status"] == "rejected"
    assert by_id["H6_ADAPTIVE_PATTERN_CANDIDATES"]["status"] == "two_candidates_under_current_controls"

    series = data["eazami_boundary"]["series"]
    assert len(series) == 4
    assert [row.split()[0] for row in series] == ["EAzami-I", "EAzami-II", "EAzami-III", "EAzami-IV"]
    assert "cannot be promoted" in data["eazami_boundary"]["rule"]
    assert "legacy lability-responsiveness and quadrant hypothesis" in data["excluded_from_active_submission"]

    print(json.dumps({
        "status": "ok",
        "questions": len(rows),
        "supported": 3,
        "rejected": 2,
        "candidate_gate": by_id["H6_ADAPTIVE_PATTERN_CANDIDATES"]["status"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
