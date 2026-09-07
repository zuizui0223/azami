import json
from pathlib import Path

import pandas as pd

from analysis.v3.build_environment_design_cohort import build_design_cohort


ROOT = Path(__file__).resolve().parents[1]


def _eligible_rows(key: int, n: int, status: str = "native") -> list[dict]:
    rows = []
    for i in range(n):
        rows.append({
            "obs_id": f"{key}-{i}",
            "accepted_key": key,
            "accepted_name": f"Taxon {key}",
            "source_taxon_name": f"Taxon {key}",
            "taxon_resolution_status": "resolved_unique_accepted_key",
            "native_range_status": status,
            "captive_state": "false",
            "date_status": "exact_day",
            "coordinate_status": "public_location_present_precision_not_gated",
            "analysis_latitude": 30 + i / 1000,
            "analysis_longitude": 130 + i / 1000,
            "observation_month": 1 + (i % 12),
        })
    return rows


def test_source_environment_manifest_contains_only_six_raw_metadata_archives():
    manifest = json.loads((ROOT / "analysis" / "v3" / "environment_source_archives.json").read_text())
    assert manifest["trait_files_in_scope"] == 0
    assert len(manifest["archives"]) == 6
    assert {row["role"] for row in manifest["archives"]} == {"raw_metadata_chunk"}
    assert 8269246732 not in {row["artifact_id"] for row in manifest["archives"]}
    assert manifest["expected"]["unique_observations"] == 665139


def test_environment_design_cohort_caps_common_taxa_and_excludes_low_support():
    frame = pd.DataFrame(
        _eligible_rows(1, 100)
        + _eligible_rows(2, 10)
        + _eligible_rows(3, 3)
        + _eligible_rows(4, 20, status="introduced")
    )
    selected, report = build_design_cohort(frame, taxon_cap=10, minimum_taxon_support=5)
    counts = selected.groupby("accepted_key").size().to_dict()
    assert counts == {"1": 10, "2": 10}
    assert report["selected_taxa"] == 2
    assert report["selected_rows"] == 20
    assert report["maximum_rows_per_taxon"] == 10
    assert report["largest_selected_taxon_fraction"] == 0.5
    sums = selected.groupby("accepted_key")["equal_taxon_weight"].sum()
    assert ((sums - 1.0).abs() < 1e-12).all()
    assert report["trait_columns_read"] == 0
    assert report["environment_columns_read"] == 0


def test_environment_design_selection_is_deterministic_and_outcome_blind():
    frame = pd.DataFrame(_eligible_rows(10, 30) + _eligible_rows(20, 30))
    first, _ = build_design_cohort(frame, taxon_cap=7, minimum_taxon_support=5)
    shuffled, _ = build_design_cohort(frame.sample(frac=1.0, random_state=99), taxon_cap=7, minimum_taxon_support=5)
    assert first[["accepted_key", "obs_id"]].equals(shuffled[["accepted_key", "obs_id"]])
