from __future__ import annotations

import pandas as pd

from analysis.v3.build_environment_design_cohort import (
    build_design_cohort,
    normalized_state as cohort_normalized_state,
)
from analysis.v3.build_native_range_join_from_source import (
    normalized_state as native_normalized_state,
)


def test_captive_state_normalization_handles_boolean_and_string_roundtrips() -> None:
    values = pd.Series([False, True, "false", "TRUE", " false "])
    expected = ["false", "true", "false", "true", "false"]
    assert cohort_normalized_state(values).tolist() == expected
    assert native_normalized_state(values).tolist() == expected


def test_environment_design_cohort_accepts_boolean_false_from_csv_inference() -> None:
    frame = pd.DataFrame(
        {
            "obs_id": ["1", "2"],
            "accepted_key": [123, 123],
            "accepted_name": ["Cirsium example", "Cirsium example"],
            "source_taxon_name": ["Cirsium example", "Cirsium example"],
            "source_taxon_rank": ["species", "species"],
            "taxon_resolution_status": [
                "resolved_unique_accepted_key",
                "resolved_unique_accepted_key",
            ],
            "native_range_status": ["native", "native"],
            "captive_state": [False, "false"],
            "date_status": ["exact_day", "exact_day"],
            "coordinate_status": [
                "public_location_present_precision_not_gated",
                "public_location_present_precision_not_gated",
            ],
            "analysis_latitude": [35.0, 36.0],
            "analysis_longitude": [135.0, 136.0],
            "observation_month": [6, 7],
        }
    )

    selected, report = build_design_cohort(
        frame,
        taxon_cap=60,
        minimum_taxon_support=1,
    )

    assert len(selected) == 2
    assert report["eligible_rows_before_taxon_support"] == 2
    assert report["selected_taxa"] == 1
    assert report["trait_columns_read"] == 0
    assert report["environment_columns_read"] == 0
