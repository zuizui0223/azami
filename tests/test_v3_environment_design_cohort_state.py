import pandas as pd

from analysis.v3.build_environment_design_cohort import build_design_cohort
from analysis.v3.build_native_range_join_from_source import normalized_state


def test_false_captive_state_survives_bool_and_string_csv_forms():
    values = pd.Series([False, "false", " FALSE "])
    assert normalized_state(values).tolist() == ["false", "false", "false"]

    frame = pd.DataFrame(
        {
            "obs_id": ["1", "2"],
            "accepted_key": [12345, 12345],
            "accepted_name": ["Cirsium test", "Cirsium test"],
            "source_taxon_name": ["Cirsium test", "Cirsium test"],
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

    assert set(selected["obs_id"].astype(str)) == {"1", "2"}
    assert report["eligible_rows_before_taxon_support"] == 2
    assert report["selected_rows"] == 2
    assert report["selected_taxa"] == 1
    assert report["trait_columns_read"] == 0
    assert report["environment_columns_read"] == 0
