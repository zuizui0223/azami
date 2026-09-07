import pandas as pd

from analysis.v3 import build_native_range_join as native


def test_contract_requires_native_primary_and_retains_introduced():
    contract = native.load_contract()
    assert contract["primary_ecology_rule"].startswith("Only observations explicitly classified native")
    assert "transportability" in contract["introduced_rule"]
    assert "not independent image-identification" in contract["misidentification_boundary"]


def test_classify_frame_keeps_all_rows_and_distinguishes_statuses():
    frame = pd.DataFrame(
        {
            "obs_id": ["a", "b", "c"],
            "source_taxon_name": ["Cirsium alpha", "Cirsium alpha", "Cirsium beta"],
            "source_taxon_rank": ["species", "species", "species"],
            "analysis_latitude": [0.5, 0.5, 0.5],
            "analysis_longitude": [0.5, 0.5, 0.5],
            "coordinate_status": ["public_location_present_precision_not_gated"] * 3,
            "captive_state": ["false"] * 3,
            "date_status": ["exact_day"] * 3,
        }
    )
    resolution = pd.DataFrame(
        {
            "input_name": ["Cirsium alpha", "Cirsium beta"],
            "resolution_status": ["resolved_unique_accepted_key", "unresolved"],
            "accepted_key": [1, ""],
            "accepted_name": ["Cirsium alpha", ""],
        }
    )
    distributions = pd.DataFrame(
        {
            "accepted_key": [1],
            "tdwg_level3_code": ["AAA"],
            "classified_status": ["native"],
        }
    )
    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"LEVEL3_COD": "AAA"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
                },
            }
        ],
    }
    joined = native.classify_frame(frame, resolution, distributions, geojson)
    assert len(joined) == 3
    assert joined.loc[joined.obs_id.eq("a"), "native_range_status"].iloc[0] == "native"
    assert joined.loc[joined.obs_id.eq("b"), "native_range_status"].iloc[0] == "native"
    assert joined.loc[joined.obs_id.eq("c"), "native_range_status"].iloc[0] == "unresolved_taxon"


def test_unusable_coordinate_is_not_guessed_native():
    frame = pd.DataFrame(
        {
            "obs_id": ["a"],
            "source_taxon_name": ["Cirsium alpha"],
            "source_taxon_rank": ["species"],
            "analysis_latitude": [None],
            "analysis_longitude": [None],
            "coordinate_status": ["restricted"],
            "captive_state": ["false"],
            "date_status": ["exact_day"],
        }
    )
    resolution = pd.DataFrame(
        {
            "input_name": ["Cirsium alpha"],
            "resolution_status": ["resolved_unique_accepted_key"],
            "accepted_key": [1],
            "accepted_name": ["Cirsium alpha"],
        }
    )
    distributions = pd.DataFrame(
        {"accepted_key": [1], "tdwg_level3_code": ["AAA"], "classified_status": ["native"]}
    )
    geojson = {"type": "FeatureCollection", "features": []}
    joined = native.classify_frame(frame, resolution, distributions, geojson)
    assert joined.loc[0, "native_range_status"] == "unmapped_or_unusable_location"
