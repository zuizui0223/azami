import json
from pathlib import Path

import pandas as pd

from analysis.v3.build_ecological_source_cohort import build


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "ecological_source_cohort_contract.json"


def load_contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_contract_keeps_full_species_level_native_source_cohort_without_caps():
    c = load_contract()
    assert c["status"] == "fixed_before_original_trait_execution_and_trait_environment_join"
    assert c["expected_phenotype_blind_denominator_from_completed_upstream"] == {"rows": 319244, "accepted_taxa": 355}
    assert c["sampling_rules"]["row_cap"] is None
    assert c["sampling_rules"]["taxon_cap"] is None
    assert c["sampling_rules"]["spatial_thinning"] is False
    assert c["sampling_rules"]["top_taxon_deletion"] is False
    assert c["separation"]["trait_values_read"] == 0
    assert c["separation"]["environment_values_used_for_cohort_selection"] == 0


def test_postmeasurement_taxon_reporting_support_is_fixed_before_coefficients():
    rules = load_contract()["post_measurement_support_rules_fixed_before_coefficients"]
    assert rules["taxon_specific_slope_reporting"]["minimum_distinct_observations_with_qualified_endpoint"] == 10
    assert rules["taxon_specific_slope_reporting"]["minimum_distinct_quarter_degree_coordinate_cells"] == 4
    assert rules["among_taxon_summary_reporting"]["minimum_distinct_observations_with_qualified_endpoint"] == 5


def test_builder_filters_only_declared_source_states(monkeypatch):
    c = load_contract()
    c["expected_phenotype_blind_denominator_from_completed_upstream"] = {"rows": 2, "accepted_taxa": 1}
    base = {
        "accepted_key": 11,
        "accepted_name": "Cirsium testum",
        "source_taxon_name": "Cirsium testum",
        "source_taxon_rank": "species",
        "taxon_resolution_status": "resolved_unique_accepted_key",
        "native_range_status": "native",
        "captive_state": False,
        "date_status": "exact_day",
        "coordinate_status": "public_location_present_precision_not_gated",
        "analysis_latitude": 35.1,
        "analysis_longitude": 135.1,
        "observation_month": 6,
        "position_accuracy_status": "reported",
        "position_accuracy_m": 30,
    }
    rows = []
    for obs, lat in [("1", 35.1), ("2", 35.4)]:
        r = dict(base, obs_id=obs, analysis_latitude=lat)
        rows.append(r)
    rows.append(dict(base, obs_id="3", native_range_status="introduced"))
    rows.append(dict(base, obs_id="4", source_taxon_rank="genus"))
    cohort, report = build(pd.DataFrame(rows), c)
    assert list(cohort["obs_id"]) == ["1", "2"]
    assert report["rows"] == 2
    assert report["accepted_taxa"] == 1
    assert report["row_cap"] is None and report["spatial_thinning"] is False
    assert report["trait_columns_read"] == 0 and report["environment_columns_read"] == 0
