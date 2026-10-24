"""Freeze real source-wide (all 350 source species) WCVP introduction screen.

The known seven native+introduced taxa are an exhaustively screened current
source subset. The long two-species climate experiment cannot claim genus-wide
plasticity or invasion biology beyond source/measurability resolution.
"""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"analysis/v4/results/20261021_unified_abcde/"
def audit():
    return json.loads((R/"S2_ALL_350_SPECIES_SOURCE_SCREEN_VERIFIED_20261021.json").read_text())


def test_genuine_source_artifact_link_and_350_species_universe():
    j=audit()
    assert j["source"]["github_actions_run"]==39544849439
    assert j["source"]["artifact_id"]==12321361815
    assert j["source"]["source_member_sha256"]==(
        "9d7019098874664611fd048a6f920432b1645957bca9af267a13f3401c92cf5b")
    assert j["source"]["original_source_observations"]==665115
    assert j["source"]["pre_detector_source_names"]==350
    assert j["evidence_census"]["source_species_rank_names"]==350
    assert j["evidence_census"]["source_species_names_with_any_accepted_wild_record"]==313
    assert j["evidence_census"]["source_species_rank_names"]==sum(
        j["evidence_census"][k] for k in (
            "wild_native_only_species_names","wild_introduced_only_species_names",
            "wild_native_and_explicit_introduced_species_names",
            "neither_wild_source_range_species_names"
        )
    )
    assert j["evidence_census"]["wild_native_and_explicit_introduced_species_names"]==7
    assert j["evidence_census"]["both_range_basic_gate_species"]==2
    assert j["evidence_census"]["both_range_stronger_gate_species"]==2
    share=j["evidence_census"]["two_qualified_source_names_dominance"]
    assert share["introduced_wild_observations_from_two_species"]==156608
    assert share["introduced_wild_observations_in_five_sparse_species"]==142
    assert share["introduced_orientation_eligible_from_two_species"]==91028
    assert share["introduced_orientation_eligible_from_five_sparse_species"]==78
    assert share["proportion_of_all_WCVP_introduced_source_wild_from_two_species"]>0.999
    assert share["proportion_of_all_WCVP_introduced_orientation_from_two_species"]>0.999


def test_literature_source_names_remaining_are_recorded_as_non_evaluable_not_null():
    j=audit()
    c=j["all_both_range_species"]
    assert len(c)==7
    assert len({z["WCVP_key"] for z in c})==7
    passnames={z["source_name"] for z in c if z["basic_pass"]}
    assert passnames=={"Cirsium arvense","Cirsium vulgare"}
    assert {z["source_name"] for z in c if z["stronger_pass"]}==passnames
    for z in c:
        min_n=j["predeclared_gates"]["basic"]["orientation_eligible_per_range"]
        min_cells=j["predeclared_gates"]["basic"]["5_degree_equal_sin_cells_per_range"]
        min_obs=j["predeclared_gates"]["basic"]["original_observers_per_range"]
        expected=all(
            z[f"{status}_orientation_eligible"]>=min_n
            and z[f"{status}_cells"]>=min_cells
            and z[f"{status}_observers"]>=min_obs
            for status in ("native","introduced")
        )
        assert expected==z["basic_pass"]
    assert {z["source_name"] for z in c if not z["basic_pass"]}=={
        "Cirsium heterophyllum","Cirsium mexicanum","Cirsium oleraceum",
        "Cirsium undulatum","Cirsium esculentum"
    }


def test_all_350_source_results_not_claimed_as_350_species_climate_tests():
    j=audit()
    assert j["source"]["all_source_taxa_checked_without_rain_or_orientation_values"]
    assert j["selection_and_interpretation"]["two_species_do_not_establish_cross_species_distribution"]
    assert j["selection_and_interpretation"]["no_photo_raw_user_coordinates_published"]
    note=(R/"S2_ALL_SPECIES_RANGE_SUPPORT_RESULT_JA.md").read_text()
    assert "350" in note and "7" in note
    assert "未実行" in note or "未検証" in note or "データ不足" in note
