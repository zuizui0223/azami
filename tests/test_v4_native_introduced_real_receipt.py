"""Cross-check source-verified S2 common-support science/claim boundary.

No private original photo/observation/user identifiers or raw climate matrices.
"""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"analysis/v4/results/20261021_unified_abcde"


def _open(name):
    return json.loads((DIR/name).read_text(encoding="utf-8"))


def test_S2_same_source_after_climate_support_gate_has_hashed_provenance():
    d=_open("S2_FIVE_EXPOSURE_COMMON_SUPPORT_VERIFIED_20261021.json")
    p=_open("S2_REAL_FIRST_PASS_SLOPES_20261021.json")
    assert d["status"].startswith("EXECUTED_SOURCE_VERIFIED")
    assert (d["source"]["github_run"],d["source"]["artifact_id"])==(
        39543173995,12320715523
    )
    assert d["source"]["common_support_member_sha256"]==(
        "1ef98d6ecd29cee86f02a29bf4b1dddc10deca70ce7d131d416d051d03edc115"
    )
    assert d["source"]["source_reference_transfer_sha256"]==(
        "a95e0eaa085cc4cbd55a8302f5d527e26e3800761f7bb043ae9e65bfecdbdb26"
    )
    assert d["source"]["original_source_status_sha256"]==p["source"]["source_status_member_sha256"]
    assert d["design"]["native_primary_observations"]==313350
    assert d["design"]["introduced_explicit_observations"]==156750
    assert d["design"]["five_D_density_overlap_confirmed"] is False
    assert d["claim_boundary"]["not_equivalent_to_zero"]
    assert d["claim_boundary"]["in_range_joint_five_climates_BIO12_contrast_crosses_zero"]


def test_selected_native_common_support_collapses_but_interval_not_positive():
    d=_open("S2_FIVE_EXPOSURE_COMMON_SUPPORT_VERIFIED_20261021.json")
    assert set(d["results"])=={"Cirsium arvense","Cirsium vulgare"}
    assert d["design"]["exposures"]==[
        "chelsa_bio01","chelsa_bio04","chelsa_bio12",
        "chelsa_bio15","chelsa_rsds_mean"
    ]
    for taxon,n_expected,i_expected in (
        ("Cirsium arvense",600,6428),
        ("Cirsium vulgare",1198,10011),
    ):
        x=d["results"][taxon]
        a=x["five_exposure_rectangle_overlap"]
        assert a["native_rows"]==n_expected
        assert a["introduced_rows"]==i_expected
        assert a["native_cells"]==23
        assert a["introduced_cells"]>=27
        assert a["native_fraction_of_observer_cell_cohort"]<.08
        assert a["difference_ci95"][0]<0<a["difference_ci95"][1]
        assert x["pre_overlap"]["difference_ci95"][0]<0<x["pre_overlap"]["difference_ci95"][1]
        assert x["rain_only_overlap"]["difference_ci95"][0]<0<x["rain_only_overlap"]["difference_ci95"][1]
        assert x["introduced_continents_evaluable"]==["North America"] if "introduced_continents_evaluable" in x else x["introduced_continent_evaluable"]==["North America"]
    assert d["results"]["Cirsium arvense"]["pre_overlap"]["native_to_intro_skill_vs_native_constant"]<0
    assert d["results"]["Cirsium arvense"]["five_exposure_rectangle_overlap"]["native_to_intro_skill_vs_native_constant"]>0
    assert d["results"]["Cirsium vulgare"]["five_exposure_rectangle_overlap"]["native_to_intro_skill_vs_native_constant"]<0


def test_current_S1_S2_T2_T3_plan_tracks_no_causal_introduction_inference():
    decision=json.loads((
        ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"
    ).read_text(encoding="utf-8"))
    s=decision["S2"]["climate_transfer"]["common_support_robustness"]
    assert s["run"]==39543173995
    assert s["arvense"]["native_support_observations"]==600
    assert s["vulgare"]["native_support_observations"]==1198
    assert s["all_slope_difference_intervals_cross_zero"]
    assert s["strong_focal_evolutionary_fixed_trait_evidence"] is False
    assert decision["quality_controls"]["introduced_native_climate_common_support_transfer_validated"] is False
    assert decision["quality_controls"]["introduced_native_five_exposure_overlap_falsification_executed"] is True
    assert decision["T2"]["existing_52_trees_sufficient_for_dated_origin_counts"] is False
    assert decision["T3"]["cross_tribe_convergence_inferred"] is False
