"""Source-verified study boundaries: all 350 introduced screening and climatic novelty.

The source archive holds only aggregate data; no row-level observations,
precise coordinates, photographer IDs or images are re-exported by this test.
"""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
E=ROOT/"analysis/v4/results/20261021_unified_abcde"


def read(name):
    return json.loads((E/name).read_text(encoding="utf-8"))


def test_niche_metric_is_real_source_run_not_a_second_species_selection():
    g=read("S2_CLIMATE_NICHE_NOVELTY_20261021.json")
    a=read("S2_ALL_350_SPECIES_SOURCE_SCREEN_VERIFIED_20261021.json")
    assert (g["source"]["action_run"],g["source"]["artifact_id"])==(
        39545903306,12320759155)
    assert g["source"]["source_member_sha256"]==(
        "8f578b6548e3ce58cbe6b3c23147544353d6aeff46b7337c20120065c01310c4")
    assert g["source"]["all_350_species_screen_in_same_run"]
    assert len(g["results"])==a["evidence_census"]["both_range_stronger_gate_species"]==2
    assert set(g["results"])=={
        x["source_name"] for x in a["all_both_range_species"] if x["stronger_pass"]
    }
    assert not g["contract"]["climate_novelty_is_realized_fundamental_niche_shift"]
    assert not g["contract"]["statistical_significance_or_phylogenetic_inference_from_these_fractions"]


def test_equal_site_cell_weighting_changes_niche_fraction_and_has_native_baseline():
    g=read("S2_CLIMATE_NICHE_NOVELTY_20261021.json")
    for name,z in g["results"].items():
        photo=z["photo_observation_weight"]
        geo=z["cell_equal_weight"]
        assert 0 < geo["introduced_outside_at_least_one_native_marginal_fraction"] < 1
        assert 0 < photo["introduced_outside_at_least_one_native_marginal_fraction"] < 1
        assert photo["introduced_outside_at_least_one_native_marginal_fraction"] > (
            geo["introduced_outside_at_least_one_native_marginal_fraction"])
        assert abs(
            geo["introduced_outside_at_least_one_native_marginal_fraction"] -
            (1-geo["native_own_fiveD_box_inside_fraction"]) -
            geo["introduced_excess_outside_over_native_internal_reference"]
        ) < 1e-12
        assert set(geo["introduced_outside_native_5_to_95_per_axis"])=={
            "bio01","bio04","bio12","bio15","rsds_mean"
        }
        assert z["naive_photo_minus_equal_cell_novelty_fraction"] > .10
    assert g["results"]["Cirsium arvense"]["cell_equal_weight"][
        "introduced_excess_outside_over_native_internal_reference"] < .05
    assert g["results"]["Cirsium vulgare"]["cell_equal_weight"][
        "introduced_excess_outside_over_native_internal_reference"] < .12


def test_v4_current_rain_question_does_not_promote_photo_niche_to_adaptation():
    registry=json.loads((ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json")
                        .read_text(encoding="utf-8"))
    data=registry["S2"]["climate_transfer"]["climate_niche_occupancy"]
    assert data["github_actions_run"]==39545903306
    assert data["novel_climate_adaptive_evolution"] is False
    assert data["per_original_geo_cell_weighting_changes_magnitude"]
    assert data["continent_global_invasion_conclusion_supported"] is False
    assert registry["T2"]["existing_52_trees_sufficient_for_dated_origin_counts"] is False
