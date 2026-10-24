"""Freeze actual two-species S2 CDF equivalence *non*-claim and source hashes."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"analysis/v4/results/20261022_S2_ECDF_EQUIVALENCE_FRONTIER_VERIFIED.json"
REG=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"

def test_source_authenticated_whole_curve_399_cell_spatial_band_precision():
    r=json.loads(RESULT.read_text())
    assert r["status"].startswith("SOURCE_VERIFIED_FUTURE_EQUIVALENCE_PRECISION_ONLY")
    assert r["source_run_id"]==39711427699
    assert r["artifact"]["id"]==12396863611
    assert r["original_full_curve_source"]["archived_successful_run_id"]==39592586205
    assert r["original_full_curve_source"]["source_member_sha256"]=="765742c0a2fdaef094f26de7999c9c322e5397a3c74ac0912edec4c2d09d9e4f"
    assert r["artifact"]["digest"].startswith("sha256:")
    s=r["result"]["species"]
    assert set(s)=={"Cirsium arvense","Cirsium vulgare"}
    expected={
      "Cirsium arvense":(21870,25657,.046184870442695036,.12346166441399255),
      "Cirsium vulgare":(16246,30287,.044243156739647405,.11459217572205463),
    }
    for taxon,(native,introduced,observed,limit) in expected.items():
        front=s[taxon]["photo_angle_frontier"]
        assert (front["n_native_original_cell_balanced_observations"],
                front["n_introduced_original_cell_balanced_observations"])==(native,introduced)
        assert front["ECDF_bootstrap_draws"]==399
        assert front["display_knots"]==401
        assert front["all_95pct_CDF_difference_bands_include_zero_at_reported_knots"] is True
        assert abs(front["observed_sup_absolute_CDF_difference"]-observed)<1e-12
        assert abs(front["minimum_delta_for_95pct_simultaneous_band_containment"]-limit)<1e-12
        assert front["NOT_a_calibrated_formal_equivalence_TOST"] is True
        assert front["NOT_a_prospectively_chosen_SESOI"] is True
        assert front["NOT_a_proof_of_unchanged_biological_head_posture"] is True
        assert s[taxon]["rain_source_distribution"]["all_95pct_CDF_difference_bands_include_zero_at_reported_knots"] is False


def test_no_posthoc_delta_cherrypicking_and_positive_wet_niche_preserved():
    data=json.loads(RESULT.read_text())
    report=data["sensitivity_without_SESOI"]
    assert report["5pct_difference_halfwidth_supported_both"] is False
    assert report["10pct_difference_halfwidth_supported_both"] is False
    assert report["15pct_difference_halfwidth_band_containment_both"] is True
    assert report["15pct_halfwidth_independently_ecologically_validated"] is False
    assert data["nonclaim"]["formally_conducted_predeclared_equivalence_TOST"] is False
    assert data["nonclaim"]["concluded_no_biologically_important_change"] is False
    assert data["nonclaim"]["S2_positive_introduced_wetter_climate_shift_retained"] is True
    assert data["nonclaim"]["S1_between_taxon_positive_BIO12_effect_retained"] is True
    x=json.loads(REG.read_text())["S2"]["source_verified_photo_angle_equivalence_feasibility_20261022"]
    assert x["true_formal_equivalence_TOST_completed"] is False
    assert x["original_S2_wetter_introduced_occurrence_niche_positive_kept"] is True
