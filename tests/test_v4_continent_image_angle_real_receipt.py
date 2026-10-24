"""Real-action evidence constraints for introduced Cirsium vulgare photo angles.

These tests freeze aggregate evidence; do not expose original observations,
human photographers, original image pixels or exact coordinates.
"""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"analysis/v4/results/20261021_unified_abcde"

def receipt():
    return json.loads((D/"S2_VULGARE_ACTUAL_CONTINENT_ORIENTATION_20261021.json").read_text(encoding="utf-8"))


def test_original_continental_source_measures_are_real_data_not_simulated():
    r=receipt()
    assert r["source"]["run"]==39550930196
    assert r["source"]["artifact_id"]==12323896390
    assert r["source"]["member_sha256"]==(
        "3ee0390d5c64de53d4b28e994007278e1fe9a5837eef20cf1ebe00761260f94b"
    )
    assert r["source"]["focal_original_orientation_eligible_observations"]==173737
    assert r["source"]["focal_after_date_original_observer_cell_balancing"]==94060
    assert not r["source"]["true_photo_individual_genet_identity_known"]


def test_continents_do_not_show_supported_uniform_shift_toward_photo_nodding():
    r=receipt()
    n=r["source_native"]["Cirsium vulgare"]
    assert n["n_balanced"]==16246 and n["cells"]==124
    assert 66<=n["cell_mean_angle"]<=67
    provinces=r["Cirsium vulgare"]
    names=("North America","Africa","South America","Australasia")
    for name in names:
        v=provinces[name]
        assert v["original_orientation_eligible"]>2000
        assert v["n_balanced"]>1000
        assert v["cells"]>=15
        lo,med,hi=v["cell_bootstrap_q025_q50_q975"]
        assert lo<=v["introduced_minus_native_image_angle"]<=hi
        assert lo<0<hi
        assert abs(v["introduced_minus_native_image_angle"])<2
    warm=provinces["intro_warm_North_America"]
    assert warm["above_native_BIO1_q95"]["original_balanced_photo_observations"]==507
    assert warm["above_native_BIO1_q95"]["cells"]==16
    assert warm["at_or_below_native_BIO1_q95"]["original_balanced_photo_observations"]==24839
    assert warm["at_or_below_native_BIO1_q95"]["cells"]==75
    assert warm["above_q95_minus_other_image_angle"]>0
    assert warm["cell_bootstrap_q025_q50_q975"][0]<0
    assert warm["cell_bootstrap_q025_q50_q975"][-1]>0
    assert warm["interval_includes_zero"]


def test_climate_missingness_blocks_Africa_and_South_America_bio15():
    r=receipt()
    d=r["Cirsium vulgare"]
    assert d["intro_precip_seasonality_Africa"]["five_original_pixel_climate_available_unbalanced"]==72
    assert d["intro_precip_seasonality_South_America"]["five_original_pixel_climate_available_unbalanced"]==0
    assert not d["intro_precip_seasonality_Africa"]["qualified_BIO15_novel_x_head_angle"]
    assert not d["intro_precip_seasonality_South_America"]["qualified_BIO15_novel_x_head_angle"]
    assert d["intro_precip_seasonality_Africa"]["verdict"].startswith("STOP_")
    assert d["intro_precip_seasonality_South_America"]["verdict"].startswith("STOP_")
    assert r["source"]["source_reference_sha256"]==(
        "bbfc9a459a1a5c230c31a405671633e3032491f92f0123634e5cbb8f4db33562"
    )


def test_same_image_orientation_not_a_gravity_nodding_certification():
    r=receipt()
    assert "image-vertical" in r["response_unit"]
    assert "no botanical gravity calibration" in r["response_unit"]
    assert any("gravitational" in x.lower() for x in r["limitations"])
    assert any("continent" in x.lower() for x in r["limitations"])
    j=(D/"S2_VULGARE_CONTINENT_IMAGE_ANGLE_RESULT_JA.md").read_text(encoding="utf-8")
    assert "気温" in j or "温暖" in j
    assert "降水季節性" in j
    assert "未判定" in j or "判断できない" in j
