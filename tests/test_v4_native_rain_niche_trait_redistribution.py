"""S2 scientific tests: rainfall opportunities and phenotype responses stay distinct."""
import json
import numpy as np
import pandas as pd
from analysis.v4.native_rain_niche_trait_redistribution import (
    RAIN,BANDS,_histogram,_composition,_bootstrap,
    _rain_classes,_background_profile,_area_and_source,_opportunity_shift
)
from analysis.v4.native_climate_trait_coupling import weighted_quantiles

def synthetic(seed=312):
    rng=np.random.default_rng(seed)
    parts={}
    # Repeated individuals within cells deliberately imbalance the samples:
    # geographic cell weighting should be respected by every stage.
    for status in ("native","introduced"):
        rows=[]
        for cell in range(36):
            for j in range(32):
                wet=cell%6==0
                rain=(
                    430+17*cell+5*j if status=="native" else
                    470+18*cell+5*j+75*wet)
                head=(67+0.015*(rain-700)
                      + 7*(status=="introduced")
                      +rng.normal(0,9))
                rows.append({
                    "_cell":cell+ (0 if status=="native" else 1000),
                    "obs_id":f"{status}-{cell}-{j}",
                    "user_id":f"{status}-u{cell}-{j}",
                    RAIN:float(rain),
                    "orientation_angle_degrees_median":float(np.clip(head,0,180)),
                })
        parts[status]=pd.DataFrame(rows)
    native=parts["native"]
    weights=1/native.groupby("_cell").obs_id.transform("size").to_numpy(float)
    cut=weighted_quantiles(native[RAIN].to_numpy(float),weights,
                           quantiles=(.25,.50,.75,.90,.95))
    return parts["native"],parts["introduced"],cut

def test_source_rainfall_cutoffs_ignore_phenotype():
    a,b,cut=synthetic()
    assert len(cut)==5 and sorted(cut)==cut
    labels=_rain_classes(b,cut)
    assert labels.max()==5 and labels.min()>=0
    inverted=a.copy()
    inverted["orientation_angle_degrees_median"]=180-inverted[
        "orientation_angle_degrees_median"]
    assert np.array_equal(_rain_classes(a,cut),_rain_classes(inverted,cut))

def test_area_opportunity_is_separate_from_occurrence_density():
    a,b,cut=synthetic()
    occ=_histogram(b,cut)
    fakebg=b[[RAIN]].copy()
    fakebg["_area_weight"]=np.ones(len(fakebg))
    bg=_background_profile(fakebg,cut)
    result=_area_and_source(occ,bg)
    assert result["status"]=="DESCRIPTIVE_SOURCE_VS_AVAILABLE_CLIMATE"
    assert abs(sum(result["source_occurrence_probabilities"])-1)<1e-10
    assert abs(sum(result["available_area_probabilities"])-1)<1e-10

def test_opportunity_shift_decomposes_wetter_occurrence_vs_background():
    a,b,cut=synthetic()
    native=_histogram(a,cut)
    introduced=_histogram(b,cut)
    nbg=a[[RAIN]].copy()
    ibg=b[[RAIN]].copy()
    nbg["_area_weight"]=1.
    ibg["_area_weight"]=1.
    out=_opportunity_shift(native,introduced,
                           _background_profile(nbg,cut),
                           _background_profile(ibg,cut))
    assert out["status"]=="DECOMPOSED_DESCRIPTIVE_WET_OPPORTUNITY_VS_OCCUPANCY"
    assert not out["not_establishment_or_selection_causal_effect"] is False
    for band in out["thresholds"].values():
        assert np.isclose(
            band["recorded_occurrence_shift_introduced_minus_native"]-
            band["available_climate_area_shift_introduced_minus_native"],
            band["occurrence_shift_minus_available_area_shift"]
        )


def test_rain_mix_decomposition_identity_and_nonzero_condition():
    a,b,cut=synthetic()
    sh1=_histogram(a,cut)
    sh2=_histogram(b,cut)
    ph1=_histogram(a,cut,angle=True)
    ph2=_histogram(b,cut,angle=True)
    r=_composition(sh1,sh2,ph1,ph2)
    assert r["status"]=="SUPPORTED_DESCRIPTIVE_RAIN_MIXTURE_DECOMPOSITION"
    assert len(r["supported_bands"])>=2
    assert abs(
        r["image_down_fraction_observed_introduced_minus_native"]-
        r["image_down_fraction_rain_composition_part"]-
        r["image_down_fraction_within_rain_band_part"])<1e-10
    assert r["angle_distribution_W1_native_vs_introduced_degrees"]>=0
    assert len(r["native_observed"]["Q10_Q50_Q90_image_degrees"])==3
    # Never export raw photographic identities in scientific receipt.
    raw=json.dumps(r)
    for field in ('"obs_id"','"user_id"','"latitude"','"longitude"'):
        assert field not in raw

def test_common_rain_bands_stop_if_no_photo_support():
    a,b,cut=synthetic()
    ph1=_histogram(a.head(8),cut,angle=True)
    ph2=_histogram(b.head(8),cut,angle=True)
    r=_composition(_histogram(a,cut),_histogram(b,cut),ph1,ph2)
    assert r["status"]=="STOP_NO_RAIN_BAND_SHARED_PHOTOGRAPHIC_SUPPORT"
    assert r["native_source_mass_supported"]<.65

def test_spatial_bootstrap_keeps_photo_source_scope_and_noncausal():
    a,b,cut=synthetic()
    sh1=_histogram(a,cut)
    sh2=_histogram(b,cut)
    ph1=_histogram(a,cut,angle=True)
    ph2=_histogram(b,cut,angle=True)
    r=_bootstrap(sh1,sh2,ph1,ph2,label="synthetic/unit")
    assert r["status"] in (
        "EXPLORATORY_FIXED_BAND_REGION_CELL_RESAMPLING",
        "NOT_ESTIMABLE_TOO_FEW_VALID_GEO_BOOTSTRAP_DRAWS")
    if r["status"].startswith("EXPLORATORY"):
        assert r["native_rain_bin_edges_fixed_before_resampling"]
        assert r["no_camera_genetic_or_climate_confounder_inference"]
        for interval in r["image_down_probability_delta_intervals_2p5_50_97p5"].values():
            assert interval[0]<=interval[1]<=interval[2]
