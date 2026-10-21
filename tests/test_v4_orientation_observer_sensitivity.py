"""Unit and invariant checks for S1 original-observer sensitivity.

All toy source identifiers are synthetic and never written to user-facing
artifacts. These do not assert an actual human-camera variance contribution.
"""
import numpy as np
import pandas as pd
import pytest

from analysis.v4.orientation_observer_sensitivity import (
    ENDPOINTS, FIXED_PREDICTORS, unbalanced_observer_component,
    residual_from_fixed_species_fit,effect_line,
)


def _toy_repeated_residuals(n_cells=30,photographers=4,repeats=5):
    rows=[]
    rng=np.random.default_rng(20261021)
    for cell in range(n_cells):
        for user in range(photographers):
            camera=rng.normal(0,3)
            for rep in range(repeats):
                rows.append({
                    "accepted_key":str(100+cell%6),
                    "_cell":cell,
                    "user_id":f"synthetic_{cell}_{user}",
                    "_residual":camera+rng.normal(0,.2),
                })
    return pd.DataFrame(rows)


def test_true_repeated_observer_residuals_have_nested_component():
    f=_toy_repeated_residuals()
    result=unbalanced_observer_component(f,replicates=40)
    assert result["status"].startswith("estimated_nested_observer")
    assert result["taxon_cell_strata_with_repeat_and_two_authors"]==30
    assert result["observation_rows_in_valid_cells"]==30*4*5
    assert result["within_observer_extra_residual_df"]==30*4*(5-1)
    assert result["unbalanced_MoM"]["observer_nested_variance"]>0
    assert result["unbalanced_MoM"]["observer_nested_fraction"]>.5
    assert len(result["source_cell_bootstrap_fraction_quantiles"])==3
    assert result["source_cell_bootstrap_fraction_quantiles"][0]>=0
    assert "not pure camera" in result["scope"] if "scope" in result else True


def test_one_obs_per_user_is_not_a_camera_measurement_replicate():
    f=_toy_repeated_residuals(repeats=1)
    result=unbalanced_observer_component(f,replicates=40)
    assert result["status"]=="insufficient_multiauthor_source_cells"
    assert not result["primary_reportable_random_camera_variance"]


def test_nontaxon_camera_effect_cannot_be_emitted_without_original_observer():
    f=_toy_repeated_residuals()
    f.loc[0,"user_id"]=None
    with pytest.raises(ValueError,match="Observer variance input"):
        unbalanced_observer_component(f,replicates=40)


def test_reference_rain_between_within_scale_is_declared_100mm_not_1_sd():
    scales={"chelsa_bio12":{"mean":700.,"sd":200.}}
    c={
        "status":"estimated","n":170066,"n_species":244,
        "variance_ratio":.3,
        "effects":{
            "chelsa_bio12__between":{
                "raw_units_per_environment_sd":.8763983845770414,
                "raw_ci95":[.4116126478448966,1.3411841213091862]
            },
            "chelsa_bio12__within":{
                "raw_units_per_environment_sd":-.08859642159991839,
                "raw_ci95":[-.43769796374799815,.26050512054816128]
            }
        },
    }
    result=effect_line(c,scales)
    assert result["between_rain_deg_per_100mm"]==pytest.approx(.4381991922885207)
    assert result["within_rain_deg_per_100mm"]==pytest.approx(-.044298210799959194)
    assert result["between_species_pointwise_95_interval_deg_per_100mm"]==pytest.approx(
        [.2058063239224483,.6705920606545931]
    )
    assert result["observations"]==170066
    assert result["accepted_taxa"]==244


def test_residual_uses_prespecified_native_between_means_not_subset_recentre():
    frame=pd.DataFrame({
        "accepted_key":["a","a"],"orientation_image_vertical_angle":[42.,47.],
    })
    scales={}
    effect={}
    for variable in FIXED_PREDICTORS:
        frame[variable+"__between"]=[10.,10.]
        frame[variable+"__within"]=[2.,3.]
        scales[variable]={"mean":10.,"sd":2.}
        effect[variable+"__between"]={"raw_units_per_environment_sd":0.}
        effect[variable+"__within"]={"raw_units_per_environment_sd":0.}
    effect["chelsa_bio12__within"]["raw_units_per_environment_sd"]=2.
    selection={"selected":FIXED_PREDICTORS,"scales":scales}
    component={"status":"estimated","reference_prediction":40.,"effects":effect}
    result=residual_from_fixed_species_fit(
        frame,component,"orientation_image_vertical_angle",selection
    )
    assert result.tolist()==pytest.approx([0.,4.])
    dropped=frame.drop(columns=["chelsa_bio12__between"])
    with pytest.raises(ValueError,match="full-native species means"):
        residual_from_fixed_species_fit(
            dropped,component,"orientation_image_vertical_angle",selection)
