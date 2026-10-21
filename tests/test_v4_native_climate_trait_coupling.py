"""Test phenotype co-variation WITH climatic novelty on the same source records.

All observations, observers, coordinates and head angles in these tests are
synthetic. Only GitHub Actions with sealed source can supply real evidence.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_climate_trait_coupling import (
    CLASSES, REF, climate_only_groups, analyse_niche_trait_pair
)
from analysis.v4.native_introduced_climate_transfer import standardization


def toy_intro(n=480,seed=51):
    rng=np.random.default_rng(seed)
    frames={}
    for state in ("native","introduced"):
        rows=[]
        for i in range(n):
            k=i%48
            v={x:float(rng.normal()) for x in VARIABLES}
            rain=800+100*rng.normal()
            if state=="introduced":
                if i%4==0:
                    rain+=600
                elif i%4==1:
                    rain-=600
                elif i%4==2:
                    v["chelsa_bio01"]+=8
            v["chelsa_bio12"]=float(rain)
            y=float(80+0.017*(rain-800)+1.1*v["chelsa_bio01"]+rng.normal(scale=2))
            row={
                "obs_id":f"{state}_{i:05d}",
                "user_id":f"author_{state}_{i%126}",
                "range_state":state,
                "taxon_name":"Cirsium test species",
                "_cell":k,
                "_month":i%12+1,
                "_year":2019+i%5,
                "orientation_angle_degrees_median":np.clip(y,0,180),
            }
            row.update(v)
            rows.append(row)
        frames[state]=pd.DataFrame(rows)
    return frames["native"],frames["introduced"]


def test_climatic_bins_are_rain_wet_rain_dry_other_novel_or_analogue_before_trait():
    n,i=toy_intro()
    labels,method=climate_only_groups(n,i)
    assert set(labels)==set(CLASSES)
    assert (labels=="rain_wetter_than_native").sum()>100
    assert (labels=="rain_drier_than_native").sum()>100
    assert (labels=="other_climate_novel").sum()>100
    assert (labels==REF).sum()>50
    assert method["climate_group_computed_without_orientation_values"]
    mutated=i.copy()
    mutated["orientation_angle_degrees_median"]=180-i.orientation_angle_degrees_median
    same,contract=climate_only_groups(n,mutated)
    assert list(labels)==list(same)
    assert contract==method


def test_observed_angle_jointly_tracks_a_constructed_novel_wetter_climate_without_claiming_evolution():
    n,i=toy_intro()
    scales=standardization(pd.concat([n,i],ignore_index=True))
    result=analyse_niche_trait_pair(n,i,scales,bootstrap=300)
    assert result["status"].startswith("EXECUTED_CLIMATE_NOVELTY_X_PHOTO_ANGLE")
    group=result["introduced_four_climate_categories"]
    assert sum(z["photo_observations"] for z in group.values())==len(i)
    assert set(group)==set(CLASSES)
    assert all(x["satisfies_predeclared_descriptive_support"] for x in group.values())
    assert group["rain_wetter_than_native"]["observed_mean_image_angle_degrees"] > (
        group[REF]["observed_mean_image_angle_degrees"])
    assert group["rain_drier_than_native"]["observed_mean_image_angle_degrees"] < (
        group[REF]["observed_mean_image_angle_degrees"])
    contrasts=result["categories_vs_climate_analogue"]
    for label in ("rain_wetter_than_native","rain_drier_than_native","other_climate_novel"):
        rec=contrasts[label]
        assert rec["status"].startswith("INTRODUCED_CLIMATE_CONTRAST_ESTIMATED")
        assert rec["bootstrap_replicates_valid"]>250
        assert rec["not_gravity_or_functional_protection_evidence"]
        assert rec["observed_image_angle_difference_degrees"]==pytest.approx(
            rec["native_model_expected_image_angle_difference_degrees"]+
            rec["observed_minus_expected_residual_difference_degrees"],abs=1e-9)
    assert contrasts["rain_wetter_than_native"][
        "image_orientation_direction_consistent_with_rain_prediction"] is True
    assert contrasts["rain_drier_than_native"][
        "image_orientation_direction_consistent_with_rain_prediction"] is True
    assert contrasts["other_climate_novel"][
        "image_orientation_direction_consistent_with_rain_prediction"] is None
    assert any("not longitudinal" in x or "not genetic" in x or "not evolution" in x
               for x in result["technical_and_causal_boundaries"])


def test_missing_rain_novelty_does_not_become_zero_trait_effect():
    n,i=toy_intro()
    i=i.copy()
    for p in VARIABLES:
        i[p]=n[p].sample(n=len(i),replace=True,random_state=10).to_numpy()
    i["chelsa_bio12"]=850
    i["_cell"]=np.arange(len(i))%48
    scales=standardization(pd.concat([n,i],ignore_index=True))
    x=analyse_niche_trait_pair(n,i,scales,bootstrap=200)
    assert x["introduced_four_climate_categories"]["rain_wetter_than_native"][
        "photo_observations"]==0
    assert x["categories_vs_climate_analogue"]["rain_wetter_than_native"][
        "status"]=="STOP_INSUFFICIENT_CATEGORY_OR_REFERENCE_SUPPORT"
    assert "observed_image_angle_difference_degrees" not in x["categories_vs_climate_analogue"][
        "rain_wetter_than_native"]


def test_report_refuses_source_nulls_and_nonindependent_ids():
    n,i=toy_intro(n=120)
    sd=standardization(pd.concat([n,i],ignore_index=True))
    bad=i.copy()
    bad.loc[0,"chelsa_bio12"]=np.nan
    with pytest.raises(ValueError,match="Missing original CHELSA"):
        climate_only_groups(n,bad)
    dup=i.copy()
    dup.loc[0,"obs_id"]=dup.loc[1,"obs_id"]
    with pytest.raises(ValueError,match="Non-independent"):
        analyse_niche_trait_pair(n,dup,sd,bootstrap=200)
    with pytest.raises(ValueError,match="original"):
        analyse_niche_trait_pair(i,n,sd,bootstrap=200)
