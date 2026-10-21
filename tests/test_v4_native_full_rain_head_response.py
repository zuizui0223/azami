"""Synthetic-only tests for annual rain -> image-head angle in newly occupied climates.

Annual BIO12 is the primary S1/S2 test. BIO15 is a secondary, previously
uncertain S1 axis. These tests never masquerade as real floral adaptation,
climatological rainfall during flowering or calibrated gravity inclinations.
"""
import json

import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_full_rain_head_response import (
    SPECIES, PRIMARY, SECONDARY, analyse_photo_orientation, climate_reference,
    compare_within_continent, _bootstrap_shared_cells, MIN_INTERVAL_CELLS
)
from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
from analysis.v4.native_spatial import cells
from analysis.v4.native_cached_climate import VARIABLES


def fake_native_intro_angles(*,ncell=40,each=4,seed=20261021):
    rng=np.random.default_rng(seed)
    records=[]
    for species in SPECIES:
        for state,cont in (("native",1),("introduced",7),
                           ("introduced",2),("introduced",8)):
            for k in range(ncell):
                for j in range(each):
                    wetter=state=="introduced" and k>=ncell//2
                    low=600.+36*np.sin(k*.28)
                    clim={
                        "chelsa_bio01":float(10+3*np.cos(k*.39)+rng.normal(scale=.15)),
                        "chelsa_bio04":float(780+30*np.sin(k*.18)+rng.normal(scale=2)),
                        "chelsa_bio12":float((920 if wetter else low)+rng.normal(scale=5)),
                        "chelsa_bio15":float((58 if wetter else 30+3*np.sin(k*.25))+rng.normal(scale=.7)),
                        "chelsa_rsds_mean":float(15+1*np.cos(k*.16)+rng.normal(scale=.1)),
                    }
                    response=55.+(8. if wetter else 0.)+rng.normal(scale=.75)
                    row={
                        "obs_id":f"{species}-{state}-{cont}-{k}-{j}",
                        "taxon_name":species,"range_state":state,
                        "user_id":f"observer-{species}-{state}-{cont}-{k}-{j}",
                        "latitude":float(-60+35*(k//20)),
                        "longitude":float(-160+(k%20)*7),
                        "observed_on":f"2024-{(k%12+1):02d}-10",
                        "tdwg1":cont,
                        "orientation_angle_degrees_median":float(response),
                    }
                    row.update(clim)
                    records.append(row)
    result=pd.DataFrame(records)
    result["_cell"]=cells(result.latitude.to_numpy(float),result.longitude.to_numpy(float))
    return result


def test_both_species_novel_wet_and_secondary_seasonality_labels_are_angle_blind():
    f=fake_native_intro_angles()
    a=analyse_photo_orientation(f,bootstrap=100)
    assert a["status"]=="DIRECT_CLIMATE_NOVELTY_AND_IMAGE_ORIENTATION_FINISHED_EXPLORATORY"
    assert set(a["focal"])==set(SPECIES)
    for species,out in a["focal"].items():
        assert out["native_source_observation_count"]==160
        assert set(out["regions"])=={"North America","Africa","South America"}
        assert out["native_reference"][PRIMARY]["angle_ignored_in_threshold"]
        assert out["native_reference"][SECONDARY]["angle_ignored_in_threshold"]
        for continent,region in out["regions"].items():
            for key in ("annual_precipitation","precipitation_seasonality"):
                x=region["hypotheses"][key]
                assert x["status"]=="REGIONAL_DESCRIPTIVE_CONTRAST_SUPPORTED"
                assert x["camera_angle_never_used_in_class_selection"]
                assert x["higher_than_native_q95_minus_native_analogue_same_continent"]["status"]=="EXPLORATORY_REGIONAL_CELL_BOOTSTRAP"
                interval=x["higher_than_native_q95_minus_native_analogue_same_continent"]["site_cell_bootstrap_q025_q50_q975"]
                assert len(interval)==3
                assert interval[0]>0
    assert a["no_original_location_observation_photographer_identifiers_in_public_result"]


def test_thresholds_unchanged_if_head_orientation_is_inverted():
    f=fake_native_intro_angles()
    first=analyse_photo_orientation(f,bootstrap=100)
    inverse=f.copy()
    inverse["orientation_angle_degrees_median"]=180.-inverse["orientation_angle_degrees_median"]
    other=analyse_photo_orientation(inverse,bootstrap=100)
    for species in SPECIES:
        x=first["focal"][species]
        y=other["focal"][species]
        assert x["native_reference"]==y["native_reference"]
        for region in x["regions"]:
            for name in ("annual_precipitation","precipitation_seasonality"):
                a=x["regions"][region]["hypotheses"][name]
                b=y["regions"][region]["hypotheses"][name]
                assert a["threshold_native_cell_equal_q95"]==b["threshold_native_cell_equal_q95"]
                assert a["categories"]["above_native_q95"]["original_observations"]==b["categories"]["above_native_q95"]["original_observations"]
                contrast=a["higher_than_native_q95_minus_native_analogue_same_continent"]
                opposite=b["higher_than_native_q95_minus_native_analogue_same_continent"]
                assert contrast["difference_wet_minus_analogue_image_degrees"]==pytest.approx(
                    -opposite["difference_wet_minus_analogue_image_degrees"]
                )


def test_missing_climate_does_not_generate_false_zero_direction():
    f=fake_native_intro_angles()
    inc=f.copy()
    inc.loc[inc.range_state.eq("introduced"),PRIMARY]=np.nan
    a=analyse_photo_orientation(inc,bootstrap=100)
    for result in a["focal"].values():
        assert result["introduced_source_observation_count"]==0
        assert result["regions"]=={}


def test_rare_regions_dont_get_significant_pointwise_bootstrap():
    f=fake_native_intro_angles()
    balanced=dated_and_balanced(f)
    sp=balanced.loc[balanced.taxon_name.eq(SPECIES[0])]
    native=sp.loc[sp.range_state.eq("native")]
    thresholds=climate_reference(native)
    introduced=sp.loc[sp.range_state.eq("introduced")&sp.tdwg1.eq(2)].copy()
    few=introduced.loc[introduced.chelsa_bio12>thresholds[PRIMARY]["native_cell_equal_q95"]].head(6)
    analogue=introduced.loc[introduced.chelsa_bio12.between(
        thresholds[PRIMARY]["native_cell_equal_q05"],thresholds[PRIMARY]["native_cell_equal_q95"]
    )].head(6)
    result=_bootstrap_shared_cells(few,analogue,name="rare",replicates=100)
    assert result["status"]=="STOP_LESS_THAN_15_CELLS_EACH_GROUP"
    assert result["interval"] is None


def test_privacy_prevents_source_ids_from_aggregate_morphological_claim():
    a=analyse_photo_orientation(fake_native_intro_angles(),bootstrap=100)
    payload=json.dumps(a)
    assert '"obs_id"' not in payload
    assert '"latitude"' not in payload
    assert '"longitude"' not in payload
    assert "observer-Cirsium" not in payload
    assert "2024-05-10" not in payload
    assert "gravity" in payload
