"""S2 simple rainfall occupancy followed by image-angle profiles.

Synthetic values: test population and observation identity, source-only rain
breakpoints, original-angle stratification, low-cell stop and privacy.
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_spatial import cells
from analysis.v4.s2_rainfall_occupancy_trait import (
    SPECIES,BANDS,RAIN,assess,save,_band,_source_rain_bands,
    _cell_bootstrap_contrast
)


def example():
    r=np.random.default_rng(20261021)
    rows=[]
    for s,name in enumerate(SPECIES):
        for state in ("native","introduced"):
            for i in range(440):
                # Native source includes genuinely broad rain variation,
                # introduced has greater upper-band occupancy.
                rain=(500.+(i%100)*10 if state=="native"
                      else 780.+(i%100)*11)
                lat=-60+(i%15)*7.
                lon=-168+((i//15)%20)*14.
                rows.append({
                    "obs_id":f"ob-{s}-{state}-{i}",
                    "taxon_name":name,
                    "range_state":state,
                    "user_id":f"photographer-{s}-{state}-{i}",
                    "_cell":int(cells([lat],[lon])[0]),
                    "tdwg1":1 if state=="native" else 7,
                    "latitude":float(lat),
                    "longitude":float(lon),
                    "observed_on":"2024-07-10",
                    RAIN:float(rain),
                    "chelsa_bio01":8+r.normal(),
                    "chelsa_bio04":700+r.normal(0,50),
                    "chelsa_bio15":35+r.normal(0,7),
                    "chelsa_rsds_mean":12+r.normal(),
                    "orientation_angle_degrees_median":float(
                        np.clip(60+(rain>1400)*12+r.normal(0,15),0,180))
                })
    photo=pd.DataFrame(rows)
    source=photo.drop(columns=["orientation_angle_degrees_median"])
    return source,photo


def test_both_species_source_shifts_and_trait_profiles_are_separate():
    source,photo=example()
    o=assess(source,photo)
    assert set(o["focal"])==set(SPECIES)
    for name in SPECIES:
        f=o["focal"][name]
        native=f["climate_occupancy"]["native"]
        introduced=f["climate_occupancy"]["introduced_all"]
        assert introduced["BIO12_median_mm_cell_balanced"]>native["BIO12_median_mm_cell_balanced"]
        for status in ("native","introduced_all"):
            freq=f["climate_occupancy"][status]["bands"]
            assert set(freq)==set(BANDS)
            assert sum(b["source_n"] for b in freq.values())==440
            assert abs(sum(b["source_fraction_cell_equal_weight"] for b in freq.values())-1)<1e-10
            assert set(f["orientation_by_native_rain_band"][status])==set(BANDS)
        assert "North America" in f["introduced_by_continent"]


def test_native_reference_and_climate_bands_are_truly_angle_blind():
    source,photo=example()
    native=source.loc[(source.taxon_name==SPECIES[0]) &
                      (source.range_state=="native")]
    q=_source_rain_bands(native)
    inversion=photo.copy()
    inversion["orientation_angle_degrees_median"]=180-inversion[
        "orientation_angle_degrees_median"]
    a=assess(source,photo)["focal"][SPECIES[0]]
    b=assess(source,inversion)["focal"][SPECIES[0]]
    assert a["native_rain_reference_BIO12_q50_q75_q95_mm"]==list(q)
    assert b["native_rain_reference_BIO12_q50_q75_q95_mm"]==list(q)
    assert a["climate_occupancy"]==b["climate_occupancy"]
    assert (_band(np.array([q[0]-1,q[0],q[1],q[2]]),q).tolist()
            ==list(BANDS))


def test_minimum_spatial_independence_and_non_identified_not_zero():
    src,photo=example()
    few=photo.iloc[:40].copy()
    one=photo.iloc[200:240].copy()
    few["_cell"]=23
    one["_cell"]=24
    z=_cell_bootstrap_contrast(few,one,seed=10)
    assert z["status"]=="NOT_IDENTIFIED_INDEPENDENT_CELLS"
    bad=src.drop(columns=["obs_id"])
    with pytest.raises(ValueError,match="source|provenance|observation"):
        assess(bad,photo)


def test_aggregate_privacy_and_no_plant_genetics_claim(tmp_path):
    src,photo=example()
    o=save(src,photo,tmp_path)
    t=(tmp_path/"S2_rain_occupancy_vs_head_angle_bands.json").read_text()
    assert json.loads(t)["status"]==o["status"]
    for token in ('"obs_id"','"user_id"','"latitude"','"longitude"'):
        assert token not in t
    assert all(v is False for v in o["privacy"].values())
    assert "not gravity" in o["trait_response"]
