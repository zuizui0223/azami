"""Synthetic, outcome-blind S2 wet-edge design and fail-closed regression tests."""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_spatial import cells
from analysis.v4.native_wet_edge import (
    SPECIES, RAIN, OTHER, evaluate, match_without_angles,
    native_envelope, _features, save,
)


def mock_original() -> tuple[pd.DataFrame,pd.DataFrame]:
    rng=np.random.default_rng(20261021)
    rows=[]
    for s,species in enumerate(SPECIES):
        for state in ("native","introduced"):
            for i in range(320):
                wet=state=="introduced" and i%4==0
                rain=(700+(i%100)*7 if state=="native"
                      else (1480+(i%5)*10 if wet else 1250+(i%7)*18))
                lat=-25.+(i%7)*8.
                lon=-155.+((i//7)%9)*12.
                rows.append({
                    "obs_id":f"source-{s}-{state}-{i}",
                    "taxon_name":species,
                    "range_state":state,
                    "user_id":f"observer-{s}-{state}-{i}",
                    "latitude":lat,
                    "longitude":lon,
                    "_cell":int(cells([lat],[lon])[0]),
                    "tdwg1":1 if state=="native" else 7,
                    "observed_on":"2023-06-15",
                    "quality_grade":"research",
                    "chelsa_bio01":5+rng.normal(),
                    "chelsa_bio04":700+rng.normal(0,50),
                    RAIN:float(rain),
                    "chelsa_bio15":30+rng.normal(0,5),
                    "chelsa_rsds_mean":12+rng.normal(),
                    "orientation_angle_degrees_median":
                        float(np.clip(70+6*wet+rng.normal(0,10),0,180))
                })
    photo=pd.DataFrame(rows)
    # Source climate independent of its phenotype! Source has NO angle column.
    source=photo.drop(columns=["orientation_angle_degrees_median"]).copy()
    return source,photo


def test_source_native_quantiles_are_not_photo_selected():
    source,photo=mock_original()
    n=source.loc[(source.taxon_name==SPECIES[0])&
                 (source.range_state=="native")]
    original=native_envelope(n)
    photo["orientation_angle_degrees_median"]=180-photo[
        "orientation_angle_degrees_median"]
    assert native_envelope(n)==original
    assert original["source_n"]==320
    assert original["BIO12_native_source_q50_q75_q90_q95_q975_mm"][3]>1300


def test_real_wet_edge_angle_distribution_not_population_mean():
    source,photo=mock_original()
    result=evaluate(source,photo)
    assert set(result["focal"])==set(SPECIES)
    for species in SPECIES:
        r=result["focal"][species]["introduced_by_continent"]["North America"]
        assert r["source_wet_edge"]["cells"]>=15
        assert r["matching"]["status"]=="ANGLE_BLIND_MATCH_COMPLETE"
        assert r["conditional_angle_distribution"]["status"]==(
            "EXPLORATORY_MATCHED_DISTRIBUTION_ESTIMATED_NOT_INDEPENDENT_CONFIRMED")
        q=r["conditional_angle_distribution"]
        assert len(q["wet_Q10_Q50_Q90"])==3
        assert len(q["wet_minus_reference_Q10_Q50_Q90"])==3
        assert q["weighted_one_dimensional_wasserstein_W1_degrees"]>=0
        uncertainty=q["geographic_block_uncertainty"]
        assert uncertainty["status"]=="GEOGRAPHIC_CELL_SENSITIVITY_INTERVAL_NOT_CAUSAL_OR_CAMERA_VALIDATED"
        assert len(uncertainty["Q10_Q50_Q90_delta_95pct_intervals"])==3
        assert 0<=uncertainty["Q90_positive_fraction_of_draws"]<=1
        assert uncertainty["not_an_independent_population_confidence_interval"]


def test_matching_is_independent_of_angle_and_bio12_amplitude():
    source,photo=mock_original()
    ref=native_envelope(source.loc[(source.taxon_name==SPECIES[0])&
                        (source.range_state=="native")])
    df=photo.loc[(photo.taxon_name==SPECIES[0])&
                 (photo.range_state=="introduced")].copy()
    df["_month"]=6
    df["_year"]=2023
    wet=df.loc[df[RAIN]>ref["BIO12_native_source_q50_q75_q90_q95_q975_mm"][3]]
    control=df.loc[df[RAIN].between(*[
        ref["BIO12_native_source_q50_q75_q90_q95_q975_mm"][i] for i in (1,3)])]
    a=match_without_angles(wet,control,ref)
    changed=wet.copy()
    changed[RAIN]=changed[RAIN]+9999
    changed["orientation_angle_degrees_median"]=180-changed[
        "orientation_angle_degrees_median"]
    b=match_without_angles(changed,control,ref)
    assert a["status"]=="ANGLE_BLIND_MATCH_COMPLETE"
    assert np.array_equal(a["wet_indices"],b["wet_indices"])
    assert np.array_equal(a["control_indices"],b["control_indices"])
    assert np.allclose(_features(wet,ref),_features(changed,ref))


def test_insufficient_independent_wet_cells_stop_as_unidentified():
    source,photo=mock_original()
    all_wet=photo.loc[(photo.taxon_name==SPECIES[1]) &
                       photo.range_state.eq("introduced") & (photo[RAIN]>1400)]
    keep_cells=sorted(all_wet._cell.unique())[:2]
    selected=(photo.taxon_name==SPECIES[1]) & photo.range_state.eq("introduced") & (
        photo[RAIN]>1400) & ~photo["_cell"].isin(keep_cells)
    photo.loc[selected,"orientation_angle_degrees_median"]=np.nan
    result=evaluate(source,photo)
    r=result["focal"][SPECIES[1]]["introduced_by_continent"]["North America"]
    assert r["conditional_angle_distribution"]["status"]=="NOT_IDENTIFIED_NOT_ZERO"
    assert r["source_wet_edge"]["cells"]>=15  # unchanged, true source!


def test_frozen_source_keys_necessary_and_public_receipt_has_no_identifiers(tmp_path):
    source,photo=mock_original()
    with pytest.raises(ValueError,match="source|climate|fields"):
        evaluate(source.drop(columns=["user_id"]),photo)
    output=save(source,photo,tmp_path)
    f=tmp_path/"S2_wet_edge_angle_blind_distribution.json"
    raw=f.read_text()
    assert json.loads(raw)["status"]==output["status"]
    for token in ('"obs_id"','"user_id"','"latitude"','"longitude"'):
        assert token not in raw
    assert "source_wild_complete_climate_n" in raw
    assert "individual_identifiers" not in raw  # No individuals in the public JSON.
