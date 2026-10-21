"""PCA-env-style source-availability science contracts; synthetic sources only.

Do not conflate a TDWG WCVP native/introduced environment with the species'
realized range, a raw occurrence histogram, or ecospat's calibrated niche
dynamic tests. Reject impossible source samples instead of fabricating niches.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_pca_env import (
    BINS, _pca_env_2d, _weighted_projected_pca, _probability_content,
    _native_and_introduced_codes, background_mesh, climatic_shift_summary,
    SPECIES_KEYS, VARIABLES
)


def fake_regions():
    # Two nonoverlapping 15x20° TDWG pseudo-regions.
    def feature(code,group,x0,x1):
        return {"type":"Feature","geometry":{
            "type":"Polygon","coordinates":[[
                [x0,-15],[x1,-15],[x1,15],[x0,15],[x0,-15]
            ]]},
            "properties":{"LEVEL3_COD":code,"LEVEL1_COD":group}}
    return {
        "geometry":{"features":[feature("NA1",1,-140,-120),
                                 feature("IN1",7,-110,-90)]},
        "distributions":{
            str(SPECIES_KEYS["Cirsium arvense"]):[
                {"locationId":"TDWG:NA1","occurrenceStatus":"PRESENT","establishmentMeans":"NATIVE"},
                {"locationId":"TDWG:IN1","occurrenceStatus":"PRESENT","establishmentMeans":"INTRODUCED"}],
            str(SPECIES_KEYS["Cirsium vulgare"]):[
                {"locationId":"TDWG:NA1","occurrenceStatus":"PRESENT","establishmentMeans":""},
                {"locationId":"TDWG:IN1","occurrenceStatus":"PRESENT","establishmentMeans":"INTRODUCED"}]
        }
    }


def fake_values(n, *, climate_shift=0., seed=20261021):
    rng=np.random.default_rng(seed)
    c=rng.normal(size=(n,len(VARIABLES)))
    c[:,0]+=climate_shift
    tab=pd.DataFrame(c,columns=VARIABLES)
    tab["_area_weight"]=np.cos(np.deg2rad(rng.uniform(-15,15,n)))
    tab["_cell"]=np.arange(n)%50
    tab["obs_id"]=[f"synth-{seed}-{j}" for j in range(n)]
    return tab


def test_species_WCVP_region_masks_are_climate_and_photo_outcome_blind():
    a=fake_regions()
    assert _native_and_introduced_codes(
        a,str(SPECIES_KEYS["Cirsium arvense"]))=={
            "NA1":"native","IN1":"introduced"
        }
    background=background_mesh(a,step=1.)
    assert set(background.TDWG_code)=={"NA1","IN1"}
    assert len(background)>600
    assert set(background.TDWG_level1)=={1,7}
    assert background._area_weight.between(0,1).all()
    assert background._cell_key.str.contains("_").all()
    a["distributions"][str(SPECIES_KEYS["Cirsium arvense"])].append(
        {"locationId":"TDWG:IN1","occurrenceStatus":"PRESENT",
         "establishmentMeans":"NATIVE"}
    )
    assert "IN1" not in _native_and_introduced_codes(
        a,str(SPECIES_KEYS["Cirsium arvense"])
    )


def test_pca_background_projection_is_one_shared_climate_system():
    env_n=fake_values(900,seed=1)
    env_i=fake_values(700,climate_shift=2,seed=2)
    o_n=fake_values(400,seed=3)
    o_i=fake_values(350,climate_shift=2,seed=4)
    d=_pca_env_2d(env_n,env_i,o_n,o_i)
    assert d["status"]=="EXECUTED_PCA_ENV_STYLE_CACHED_BACKGROUND_NOT_EXACT_R_ECOSPAT"
    assert len(d["pca_axes_explained_background_variance"])==2
    assert 0<sum(d["pca_axes_explained_background_variance"])<=1
    assert set(d["two_climate_PC_loadings"])==set(VARIABLES)
    assert 0<=d["Schoener_D_available_environment_corrected_density"]<=1
    assert 0<=d["Schoener_D_uncorrected_photo_observation_density"]<=1
    assert 0<=d["introduced_probability_outside_native_95pct_probability_content"]<=1
    assert d["not_ecospat_exact_implementation"]
    assert d["range_native_projection"]["photographed_occurrences_projected_inside_common_grid"]>.8


def test_climate_PCA_env_does_not_read_orientation_or_photo_colour():
    en=fake_values(700,seed=50)
    ei=fake_values(700,climate_shift=1.4,seed=51)
    n=fake_values(300,seed=52)
    i=fake_values(300,climate_shift=1.4,seed=53)
    a=_pca_env_2d(en,ei,n,i)
    n["orientation_angle_degrees_median"]=0.
    i["orientation_angle_degrees_median"]=180.
    n["corolla_lab_chroma"]=40.
    i["corolla_lab_chroma"]=0.
    b=_pca_env_2d(en,ei,n,i)
    assert a==b
    summary=climatic_shift_summary(n,i)
    assert set(summary)==set(VARIABLES)
    assert summary["chelsa_bio01"]["introduced_minus_native_median_original_units"]>0


def test_PCA_env_exact_missing_environment_stops_not_zero_shift():
    a=fake_values(600,seed=5)
    b=fake_values(600,seed=6)
    n=fake_values(300,seed=7)
    i=fake_values(300,seed=8)
    x=_pca_env_2d(a.head(100),b,n,i)
    assert x["status"]=="STOP_BACKGROUND_INSUFFICIENT_OR_RASTER_COVERAGE"
    n.loc[0,"chelsa_bio15"]=np.nan
    with pytest.raises(ValueError,match="non-finite|Invalid"):
        _pca_env_2d(a,b,n,i)


def test_probability_content_is_true_95pct_density_mass_not_raw_species_fraction():
    z=np.array([[.85,.1],[.03,.02]])
    region=_probability_content(z,.95)
    assert region.shape==(2,2)
    assert region[0,0] and region[0,1]
    assert not region[1,1]
    assert abs(z[region].sum()-.95)<1e-12
    with pytest.raises(ValueError,match="not normalized"):
        _probability_content(z*3,.95)
