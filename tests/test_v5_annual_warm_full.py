"""Source-free checks for V5 annual/warm + seasonality comparison."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v5.annual_warm_s1 import (
    CLIMATE, CLIMATE_NAMES, EXTRA, THRESH, _bootstrap_moments
)
from analysis.v5.annual_warm_s2 import (
    MODELS, full_climate_pca, REFERENCE_ANNUAL_D
)
from analysis.v5.annual_warm_sources import (
    EXTRA as WARM_EXTRA, RECEIPT
)
from analysis.v4.s1_climate_trait_hypervolume import AXES, PHENOTYPE


def test_pairwise_five_axis_contract():
    assert set(CLIMATE)==set(MODELS)
    assert len(MODELS)==2
    assert MODELS["annual_plus_seasonality"]==(
        "chelsa_bio01","chelsa_bio04","chelsa_bio12",
        "chelsa_bio15","chelsa_rsds_mean")
    assert MODELS["warm_plus_seasonality"]==(
        "chelsa_bio05","chelsa_bio04","chelsa_bio18",
        "chelsa_bio15","chelsa_rsds_mean")
    assert set(MODELS["annual_plus_seasonality"])&set(
        MODELS["warm_plus_seasonality"])=={
            "chelsa_bio04","chelsa_bio15","chelsa_rsds_mean"}
    assert set(EXTRA[2:])==set(WARM_EXTRA)
    assert len(PHENOTYPE)==3
    assert len(AXES)==6
    assert np.isfinite(THRESH) and THRESH>1


def test_warm_raster_sha_and_physical_units_pinned_in_v4_receipt():
    assert RECEIPT.is_file()
    receipts={x["variable"]:x for x in json.loads(RECEIPT.read_text())["sources"]}
    assert set(WARM_EXTRA)<=set(receipts)
    for var in WARM_EXTRA:
        entry=receipts[var]
        assert len(entry["download_sha256"])==64
        assert entry["source_dtype"]=="uint16"
        assert entry["source_nodata"]==65535
        assert entry["scale"]==.1
        assert entry["sampled_from_original_raster_pixels"]
        assert entry["crs"]=="EPSG:4326"
        assert entry["width"]==43200 and entry["height"]==20880


def test_weighted_pca_comparable_and_normalised_in_any_variable_order():
    rng=np.random.default_rng(41)
    n=400
    orig=rng.normal(size=(n,5))
    frame=pd.DataFrame(orig,columns=MODELS["annual_plus_seasonality"])
    frame["range_state"]=["native"]*(n//2)+["introduced"]*(n//2)
    frame["_area_weight"]=np.linspace(.4,1.3,n)
    native_proj,v=full_climate_pca(frame,MODELS["annual_plus_seasonality"])
    pts=native_proj(frame)
    assert pts.shape==(n,2)
    assert np.isfinite(pts).all()
    assert 0<v<1


def test_temporal_model_comparison_preserves_shared_phenotype_axes():
    annual=CLIMATE["annual_plus_seasonality"]
    warm=CLIMATE["warm_plus_seasonality"]
    assert annual[2:]==warm[2:]
    assert annual[:2]!=warm[:2]
    assert annual[0]=="_z0" and warm[0]=="_z8"
    assert annual[1]=="_z1" and warm[1]=="_z9"


def test_taxon_block_bootstrap_pairwise_moments():
    rng=np.random.default_rng(32)
    x=rng.normal(size=(3,8,5))
    y=rng.normal(size=(3,8,3))
    taxa_i=np.array([0,0])
    taxa_j=np.array([1,2])
    bootstrap=np.ones((5,3))
    result=_bootstrap_moments(x,y,taxa_i,taxa_j,bootstrap)
    assert set(result)=={"between_pair","within_pair"}
    for xt,yt in result.values():
        assert xt.shape==(5,25)
        assert yt.shape==(5,15)
        assert np.isfinite(xt).all()
        assert np.isfinite(yt).all()


def test_annual_known_overlap_is_pinned_only_for_annual_definition():
    assert set(REFERENCE_ANNUAL_D)=={
        "Cirsium arvense","Cirsium vulgare"}
    assert REFERENCE_ANNUAL_D["Cirsium arvense"]==.474
    assert REFERENCE_ANNUAL_D["Cirsium vulgare"]==.442


def test_climate_axis_names_include_both_seasonality_terms():
    for name,axes in CLIMATE_NAMES.items():
        assert len(axes)==5
        assert axes[2:]==("RSDS","BIO04","BIO15")
        assert axes[0] in {"BIO12","BIO18"}
        assert axes[1] in {"BIO01","BIO05"}
