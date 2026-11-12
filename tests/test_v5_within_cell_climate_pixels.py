"""Synthetic geographic-pixel tests: no private S1 observations required."""
import numpy as np
import pandas as pd
import pytest

from analysis.v5.within_cell_climate_pixels import (
    CLIMATE, OUTCOME, FOCAL, PIX_DEG, X_LEFT, Y_TOP,
    _climate_pixel_index, _group_demean, build_pixel_cells,
    fit_within_cells, _cluster_moments, _bootstrap_group
)


def fake_photos():
    rows=[]
    rng=np.random.default_rng(189)
    # 3 taxa with 8 populated 1° cells, 6 independent CHELSA pixels per
    # taxon×cell, 2 photos per pixel.
    for species in ("taxon_a","taxon_b","taxon_c"):
        for cell in range(8):
            south=20+cell
            for pixel in range(6):
                lat=south+.11+pixel*.009
                lon=35+.11+pixel*.009
                rain=600+cell*90+rng.normal(0,28)
                temp=8+cell*.6+rng.normal(0,1)
                radiation=15+rng.normal(0,1)
                tempseason=150+rng.normal(0,13)
                rainseason=30+cell*.7+rng.normal(0,5)
                angle=70+.04*(rain-600)+rng.normal(0,.9)+8*(species=="taxon_c")
                chroma=25+.12*(rainseason-30)+rng.normal(0,.2)+4*(species=="taxon_b")
                for dup in range(2):
                    rows.append({
                        "accepted_key":species,"latitude":lat,"longitude":lon,
                        "chelsa_bio12":rain,"chelsa_bio01":temp,
                        "chelsa_rsds_mean":radiation,
                        "chelsa_bio04":tempseason,"chelsa_bio15":rainseason,
                        "orientation_image_vertical_angle":angle,
                        "corolla_lab_chroma":chroma,
                        "capitulum_outline_aspect_ratio":.7,
                        "obs_id":f"{species}_{cell}_{pixel}_{dup}",
                    })
    return pd.DataFrame(rows)


def test_30arcsec_original_pixel_index_is_not_one_degree_grid():
    d=pd.DataFrame({"latitude":[20.15,20.159,21.15],
                    "longitude":[35.15,35.159,35.15]})
    row,col=_climate_pixel_index(d)
    assert row[0]!=row[1] and col[0]!=col[1]
    assert row[0]!=row[2] and col[0]==col[2]
    assert np.isclose(PIX_DEG,1/120,rtol=0,atol=1e-8)


def test_photo_replication_does_not_increase_pixel_support_or_change_slopes():
    d=fake_photos()
    p,s=build_pixel_cells(d,1,min_photos=10,min_pixels=5)
    assert s["supported_taxa"]==3
    assert s["supported_taxon_cells"]==24
    assert s["supported_unique_taxon_cell_pixels"]==144
    assert s["photos_supported_by_independent_pixels"]==288
    baseline=p[list(CLIMATE)+list(OUTCOME)].to_numpy()
    # Two extra copies of the same measured photograph do not create
    # additional geographic independent pixels or change their medians.
    doubled=pd.concat([d,d],ignore_index=True)
    p2,s2=build_pixel_cells(doubled,1,min_photos=10,min_pixels=5)
    assert len(p)==len(p2)
    np.testing.assert_allclose(baseline,p2[list(CLIMATE)+list(OUTCOME)].to_numpy())
    assert s2["photos_supported_by_independent_pixels"]==576


def test_source_taxon_cell_centering_removes_between_cell_intercepts():
    p,s=build_pixel_cells(fake_photos(),1)
    x,y,w,sd=_group_demean(p)
    assert len(x)==s["supported_unique_taxon_cell_pixels"]
    assert x.shape[1]==5 and y.shape[1]==2
    assert np.all(np.isfinite(sd)) and np.all(sd>0)
    np.testing.assert_allclose(x.T@w,np.zeros(5),atol=1e-12)
    # Rain–angle and seasonality–chroma invented positive within each
    # cell, beyond strong between-group differences.
    b=np.linalg.solve(x.T@(w[:,None]*x),x.T@(w[:,None]*y))
    assert b[0,0]>.15
    assert b[4,1]>.15


def test_region_cluster_and_taxon_cluster_bootstrap_no_pseudo_replicates():
    p,s=build_pixel_cells(fake_photos(),1)
    x,y,w,_=_group_demean(p)
    xt,yt,ids=_cluster_moments(x,y,w,p._taxon)
    assert len(ids)==3
    assert np.isclose(xt.sum(axis=0)[0,0],1.0)
    one,multi=_bootstrap_group(xt,yt,123,nboot=20)
    # Fewer than 10 source clusters is deliberately not used to
    # certify biological inference, though moments are still available.
    assert len(one)==0
    assert len(multi)==0


def test_same_raster_pixel_must_not_have_conflicting_source_climate():
    d=fake_photos()
    dup=d.iloc[0].copy()
    dup["chelsa_bio12"]+=100
    with pytest.raises(ValueError,match="INCONSISTENT_ORIGINAL_RASTER_PIXEL_CLIMATE"):
        build_pixel_cells(pd.concat([d,pd.DataFrame([dup])],ignore_index=True),1)


def test_only_two_target_hypotheses_and_original_five_climates():
    assert len(FOCAL)==2
    assert FOCAL[0][1:] == ("chelsa_bio12","orientation_image_vertical_angle")
    assert FOCAL[1][1:] == ("chelsa_bio15","corolla_lab_chroma")
    assert set(CLIMATE)=={"chelsa_bio12","chelsa_bio01","chelsa_rsds_mean",
                          "chelsa_bio04","chelsa_bio15"}
