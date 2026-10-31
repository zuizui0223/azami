"""Synthetic-only DEM geometry; never access original thistle coordinates."""
import numpy as np
import pandas as pd
import pytest

from analysis.v4.rf_worldclim_terrain import (
    _terrain_from_pixels, sample_dem_terrain, MODEL_SET, TERRAIN5,
)
from analysis.v4.rf_multitrait_predictability import sha256


def test_model_has_five_distinct_terrain_features():
    left, right = list(MODEL_SET.values())
    assert len(left) == 23 and len(right) == 28
    assert set(right)-set(left) == set(TERRAIN5)


def test_planar_surface_known_slope_tpi_and_downhill_direction():
    # increase elevation towards NORTH, equally flat in E-W: downhill = south
    rr,cc = np.meshgrid(np.arange(70), np.arange(70), indexing='ij')
    z = 2000 - 8 * rr
    pts_r=np.array([20, 25, 30]);pts_c=np.array([20,30,35])
    result,center,neighbors = _terrain_from_pixels(z,pts_r,pts_c,
                                                    np.array([45.]*3),
                                                    1/120,1/120)
    assert np.isfinite(result).all()
    assert np.max(np.abs(result[:,2])) < 1e-6  # linear plane: TPI 0
    assert (result[:,0] > 0).all()
    assert np.allclose(result[:,1],np.sqrt((8**2*6)/8))  # 3 northern + 3 southern neighbours, 2 lateral zero
    assert np.allclose(result[:,3],-1.,atol=1e-6)  # downhill to south
    assert np.allclose(result[:,4],0.,atol=1e-6)
    assert np.all(neighbors==120)
    assert np.all(center==z[pts_r,pts_c])


def test_flat_and_ocean_mask_propagate():
    z = np.full((60,60),100.,dtype=float)
    z[10,10] = np.nan
    r=np.array([15,35, 9]);c=np.array([15,35,9])
    metrics,_,support=_terrain_from_pixels(z,r,c,np.array([35.,35.,35.]),1/120,1/120)
    assert (metrics[1]==0).all()
    assert np.isnan(metrics[2]).all()  # within 3x3 missing neighbor
    assert support[0]==119


def test_real_tif_source_sha_and_tile_borders(tmp_path):
    rasterio=pytest.importorskip('rasterio')
    from rasterio.transform import from_origin
    n=300
    yy, xx = np.meshgrid(np.arange(n),np.arange(n),indexing='ij')
    data=(3000+3*xx-2*yy).astype('float32')
    tif=tmp_path/'dem.tif'
    with rasterio.open(tif,'w',driver='GTiff',width=n,height=n,count=1,
                       dtype='float32',crs='EPSG:4326',nodata=-9999,
                       transform=from_origin(90,40,1/120,1/120)) as dst:
        dst.write(data,1)
    # spans the 256x256 tile boundary
    rr=np.array([12,250,255,260,290]);cc=np.array([12,250,255,260,290])
    lat=40-(rr+0.5)/120;lon=90+(cc+0.5)/120
    values,center,rec=sample_dem_terrain(tif,sha256(tif),lat,lon,require_global=False)
    assert np.isfinite(values.to_numpy()).all()
    assert np.allclose(center,data[rr,cc])
    assert rec['finite_all_terrain_records']==5
    # Manifests require meaningful, non-placeholder strings. A one-character
    # SI symbol "m" is scientifically valid but rejected by the source gate;
    # the published contract uses the explicit label "metres".
    assert all(len(rec['unit_by_variable'][name]) > 2 for name in TERRAIN5)
    assert rec['unit_by_variable']['topo_tpi'] == 'metres'
    with pytest.raises(ValueError,match='SHA_MISMATCH'):
        sample_dem_terrain(tif,'0'*64,lat,lon,require_global=False)


def test_invalid_dem_grain_refused(tmp_path):
    tif=tmp_path/'bad.tif';tif.write_text('not a raster')
    with pytest.raises(Exception):
        sample_dem_terrain(tif,sha256(tif),np.array([20]),np.array([90]),require_global=False)