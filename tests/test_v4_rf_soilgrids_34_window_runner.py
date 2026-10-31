"""Window fallback uses only synthetic GeoTIFF, never S1 coordinates."""
import numpy as np
import rasterio
import pytest
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from analysis.v4.rf_soilgrids_34_window_runner import soil_reader_windowed


def test_34_fallback_reader_exact_pixels_mask_and_large_array(monkeypatch):
    rng=np.random.default_rng(20261031)
    grid=rng.integers(1,30000,size=(1024,1024),dtype='int16')
    grid[100:140,100:180]=-9999
    trans=from_origin(-180,80,0.01,0.01)
    rows=rng.integers(0,1024,size=2300)
    cols=rng.integers(0,1024,size=2300)
    lat=trans.f + (rows+.5)*trans.e
    lon=trans.c + (cols+.5)*trans.a
    lat=np.r_[lat,trans.f-(110+.5)*.01]
    lon=np.r_[lon,trans.c+(120+.5)*.01]
    with MemoryFile() as mem:
        with mem.open(driver='GTiff',width=1024,height=1024,count=1,
                      dtype='int16',nodata=-9999,transform=trans,crs='EPSG:4326',
                      tiled=True,blockxsize=256,blockysize=256) as dst:
            dst.write(grid,1)
        with mem.open() as original:
            expect=np.array([np.nan if np.ma.is_masked(p[0]) else float(p[0])
                             for p in original.sample(zip(lon,lat),indexes=1,masked=True)])
        monkeypatch.setattr(rasterio,'open',lambda *args,**kw: mem.open())
        got=soil_reader_windowed('https://files.isric.org/soilgrids/latest/data/soc/soc_0-5cm_Q0.5.vrt',lat,lon)
    assert got.shape==(2301,)
    assert np.allclose(got,expect,equal_nan=True)
    assert np.isnan(got[-1])

@pytest.mark.parametrize('lat,lon',[(np.array([10.,np.nan]),np.array([20.,30.])),
                                          (np.array([90.1]),np.array([0.])),
                                          (np.array([1.,2.]),np.array([1.]))])
def test_fallback_rejects_invalid_private_coordinates_before_any_network(lat,lon):
    with pytest.raises(ValueError,match='INVALID_PRIVATE_SOURCE_WGS84'):
        soil_reader_windowed('https://files.isric.org/soilgrids/latest/data/soc/soc_0-5cm_Q0.5.vrt',lat,lon)