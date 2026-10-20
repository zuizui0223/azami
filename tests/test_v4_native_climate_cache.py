import numpy as np
import pytest
from analysis.v4.native_climate_cache import tile_window, choose_url_paths, extract_tiles, sample_tiles


def test_source_names_and_unchanged_version():
    urls=choose_url_paths('rsds_mean')
    assert any(u.endswith('CHELSA_rsds_1981-2010_mean_V.2.1.tif') for u in urls)
    assert all('V.2.1' in u for u in urls)
    assert all('1981-2010' in u for u in choose_url_paths('bio01'))
    with pytest.raises(ValueError):choose_url_paths('invented')


def test_exact_tile_and_sampling_no_resampling(tmp_path):
    rasterio=pytest.importorskip('rasterio')
    from rasterio.transform import from_origin
    p=tmp_path/'toy.tif'
    a=np.arange(32*72,dtype='int16').reshape(32,72)
    with rasterio.open(p,'w',driver='GTiff',height=32,width=72,count=1,dtype='int16',
        crs='EPSG:4326',transform=from_origin(-180,90,5,5),nodata=-32768) as dst:dst.write(a,1)
    with rasterio.open(p) as src:
        w=tile_window(src,20,30)
        assert w.width==1 and w.height==1
    out=tmp_path/'cache';extract_tiles(p,[[20,30],[21,30]],out)
    got=sample_tiles(out,np.array([12.2,17.4,0]),np.array([-27.2,-27.4,0]))
    assert got[0]==a[15,30] and got[1]==a[14,30]
    assert np.isnan(got[2])
