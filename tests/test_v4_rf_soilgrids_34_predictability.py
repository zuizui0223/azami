"""NO live soil or original observations; test sampling and no-leak contracts only."""
import hashlib
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4.rf_soilgrids_34_predictability import (
    MODEL28, MODEL34, MODELS, SOIL6, stable_hash,
    selected_source_positions, sample_six_source_sites,
    attach_soil_to_cells, score_models,
)
from analysis.v4.rf_soilgrids_webdav_probe import CONVERSIONS


def synthetic_source(n_taxa=30, n_cells=14, photos=5):
    rng=np.random.default_rng(83)
    rows=[]
    for t in range(n_taxa):
        for cell in range(n_cells):
            lon=(-110, 0, 110)[t%3] + cell*1.1 + .2
            lat=-45 + cell*1.1 + .2
            for i in range(photos):
                data={n: float(rng.normal()) for n in MODEL28}
                rows.append({'obs_id':f'obs_{t}_{cell}_{i}',
                             'accepted_key':f'species_{t}',
                             'latitude':lat+i/10000,
                             'longitude':lon+i/10000,
                             'orientation_image_vertical_angle':float(rng.normal(50,4)),
                             'orientation_image_vertical_angle__eligible':True,
                             **data})
    return pd.DataFrame(rows)


def fake_vrt(code):
    return (f'<VRTDataset rasterXSize="40" rasterYSize="25">'
            '<SRS>EPSG:152160</SRS><GeoTransform>0,250,0,0,0,-250</GeoTransform>'
            '<VRTRasterBand dataType="Int16" band="1">'
            f'<SimpleSource><SourceFilename relativeToVRT="1">{code}_0-5cm_Q0.5/tile000.tif'
            '</SourceFilename></SimpleSource></VRTRasterBand></VRTDataset>').encode()


def test_registration_and_fixed_predictor_scope():
    assert len(MODEL28)==28 and len(MODEL34)==34
    assert set(MODEL34)-set(MODEL28)==set(SOIL6)
    assert set(MODELS)=={'terrain28_soilmatched','terrain28_plus_soil34'}


def test_selected_one_degree_sites_deterministic_hash_and_no_leak():
    s=synthetic_source()
    a,x,loc,report=selected_source_positions(s)
    b,y,loc2,report2=selected_source_positions(s)
    assert a.equals(b)
    assert x.equals(y)
    assert loc.equals(loc2)
    assert report==report2
    assert report['n_outcome_blind_selected_taxa']==30
    assert report['n_selected_taxon_cells']==30*14
    assert report['n_selected_source_photo_positions']==30*14*3
    assert loc.obs_id.nunique()==len(loc)
    assert x.groupby(['taxon','lat_cell','lon_cell']).size().eq(3).all()
    assert 'species_0' not in json.dumps(report)
    assert 'obs_0' not in json.dumps(report)
    # Angles changed without eligibility change: selection must be invariant.
    s['orientation_image_vertical_angle']+=np.linspace(-3,3,len(s))
    _,x3,_,q=selected_source_positions(s)
    assert x.equals(x3) and q==report


def test_missingness_two_of_three_without_outcome_cherry_pick():
    s=synthetic_source()
    cells,photo,locations,_=selected_source_positions(s)
    good=locations[['obs_id']].copy()
    for k in SOIL6:good[k]=1.
    result,support,elig=attach_soil_to_cells(cells,photo,good)
    assert elig and support['n_taxa_after_soil_requalification']==30
    assert len(result)==len(cells)
    assert not set(SOIL6).difference(result.columns)
    assert not set(['obs_id','latitude','longitude']).intersection(result.columns)
    # Lose 1 pixel in every group; two others suffice.
    first=photo.groupby(['taxon','lat_cell','lon_cell']).head(1).obs_id
    with_one_missing=good.copy()
    with_one_missing.loc[with_one_missing.obs_id.isin(first),list(SOIL6)]=np.nan
    result,support,elig=attach_soil_to_cells(cells,photo,with_one_missing)
    assert elig and len(result)==len(cells)
    # Lose two photo positions in each group => no cell should pass.
    first_two=photo.groupby(['taxon','lat_cell','lon_cell']).head(2).obs_id
    bad=good.copy()
    bad.loc[bad.obs_id.isin(first_two),list(SOIL6)]=np.nan
    result,support,elig=attach_soil_to_cells(cells,photo,bad)
    assert not elig and support['n_taxa_after_soil_requalification']==0


def test_original_vrt_sha_pin_rejects_fake_remote_soil_bytes():
    sites=pd.DataFrame({'obs_id':[f'x{i}' for i in range(30)],
                        'latitude':np.ones(30)*43.,'longitude':np.ones(30)*-80.})
    def bad(url):return b'<html>pretend source</html>',{}
    with pytest.raises(ValueError,match='VRT_SHA_MISMATCH'):
        sample_six_source_sites(sites,fetcher=bad,
                                reader=lambda u,lat,lon: np.ones(len(lat)))


def test_realistic_six_source_pixels_and_full_pair_fit(monkeypatch):
    import analysis.v4.rf_soilgrids_34_predictability as mod
    for k in SOIL6:
        code=CONVERSIONS[k][0]
        monkeypatch.setitem(mod.VRT_PINNED_SHA256,k,hashlib.sha256(fake_vrt(code)).hexdigest())
    def fetch(url):return fake_vrt(url.split('/')[-1].split('_')[0]),{}
    def reader(url,lat,lon):
        code=url.split('/')[-1].split('_')[0]
        physical={'phh2o':6.3,'soc':25,'clay':30,'sand':52,'bdod':1.1,'cec':20}[code]
        return np.full(len(lat),physical*CONVERSIONS['soil_'+code+'_0_5cm'][1],dtype=float)
    cells,photo,sites,_=selected_source_positions(synthetic_source())
    soil,provenance=sample_six_source_sites(sites,fetcher=fetch,reader=reader)
    assert len(provenance['variables'])==6
    assert provenance['n_six_layer_complete_original_photo_locations']==len(sites)
    merged,support,elig=attach_soil_to_cells(cells,photo,soil)
    assert elig and support['n_taxa_after_soil_requalification']==30
    for scheme in ('unseen_taxa_balanced','unseen_longitude_sectors'):
        fit=score_models(merged,'angle',scheme,trees=8,repeats=1,
                         predictor_sets=MODELS,predictor_groups=mod.GROUPED_IMPORTANCE)
        assert fit['n_taxa']==30
        assert set(fit['models'])==set(MODELS)
        assert all(x['n_folds'] in (3,5) for x in fit['models'].values())
        assert 'species_' not in str(fit)
        assert 'obs_' not in str(fit)


def test_budget_rejects_too_many_selected_source_positions():
    with pytest.raises(ValueError,match='SOURCE_POINT_BUDGET'):
        selected_source_positions(synthetic_source(n_taxa=80,n_cells=14,photos=5))

def test_real_geotiff_vectorized_pixel_index_regression(tmp_path, monkeypatch):
    """The source-wide 3K point call must work with array coordinates."""
    import rasterio
    from rasterio.transform import from_origin
    from analysis.v4.rf_soilgrids_34_predictability import soil_reader
    path=tmp_path/'soil-grid-mini.tif'
    with rasterio.open(path,'w',driver='GTiff',count=1,crs='EPSG:4326',
                       height=4,width=4,dtype='uint16',nodata=65535,
                       transform=from_origin(0,4,1,1)) as writer:
        writer.write(np.arange(1,17,dtype=np.uint16).reshape(4,4),1)
    genuine_open=rasterio.open
    monkeypatch.setattr(rasterio,'open',lambda *_args,**_kwargs:genuine_open(path))
    result=soil_reader('https://files.isric.org/soilgrids/latest/data/phh2o/phh2o_0-5cm_Q0.5.vrt',
                       np.array([3.5,2.5,1.5]),np.array([0.5,1.5,2.5]),
                       max_seconds=10)
    assert result.tolist()==[1.0,6.0,11.0]
