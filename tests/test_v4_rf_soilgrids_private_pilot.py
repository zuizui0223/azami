"""Synthetic-only tests. No actual S1 data, locations or real soil maps."""
import hashlib
import json
import numpy as np
import pandas as pd
import pytest
from analysis.v4.rf_soilgrids_private_pilot import (
    VRT_PINNED_SHA256, SOIL6, source_sentinels,
    audit_pilot, validated_roster,
)
from analysis.v4.rf_soilgrids_webdav_probe import CONVERSIONS
from analysis.v4.rf_multitrait_predictability import sha256


def synthetic_source(n=100):
    rows=[]
    for s, lon in enumerate((-110.,15.,130.)):
        for i in range(n):
            rows.append({'obs_id':f'source_{s}_{i}',
                         'latitude':30.+i/1000,
                         'longitude':lon+i/1000})
    return pd.DataFrame(rows)


def fake_vrt(code):
    return (f'<VRTDataset rasterXSize="40" rasterYSize="25">'
            '<SRS>EPSG:152160</SRS><GeoTransform>0,250,0,0,0,-250</GeoTransform>'
            '<VRTRasterBand dataType="Int16" band="1">'
            f'<SimpleSource><SourceFilename relativeToVRT="1">{code}_0-5cm_Q0.5/tile000.tif'
            '</SourceFilename></SimpleSource></VRTRasterBand></VRTDataset>').encode()


def test_deterministic_outcome_blind_three_sector_cohort():
    df=synthetic_source()
    a,rec=source_sentinels(df,per_sector=12)
    b,_=source_sentinels(df,per_sector=12)
    assert rec['sample_count']==36 and rec['source_sector_observations']=={'0':12,'1':12,'2':12}
    assert a.obs_id.tolist()==b.obs_id.tolist()
    assert a[['latitude','longitude']].drop_duplicates().shape[0]==36
    assert 'source_0_0' not in str(rec)
    with pytest.raises(ValueError,match="BUDGET"):
        source_sentinels(df,per_sector=33)


def test_roster_sha_decode_and_unique_obs_gate(tmp_path):
    df=synthetic_source()
    src=tmp_path/'original.csv.gz'; df.to_csv(src,index=False,compression='gzip')
    receipt=tmp_path/'receipt.json'
    r={'prepared_table_sha256':sha256(src),'native_observations':len(df),
       'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}
    receipt.write_text(json.dumps(r))
    assert len(validated_roster(src,receipt,expected_n=len(df)))==len(df)
    r['prepared_table_sha256']='0'*64;receipt.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="SHA_DECODE"):
        validated_roster(src,receipt,expected_n=len(df))


def test_original_vrt_sha_gate_rejects_any_impersonation():
    df=synthetic_source()
    def fake_fetch(url):
        code=url.split('/')[-1].split('_')[0]
        return fake_vrt(code),{}
    o=audit_pilot(df,per_sector=8,fetcher=fake_fetch,
        pixel_reader=lambda u,lat,lon:np.ones(len(lat)))
    assert o['status']=='insufficient_or_inconsistent_soil_sample'
    assert o['complete_six_original_pixel_pilot_n']==0
    assert all(v['failure_code']=='SOILGRIDS_SOURCE_VRT_SHA_CHANGED' for v in o['layers'].values())
    assert 'source_0_0' not in json.dumps(o)


def test_realistic_original_pin_and_sample_coverage_with_mock(monkeypatch):
    import analysis.v4.rf_soilgrids_private_pilot as mod
    for k in SOIL6:
        code=CONVERSIONS[k][0]
        monkeypatch.setitem(mod.VRT_PINNED_SHA256,k,hashlib.sha256(fake_vrt(code)).hexdigest())
    def fetch(url):
        return fake_vrt(url.split('/')[-1].split('_')[0]),{}
    def values(url, lat, lon):
        name=url.split('/')[-1].split('_')[0]
        d={'phh2o':6.5,'soc':25,'clay':32,'sand':45,'bdod':1.2,'cec':20}
        v=np.full(len(lat),d[name]*CONVERSIONS['soil_'+name+'_0_5cm'][1],dtype=float)
        v[0]=np.nan
        return v
    result=audit_pilot(synthetic_source(),per_sector=6,fetcher=fetch,pixel_reader=values)
    assert result['status']=='pilot_source_pixels_validated_not_full_cohort'
    assert result['complete_six_original_pixel_pilot_n']==17
    assert result['six_soil_properties_covered']==6
    assert 'latitude' not in json.dumps(result)
    assert 'source_0_0' not in json.dumps(result)


def test_no_silent_fill_for_invalid_soil_scale_or_missing(monkeypatch):
    import analysis.v4.rf_soilgrids_private_pilot as mod
    for k in SOIL6:
        code=CONVERSIONS[k][0]
        monkeypatch.setitem(mod.VRT_PINNED_SHA256,k,hashlib.sha256(fake_vrt(code)).hexdigest())
    result=audit_pilot(synthetic_source(),per_sector=4,
        fetcher=lambda u:(fake_vrt(u.split('/')[-1].split('_')[0]),{}),
        pixel_reader=lambda u,lat,lon:np.full(len(lat),np.nan))
    assert result['status']=='insufficient_or_inconsistent_soil_sample'
    assert result['complete_six_original_pixel_pilot_n']==0
    assert all(x['finite_original_pixel_n']==0 for x in result['layers'].values())