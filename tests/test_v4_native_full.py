"""Native-first master-frame invariants; no image/trait result enters inclusion."""
import importlib
import importlib.util
import numpy as np
import pandas as pd
import pytest


def module():
    assert importlib.util.find_spec('analysis.v4.native_full') is not None, 'native-full executor not implemented'
    return importlib.import_module('analysis.v4.native_full')


def source():
    return {'resolutions':[{'input_name':'Cirsium example','resolution_status':'resolved_unique_accepted_key','accepted_key':1,'accepted_name':'Cirsium example A.'}], 'distributions':{'1':[{'locationId':'TDWG:NAT','establishmentMeans':'','occurrenceStatus':'PRESENT'}, {'locationId':'TDWG:INT','establishmentMeans':'INTRODUCED'}]},'geometry':{'features':[{'properties':{'LEVEL3_COD':key},'geometry':{'type':'Polygon','coordinates':[[[x,0],[x+2,0],[x+2,2],[x,2],[x,0]]]}} for x,key in [(0,'NAT'),(3,'INT')]]}}


def roster():
    return pd.DataFrame({'obs_id':['n','i','unlisted','bad','genus'],'taxon_name':['Cirsium example']*4+['Cirsium'],'taxon_rank':['species']*4+['genus'],'latitude':[1,1,10,np.nan,1],'longitude':[1,4,10,np.nan,1],'coordinate_usable_for_environment':[True]*5,'obscured':[False]*5,'captive':[False]*5,'geoprivacy':['']*5})


def test_all_source_rows_retained_and_only_native_is_primary():
    a=module().classify_all(roster(),source())
    assert len(a)==5 and a.obs_id.is_unique
    assert a.loc[a.native_primary,'obs_id'].tolist()==['n']
    assert a.set_index('obs_id').loc['i','native_range_status']=='introduced'
    assert a.set_index('obs_id').loc['genus','eligibility_reason']=='rank_not_species'
    assert not a.set_index('obs_id').loc['unlisted','native_primary']


def test_uncaptured_species_must_stop_not_be_silently_dropped():
    a=roster();a.loc[0,'taxon_name']='Cirsium unknown'
    with pytest.raises(ValueError,match='not completely queried'):
        module().classify_all(a,source())


def test_non_native_status_never_promoted_by_climate_or_phenotype():
    a=roster();a['chroma']=100;a['yolo_score']=1;a['bio01']=20
    b=module().classify_all(a,source())
    assert b.native_primary.sum()==1


def test_private_locations_excluded_even_if_other_flags_are_inconsistent():
    a=roster();a.loc[0,'geoprivacy']='private'
    b=module().classify_all(a,source())
    assert b.native_primary.sum()==0


def test_observation_duplicates_identical_allowed_conflicts_rejected():
    a=roster();same=pd.concat([a,a.iloc[[0]]])
    assert len(module().unique_roster(same))==5
    bad=a.iloc[[0]].copy();bad['longitude']=4
    with pytest.raises(ValueError,match='Conflicting'):
        module().unique_roster(pd.concat([a,bad]))


def test_scale_offset_preserved_without_turning_nodata_into_a_measurement():
    actual=module().physical_values(np.array([3000,np.nan,0]),{'scale':.1,'offset':-273.15})
    assert np.isclose(actual[0],26.85) and np.isnan(actual[1])


def test_downstream_join_rejects_non_native_or_duplicate_observations():
    m=module();native=m.classify_all(roster(),source());t=pd.DataFrame({'obs_id':['n','i'],'value':[1,2]})
    a=m.native_measurement_join(native,t)
    assert a.obs_id.tolist()==['n']
    with pytest.raises(ValueError,match='unique'):
        m.native_measurement_join(native,pd.concat([t,t.iloc[[0]]]))
