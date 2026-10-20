from analysis.v4.native_scope import distribution_status, resolve_exact, inclusion_reason


def test_distribution_is_source_specific_and_unlisted_never_native():
    assert distribution_status({'locationId':'TDWG:FRA'}) == 'native'
    assert distribution_status({'locationId':'TDWG:FRA','establishmentMeans':'INTRODUCED'}) == 'introduced'
    assert distribution_status({'locationId':'TDWG:FRA','occurrenceStatus':'DOUBTFUL'}) == 'uncertain'
    assert distribution_status({'locationId':'TDWG:FRA','occurrenceStatus':'ABSENT'}) == 'uncertain'
    assert distribution_status({'locationId':'TDWG:FRA','establishmentMeans':'UNKNOWN'}) == 'uncertain'
    assert distribution_status({}) == 'uncertain'


def test_names_must_be_exact_with_one_convergent_accepted_target():
    rows=[{'canonicalName':'Cirsium x','key':1,'taxonomicStatus':'ACCEPTED'},
          {'canonicalName':'Cirsium x','key':2,'acceptedKey':1,'taxonomicStatus':'SYNONYM'}]
    assert resolve_exact('Cirsium x',rows)['accepted_key']==1
    assert resolve_exact('Cirsium y',rows)['resolution_status']=='unresolved_taxon'
    rows.append({'canonicalName':'Cirsium x','key':3,'taxonomicStatus':'ACCEPTED'})
    assert resolve_exact('Cirsium x',rows)['resolution_status']=='unresolved_taxon'


def test_native_range_and_wild_status_are_separate_gates():
    assert inclusion_reason('native',True,False,False)=='include'
    assert inclusion_reason('introduced',True,False,False)=='range_introduced'
    assert inclusion_reason('unlisted',True,False,False)=='range_unlisted'
    assert inclusion_reason('native',True,False,True)=='captive'
    assert inclusion_reason('native',True,False,None)=='captive_unknown'
    assert inclusion_reason('native',True,True,False)=='obscured'
    assert inclusion_reason('native',False,False,False)=='coordinates_unusable'


def test_geometry_and_missing_distribution_are_fail_closed():
    import pytest
    pytest.importorskip('shapely')
    import pandas as pd
    from analysis.v4.native_scope import classify_roster
    source={'resolutions':[{'input_name':'Cirsium x','accepted_key':1,'accepted_name':'Cirsium x','resolution_status':'resolved_unique_accepted_key'}],
            'distributions':{'1':[{'locationId':'TDWG:AAA'}]},
            'geometry':{'features':[{'properties':{'LEVEL3_COD':'AAA'},'geometry':{'type':'Polygon','coordinates':[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}}]}}
    f=pd.DataFrame({'obs_id':['a','b','c'],'taxon_name':['Cirsium x']*3,'longitude':[1,0,5],'latitude':[1,1,5],
                    'coordinate_usable_for_environment':[True]*3,'obscured':[False]*3,'captive':[False]*3})
    result=classify_roster(f,source)
    assert result.native_primary.tolist()==[True,False,False]
    assert result.native_range_status.tolist()==['native','unmapped_tdwg','unmapped_tdwg']
    source['distributions']={}
    assert classify_roster(f,source).native_range_status.iloc[0]=='unlisted'
    source['distributions']={'1':[{'locationId':'TDWG:AAA'},{'locationId':'TDWG:AAA','establishmentMeans':'INTRODUCED'}]}
    assert classify_roster(f,source).native_range_status.iloc[0]=='ambiguous'
