import pandas as pd
import pytest
from analysis.v4.measurement_coverage_audit import summarize_heads


def test_attrition_preserves_overlapping_failures_and_native_population():
    d=pd.DataFrame({'obs_id':['a','a','b','c'],
                    'annotation_unit_id':['1','2','3','4']})
    for metric in ('involucre_projection_p95','involucre_spread_fraction'):
        d[metric+'_status']=['usable','low_resolution;measurement_missing',
                            'low_resolution','usable']
    a=summarize_heads(d,{'a','b'})['involucre_projection_p95']
    assert a['all_measured']['usable_observations']==2
    assert a['current_native']['usable_observations']==1
    assert a['current_native']['heads']==3
    assert a['current_native']['unusable_heads']==2
    assert a['current_native']['overlapping_failure_reasons']=={
        'low_resolution':2,'measurement_missing':1}
    assert sum(a['current_native']['exclusive_status_combinations'].values())==3
    d.loc[0,'involucre_projection_p95_status']=None
    with pytest.raises(ValueError,match='Unrecorded'):
        summarize_heads(d,{'a','b'})
