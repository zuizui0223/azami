import pandas as pd
import pytest
import hashlib
import io
import json
import zipfile
from analysis.v4 import measurement_coverage_audit as audit
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


def test_clean_runner_restores_archive_before_auditing(tmp_path,monkeypatch):
    metrics=('involucre_projection_p95','involucre_spread_fraction')
    heads=pd.DataFrame({'obs_id':['a'],'annotation_unit_id':['head1'],
                       **{m+'_status':['usable'] for m in metrics}})
    payload=io.BytesIO()
    with zipfile.ZipFile(payload,'w') as z:
        z.writestr('extended_cohort_audit.json',json.dumps({'n_selected_heads':1}))
        z.writestr('rebuilt/continuous_trait_rebuild_report.json',
                   json.dumps({'n_rebuilt_heads':1,'n_failed_heads':0}))
        z.writestr('extended_measurements/extended_continuous_provenance.json',
                   json.dumps({'n_heads':1,'endpoint_usable_heads':dict.fromkeys(metrics,1)}))
        z.writestr('extended_measurements/extended_continuous_head_level.csv',heads.to_csv(index=False))
        z.writestr('git_sha.txt','fixture')
    raw=payload.getvalue()
    monkeypatch.setattr(audit,'INPUTS',{'continuous':(1,hashlib.sha256(raw).hexdigest(),{})})
    def restore(archives,inputs):
        assert archives.is_dir() and inputs.is_dir()
        (archives/'artifact-1-continuous.zip').write_bytes(raw)
    monkeypatch.setattr(audit,'restore_zenodo',restore)
    native=tmp_path/'native.csv'
    native.write_text('obs_id\na\nb\n')
    result=audit.execute(native,tmp_path/'fresh',tmp_path/'out/audit.json')
    assert result['native_population_observations']==2
    assert result['full_native_measurement_completed'] is False
    assert result['endpoint_dispositions'][metrics[0]]['current_native']['usable_observations']==1
