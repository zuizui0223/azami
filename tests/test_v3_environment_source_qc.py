import json
import pandas as pd
import pytest
from analysis.v3.environment_source_qc import CONTRACT, apply_qc
from analysis.v3.workflow import digest


def setup(tmp_path):
    raw = tmp_path/'raw.csv'
    pd.DataFrame({'obs_id':['a','b','c'],'accepted_key':['x','x','y'],
                  'GSP':[429496729.5,20.,30.], 'BIO12':[100.,6553.5,120.],
                  'rsds_month':[150.,200.,180.]}).to_csv(raw,index=False)
    spec = json.loads(CONTRACT.read_text())
    spec.update(input_matrix_sha256=digest(raw), observations=3, taxa=2)
    for r in spec['rules']: r['expected_rows']=1
    path = tmp_path/'contract.json'
    path.write_text(json.dumps(spec))
    return raw,path


def test_qc_changes_only_flagged_working_values_and_preserves_every_raw_row(tmp_path):
    raw,spec = setup(tmp_path)
    before = digest(raw)
    result = apply_qc(raw,tmp_path/'out',spec)
    data = pd.read_csv(tmp_path/'out/environment_candidate_matrix_qc_private.csv')
    assert result['source_rows_deleted']==0 and digest(raw)==before
    assert data.GSP.isna().tolist()==[True,False,False]
    assert data.BIO12.isna().tolist()==[False,True,False]
    assert data.GSP_source_raw.tolist()==[429496729.5,20.,30.]
    assert data.rsds_month.tolist()==[150.,200.,180.]
    assert result['ecological_fitting_authorized'] is False


def test_qc_requires_exact_input_and_does_not_overwrite(tmp_path):
    raw,spec = setup(tmp_path)
    apply_qc(raw,tmp_path/'out',spec)
    with pytest.raises(ValueError,match='Preserve'):
        apply_qc(raw,tmp_path/'out',spec)
    with raw.open('a') as f: f.write('\n')
    with pytest.raises(ValueError,match='identity'):
        apply_qc(raw,tmp_path/'bad',spec)
    assert not (tmp_path/'bad').exists()


def test_unexpected_anomaly_count_stops_without_any_qc_view(tmp_path):
    raw,spec = setup(tmp_path)
    plan=json.loads(spec.read_text()); plan['rules'][0]['expected_rows']=2
    spec.write_text(json.dumps(plan))
    with pytest.raises(ValueError,match='counts differ'):
        apply_qc(raw,tmp_path/'out',spec)
    assert not (tmp_path/'out').exists()
