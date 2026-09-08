import numpy as np
import pandas as pd
import pytest

from analysis.v3 import source_environment_identification as audit
from analysis.v3.workflow import digest


def test_source_audit_does_not_read_traits_or_reselect_missing_predictors(tmp_path,monkeypatch):
    rng=np.random.default_rng(18087)
    n=80
    variables=['pr_month','rsds_month','vpd_month','sfcWind_month']
    data=pd.DataFrame(rng.normal(size=(n,4)),columns=variables)
    data['obs_id']=[str(i) for i in range(n)]
    data['accepted_key']=np.repeat(['a','b'],40)
    data['native_range_status']='native'
    data['analysis_latitude']=rng.uniform(10,70,n)
    data['analysis_longitude']=rng.uniform(-170,170,n)
    data['observed_year']=np.r_[1950,rng.integers(2000,2026,n-1)]
    data['sin_doy']=np.sin(np.arange(n))
    data['cos_doy']=np.cos(np.arange(n))
    for key in ('south_indicator','south_sin','south_cos'): data[key]=0
    data['trait_must_not_be_read']='not a numeric trait'
    data.loc[4,'pr_month']=np.nan
    path=tmp_path/'source.csv'; data.to_csv(path,index=False)
    monkeypatch.setattr(audit,'definition',lambda:{'source_qc_matrix_sha256':digest(path),'variables':variables})
    report=audit.run(path,tmp_path/'out')
    assert report['source_rows']==80
    assert report['excluded_missing_required_exposure_or_calendar']==1
    assert report['diagnostics']['observations']==79
    assert not report['source_variables_reselected'] and report['trait_values_read']==0
    membership=pd.read_csv(tmp_path/'out/source_membership_private.csv')
    assert 0 in set(membership.obs_id) and 4 not in set(membership.obs_id)
    assert len(data)==80 and digest(path)==report['source_matrix_sha256']
    with pytest.raises(ValueError,match='Preserve prior'):
        audit.run(path,tmp_path/'out')
