import numpy as np
import pandas as pd

from analysis.v4.environment_candidate_audit import CANDIDATES, audit, matrix_summary
from analysis.v4.observation_models import MODULES


def test_singular_predictors_are_reported_not_silently_deleted():
    x=np.arange(30,dtype=float)
    result=matrix_summary(pd.DataFrame({'a':x,'b':2*x}),['a','b'])
    assert result['status']=='singular'
    assert result['columns']==['a','b']


def test_audit_is_invariant_to_trait_values_and_tracks_new_missingness():
    rng=np.random.default_rng(391)
    d=pd.DataFrame({'accepted_key':np.repeat(np.arange(20),10)})
    for v in CANDIDATES:
        d[v]=rng.normal(size=len(d))
    for outcomes in MODULES.values():
        for name in outcomes:
            d[name]=rng.normal(size=len(d))
            d[name+'__eligible']=True
    d.loc[0,CANDIDATES[-1]]=np.nan
    before=audit(d)
    for outcomes in MODULES.values():
        for name in outcomes:
            d[name]=rng.normal(size=len(d))*100
    assert audit(d)==before
    entry=before['populations']['orientation_image_vertical_angle']
    assert entry['four_environment_n']==200
    assert entry['six_environment_n']==199
    assert entry['lost_from_four']==1
