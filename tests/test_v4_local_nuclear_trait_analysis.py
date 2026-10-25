from io import StringIO
import numpy as np
import pandas as pd
import pytest
from Bio import Phylo
from analysis.v4 import local_nuclear_trait_analysis as m

def sample_inputs():
    manifest=pd.DataFrame({'run':['r1','r2'],'accepted_key':['1','2']})
    moments=pd.DataFrame([dict(cohort='common8',accepted_key=k,coordinate=c,n_observations=5,mean=v,sample_sd=1.) for k,v in [('1',2.),('2',4.)] for c in ['x','env']])
    return moments,manifest

def test_join_uses_accepted_key_not_row_order_and_preserves_counts():
    moments,manifest=sample_inputs()
    result=m.join_moments(moments.iloc[::-1],manifest.iloc[::-1],'common8',['x','env'])
    assert result.loc['r1','x']==2 and result.loc['r2','env']==4
    assert result.n_original_observations.to_dict()=={'r2':5,'r1':5}

def test_join_rejects_different_photo_cohorts_duplicates_and_missing_values():
    moments,manifest=sample_inputs()
    bad=moments.copy();bad.loc[0,'n_observations']=6
    with pytest.raises(ValueError):m.join_moments(bad,manifest,'common8',['x','env'])
    with pytest.raises(ValueError):m.join_moments(pd.concat([moments,moments.iloc[:1]]),manifest,'common8',['x','env'])
    bad=moments.copy();bad.loc[0,'mean']=np.nan
    with pytest.raises(ValueError):m.join_moments(bad,manifest,'common8',['x','env'])

def test_outgroup_rooting_and_tip_removal_preserve_pairwise_distances():
    tree=Phylo.read(StringIO('((a:1,b:2):3,(c:2,d:3):4,OUTGROUP_saff:5);'),'newick')
    result=m.rooted_ingroup(tree,['a','b','c','d'],['a','b','c'])
    assert {t.name for t in result.get_terminals()}=={'a','b','c'}
    assert {t.name for t in tree.get_terminals()}=={'a','b','c','d','OUTGROUP_saff'}
    for x in ['a','b','c']:
        for y in ['a','b','c']:assert result.distance(x,y)==pytest.approx(tree.distance(x,y))
    with pytest.raises(ValueError):m.rooted_ingroup(tree,['a','b','c'],['a','b','c'])

def test_continuous_unit_returns_analytic_ancestor_and_paired_hue_test():
    tree=Phylo.read(StringIO('('+','.join(f't{i}:1' for i in range(8))+');'),'newick')
    x=np.arange(8,dtype=float)
    frame=pd.DataFrame({'x':x,'sin':np.sin(x),'cos':np.cos(x),
        'chelsa_bio18':x**2,'chelsa_bio05':x+1,'chelsa_rsds_mean':np.sqrt(x+1),
        'elevation_m':x[::-1]},index=[f't{i}' for i in range(8)])
    result=m.evaluate_unit(tree,frame,['x'],nperm=999,seed=7)
    assert result['ancestry']['nodes'][0]['estimate']['x']==pytest.approx(3.5)
    assert len(result['scopes']['all_edges']['axes'])==4
    assert result['scopes']['internal_edges']['status']=='INSUFFICIENT_INTERNAL_EDGES'
    hue=m.evaluate_unit(tree,frame,['sin','cos'],nperm=999,seed=7)
    assert all(a['score_type']=='multivariate_nondirectional_norm' for a in hue['scopes']['all_edges']['axes'].values())
