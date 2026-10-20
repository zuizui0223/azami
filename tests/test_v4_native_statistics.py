import numpy as np
import pandas as pd
from analysis.v4.native_statistics import moments, integration
from analysis.v4.native_environment import partition, ridge_predict, fold_ids


def toy():
    rng=np.random.default_rng(8);rows=[]
    for g in range(8):
        for j in range(10+g):
            v=rng.normal(size=5)+np.array([g,g,-g,0,0])
            rows.append([str(g),str(j%3),*v])
    return pd.DataFrame(rows,columns=['taxon_name','cell','angle','light','chroma','hs','hc'])


def test_total_covariance_law_and_equal_taxon_weight():
    f=toy();m=moments(f,['angle','light','chroma','hs','hc'])
    assert m['n_taxa']==8
    w=m['weights'];y=m['Y'];mu=np.average(y,axis=0,weights=w)
    expected=((y-mu)*w[:,None]).T@(y-mu)
    assert np.allclose(m['within']+m['among'],expected)
    assert np.allclose(pd.Series(w).groupby(f.taxon_name).sum(),1/8)


def test_joint_hue_is_rotation_invariant_and_repeatable():
    f=toy();cols=['angle','light','chroma','hs','hc']
    a=integration(f,cols,replicates=50,seed=4,colour_metric='joint_hue')
    theta=.71;rot=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    f[['hs','hc']]=f[['hs','hc']].to_numpy()@rot
    b=integration(f,cols,replicates=50,seed=4,colour_metric='joint_hue')
    assert np.isclose(a['within_rv'],b['within_rv'])
    assert np.isclose(a['among_rv'],b['among_rv'])
    assert np.allclose(a['nested_delta_interval'],b['nested_delta_interval'])
    assert a['interval_definition'].startswith('taxon bootstrap')


def test_unique_rows_and_minimum_replication():
    f=toy();f=f[f.taxon_name!='0']
    a=integration(f,['angle','light','chroma','hs','hc'],replicates=20)
    assert a['n_taxa']==7
    assert 0<=a['nested_positive_fraction']<=1


def test_partition_adds_back_and_exact_linear_signal():
    rng=np.random.default_rng(3);e=rng.normal(size=(60,2));s=rng.normal(size=(60,3));y=e[:,:1]+s[:,:1]
    r=partition(e,s,y)
    assert np.isclose(r['pure_environment']+r['pure_space']+r['shared']+r['unexplained'],1)
    assert r['combined']>.99999


def test_prediction_uses_training_only_scaling():
    x=np.arange(20.)[:,None];y=2*x+3;test=np.array([[25.],[100.]])
    p=ridge_predict(x,y,test,alpha=0)
    assert np.allclose(p,2*test+3)
    assert np.isclose(p[0],ridge_predict(x,y,test[:1],alpha=0)[0]).all()


def test_geographic_groups_not_split():
    g=np.repeat(np.arange(12),np.arange(1,13));ids=fold_ids(g)
    assert len(set(ids))==5
    for a in np.unique(g):assert len(set(ids[g==a]))==1
