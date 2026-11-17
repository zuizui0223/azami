"""Synthetic source and null-calibration tests; NO protected original V4 records."""
import numpy as np
import pandas as pd
import pytest

from analysis.v4 import species_environment_predictive_permutation as t
from analysis.v4 import species_environment_spls2 as pls
from analysis.v4.species_environment_scan import FEATURES


def demo(n_taxa=54,photos=5,signal=True):
    rng=np.random.default_rng(11)
    rows=[]
    for i in range(n_taxa):
        cl=rng.normal(size=22)
        sector=i%3
        for j in range(photos):
            x=cl+rng.normal(0,.04,size=22)
            row={'obs_id':f'{i}-{j}','accepted_key':f't{i:04d}',
                 'longitude':[-85,10,125][sector] + rng.normal(0,.03)}
            row.update(zip(pls.CLIMATE,x))
            effects={'angle':(.3*cl[11] if signal else 0),
                     'lightness':(.2*cl[19] if signal else 0),
                     'chroma':(.3*cl[19] if signal else 0),
                     'aspect':0,'circularity':0,'solidity':0,'width_profile':0}
            for trait,cols in FEATURES.items():
                if trait=='hue':
                    theta=.1+(cl[7]*.3 if signal else 0)+rng.normal(0,.04)
                    row.update({cols[0]:float(np.sin(theta)),cols[1]:float(np.cos(theta))})
                else:
                    row[cols[0]]=float(effects[trait]+rng.normal(0,.13))
                for k in cols:row[k+'__eligible']=True
            rows.append(row)
    return pd.DataFrame(rows)


def test_holm_monotone():
    r=t.holm({'a':.01,'b':.02,'c':.6})
    assert r=={'a':.03,'b':.04,'c':.6}


def test_holm_ties():
    r=t.holm({'a':.2,'b':.2,'c':.01})
    assert r['a']==r['b']==.4 and r['c']==.03


def test_shuffled_multitrait_sectors_and_whole_vectors():
    Y=np.arange(180,dtype=float).reshape(60,3)
    sector=np.repeat([0,1,2],20)
    Z=t.shuffled_labels(Y,sector,19,within_sector=True)
    assert sorted(map(tuple,Z))==sorted(map(tuple,Y))
    for s in range(3):
        assert sorted(map(tuple,Y[sector==s]))==sorted(map(tuple,Z[sector==s]))


def test_shuffled_all_taxa_preserves_whole_response_vectors():
    y=np.arange(180.).reshape(60,3);s=np.repeat([0,1,2],20)
    z=t.shuffled_labels(y,s,2,within_sector=False)
    assert sorted(map(tuple,z))==sorted(map(tuple,y))


def test_regional_shuffle_requires_three_sectors():
    with pytest.raises(ValueError,match='THREE_SECTORS'):
        t.shuffled_labels(np.zeros((20,2)),np.zeros(20,int),1,within_sector=True)


def test_rf_uses_original_22_without_elevation():
    source=demo()
    scores=t.prep_rf_species(source)
    assert len(scores)==3
    for x in scores.values():
        assert x['taxa']==54 and x['x'].shape==(54,22)
        assert 'topo_elevation' not in x['features']
    assert len(scores['precipitation__angle']['group'])==8


def test_rf_permutation_importance_not_signed_beta():
    a=t.prep_rf_species(demo())['precipitation__angle']
    s=t.rf_group_stat(a,a['y'],trees=15,n_folds=3,permutations_per_fold=2)
    assert np.isfinite(s['group_importance_normalized'])
    assert np.isfinite(s['pooled_skill'])
    assert len(s['folds'])==3


def test_wrong_response_grain_rejected():
    a=t.prep_rf_species(demo())['precipitation__angle']
    with pytest.raises(ValueError,match='SHAPE_MISMATCH'):
        t.rf_group_stat(a,np.ones((3,1)),trees=10)


def test_spls_stats_are_component_sign_invariant():
    rng=np.random.default_rng(19)
    X=rng.normal(size=(60,22));Y=rng.normal(size=(60,9))
    m=pls.SparsePLS2(22,2).fit(X,Y)
    a=t.pls_stats(m)
    m.weights[:,0]*=-1
    assert t.pls_stats(m)==a  # test derives from beta, not component sign


def test_spls_fitted_stats_and_nested_cv():
    rng=np.random.default_rng(77)
    X=rng.normal(size=(54,22));Y=rng.normal(size=(54,9))
    r=t.pls_cv_and_coeffs(X,Y,nested=True,k_grid=(2,22),h_grid=(1,2))
    assert r['keepX'] in (2,22)
    assert len(r['coefficient_stats']['BIO08__hue_vector'])==2
    assert np.isfinite(r['nested_cv_skill'])


def test_full_pipeline_null_small_produces_three_adjusted_p_values(monkeypatch):
    # Keep full RF fitting and output identical, abbreviate repeated nested CV
    # folds exclusively to stay within fast synthetic test budget.
    rng=np.random.default_rng(42)
    source=demo(54,5)
    original_pls=t.pls_cv_and_coeffs
    def short_pls(X,Y,*,seed=t.SEED,nested=True,k_grid=None,h_grid=None):
        return original_pls(X,Y,seed=seed,nested=nested,k_grid=(2,22),h_grid=(1,))
    monkeypatch.setattr(t,'pls_cv_and_coeffs',short_pls)
    # Use a small null for the test; the production command requires >=19.
    result=t.observed_and_null(source,unrestricted=19,regional=19,trees=10)
    assert result['data']['original_source_observation_n']==270
    assert result['data']['sPLS2_joint_taxa']==54
    for method in result['permutation_inference'].values():
        assert 0<method['sPLS2_global_nested_cv']['null_p_one_sided']<=1
        for item in ('RF_group_importance','sPLS2_selected_coefficient'):
            assert set(method[item])==set(t.CANDIDATES)
            for v in method[item].values():
                assert 0<v['null_p_one_sided']<=v['posthoc_family3_holm']<=1
