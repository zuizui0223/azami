"""Synthetic only; genuine V4 requires protected source and SHA receipts."""
import json
import numpy as np
import pandas as pd
import pytest
from analysis.v4 import species_environment_partial_rda as rd
from analysis.v4.species_environment_spls2 import CLIMATE, YCOLS


def fixture(n_taxa=85, photos=6, seed=27, signal='climate'):
    rng=np.random.default_rng(seed)
    taxa=[f'synthetic_{i:03}' for i in range(n_taxa)]
    base=rng.normal(size=(n_taxa,22)); lat=rng.uniform(-52,64,n_taxa)
    lon=np.r_[rng.uniform(-169,-45,n_taxa//3),rng.uniform(-23,50,n_taxa//3),rng.uniform(65,170,n_taxa-2*(n_taxa//3))]
    shift=(lon-20)/170
    base[:,0]=shift+rng.normal(size=n_taxa)*.35
    base[:,1]=base[:,0]*.8 + rng.normal(size=n_taxa)*.24
    responses=rng.normal(size=(n_taxa,9))*.2
    if signal=='climate':responses[:,0]=base[:,4]*2.5+rng.normal(size=n_taxa)*.15
    elif signal=='space':responses[:,0]=np.cos(np.deg2rad(lon))*3+rng.normal(size=n_taxa)*.15
    if signal=='climate':responses[:,1]=base[:,5]*1.4+rng.normal(size=n_taxa)*.25
    responses[:,3:5] = rng.normal(size=(n_taxa,2))
    responses[:,3:5] /= np.linalg.norm(responses[:,3:5],axis=1)[:,None]
    out=[]
    for i,tax in enumerate(taxa):
        for j in range(photos):
            row={'obs_id':f'fake_{i:03}_{j}','accepted_key':tax,'longitude':float(lon[i]+rng.normal()*0.01),'latitude':float(lat[i]+rng.normal()*0.01)}
            for k,p in enumerate(CLIMATE):row[p]=float(base[i,k]+rng.normal()*.025)
            for k,p in enumerate(YCOLS):row[p]=float(responses[i,k]+rng.normal()*.02);row[p+'__eligible']=True
            out.append(row)
    return pd.DataFrame(out)


def test_coordinate_basis_invariant_to_dateline():
    a=rd.spherical_geo(np.array([32.]),np.array([-180.]))
    b=rd.spherical_geo(np.array([32.]),np.array([180.]))
    assert a.shape==(1,8)
    np.testing.assert_allclose(a,b,rtol=0,atol=1e-12)


def test_joint_aggregation_matches_spls2_cohort_and_no_8_cells_required():
    d=fixture(n_taxa=84)
    taxa,audit=rd.aggregate_joint_geo(d)
    assert len(taxa)==84
    assert audit['n_joint_supported_taxa']==84
    assert audit['source_photo_contributions_to_taxon_means']==84*6
    assert set(CLIMATE).issubset(taxa.columns)
    assert all(p in taxa for p in rd.GEO_COLS)
    assert len(taxa._taxon_private.unique())==84


def test_joint_hue_variance_is_one_domain_not_two():
    taxa,_=rd.aggregate_joint_geo(fixture())
    y=rd.make_response_matrix(taxa,list(YCOLS))
    h=[list(YCOLS).index(x) for x in ('corolla_hue_sin','corolla_hue_cos')]
    assert np.isclose(y[:,h].var(axis=0,ddof=1).sum(),1,atol=1e-12)
    assert np.isclose(y[:,0].var(ddof=1),1,atol=1e-12)


def test_species_mean_outcome_blind_pc_selection_and_not_22_regressions():
    taxa,_=rd.aggregate_joint_geo(fixture())
    x,diag=rd.climate_pcs(taxa,6)
    assert x.shape==(85,6)
    assert len(diag['variance_ratio'])==6
    assert 0<sum(diag['variance_ratio'])<=1.0000000001
    assert diag['learned_from_traits'] is False


def test_climate_signal_stays_after_space_control():
    d=fixture(signal='climate')
    out=rd.analyze(d,n_permutations=29,seed=20261107)
    row=out['families']['all_eight']['main_6_climate_pcs']
    assert row['n_taxa']==85
    assert row['adjR2_unique_climate']>0.01
    assert row['rawR2_climate_only']>0
    assert row['exploratory_p_climate_beyond_geo']<=1
    assert row['n_geo_dimensions']==8
    assert out['historical_phylogeny_status']=='not_estimated_no_authenticated_V4_matched_tree'


def test_space_signal_is_captured_by_sphere_and_response_composition():
    d=fixture(signal='space')
    out=rd.analyze(d,n_permutations=19,seed=2026)
    row=out['families']['all_eight']['main_6_climate_pcs']
    assert row['rawR2_geography_only']>0.01
    assert np.isfinite(row['adjR2_unique_geography'])


def test_reordered_original_observation_rows_dont_change_result():
    d=fixture(n_taxa=80)
    a=rd.analyze(d,n_permutations=19,seed=11)
    b=rd.analyze(d.sample(frac=1,random_state=2026),n_permutations=19,seed=11)
    assert a==b


def test_identical_taxa_for_main_and_sensitivity_and_all_families():
    d=fixture(n_taxa=83)
    r=rd.analyze(d,n_permutations=9)
    assert set(r['families'])=={'all_eight','color','angle','outline'}
    assert set(r['families']['all_eight'])=={'main_6_climate_pcs','sensitivity_4_climate_pcs','sensitivity_8_climate_pcs','sensitivity_22_original_climate_axes'}
    for fam in r['families'].values():
        for val in fam.values():assert val['n_taxa']==83


def test_null_signal_does_not_get_confirmed_status():
    r=rd.analyze(fixture(signal='noise'),n_permutations=19)
    assert r['independent_replication'] is False
    assert r['causal_environment_effects_identified'] is False
    assert r['multivariate_permutation_assumptions_are_verified'] is False


def test_missing_coordinate_and_unmatched_source_ids_fail_closed():
    d=fixture()
    d.loc[0,'longitude']=np.nan
    with pytest.raises(ValueError,match='COORDINATE'):rd.aggregate_joint_geo(d)
    d=fixture();d.loc[0,'obs_id']=d.obs_id.iloc[1]
    with pytest.raises(ValueError,match='IDENTIT'):rd.aggregate_joint_geo(d)


def test_single_endpoint_missing_with_qc_excludes_only_same_photo_not_fake_value():
    d=fixture(n_taxa=85)
    d.loc[0,'corolla_lab_chroma__eligible']=False
    d.loc[0,'corolla_lab_chroma']=-999999.
    taxa,audit=rd.aggregate_joint_geo(d)
    assert len(taxa)==85
    assert audit['all_eight_joint_photo_eligible']==len(d)-1


def test_diagnostics_and_output_are_privacy_safe():
    out=rd.analyze(fixture(),n_permutations=9)
    blob=json.dumps(out)
    for s in ('synthetic_000','fake_000_0','_taxon_private','taxon_name','longitude','latitude','obs_id','photo_id'):
        assert s not in blob


def test_adjusted_partition_matches_direct_linear_algebra():
    rng=np.random.default_rng(1)
    Y=rng.normal(size=(100,4));C=rng.normal(size=(100,4));G=rng.normal(size=(100,3))
    row=rd.varpart(Y,C,G,n_permutations=9,seed=15,sector=None)
    def adj(x):
        x=np.column_stack([np.ones(len(x)),x]); B=np.linalg.lstsq(x,Y,rcond=None)[0];res=Y-x@B
        ss=np.square(Y-Y.mean(axis=0)).sum();n=len(x);p=x.shape[1]-1
        return 1-(np.square(res).sum()/ss)*(n-1)/(n-p-1)
    assert np.isclose(row['adjR2_unique_climate'], adj(np.column_stack([C,G]))-adj(G))
    assert np.isclose(row['adjR2_unique_geography'], adj(np.column_stack([C,G]))-adj(C))


def test_invalid_or_missing_env_species_is_not_silently_imputed():
    d=fixture(n_taxa=85)
    d.loc[d.accepted_key=='synthetic_000','chelsa_bio08']=np.nan
    result=rd.analyze(d,n_permutations=9)
    assert result['n_joint_supported_taxa']==84
    assert result['n_joint_photo_supported_before_all22_climate_complete_case']==85
