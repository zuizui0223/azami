from importlib.util import find_spec
import json
import numpy as np
import pandas as pd
import pytest

MOD='analysis.v4.species_environment_scan'

def module():
    assert find_spec(MOD) is not None, 'species-level scan has not been implemented'
    from analysis.v4 import species_environment_scan as m
    return m

def observations():
    rng=np.random.default_rng(13); rows=[]
    for t in range(25):
        for k in range(6):
            x=(t-12)/6+rng.normal(0,.1)
            rows.append(dict(obs_id=f'o{t}_{k}', accepted_key=f't{t:02}',
                latitude=35.,longitude=140.,chelsa_bio01=x,chelsa_bio12=100*x+1500,
                orientation_image_vertical_angle=50+4*x+rng.normal(0,1),
                orientation_image_vertical_angle__eligible=True))
    return pd.DataFrame(rows)

def test_registry_all23_and_eight_not_seven():
    m=module(); assert len(m.PREDICTORS)==23 and len(m.FEATURES)==8
    assert {'chelsa_bio18','chelsa_bio05','topo_elevation'} <= set(m.PREDICTORS)

def test_species_with_only_one_geographic_cell_and_six_photos_are_retained():
    m=module(); p=m.aggregate_pair(observations(),'angle','chelsa_bio01')
    assert len(p)==25 and p.n_observations.min()==6

def test_joint_missingness_uses_same_rows_for_x_and_y_means():
    m=module(); d=observations(); mask=d.accepted_key.eq('t00')
    i=d.index[mask][0];d.loc[i,'chelsa_bio01']=np.nan
    d.loc[i,'orientation_image_vertical_angle']=1e8
    p=m.aggregate_pair(d,'angle','chelsa_bio01').set_index('taxon')
    assert p.loc['t00','n_observations']==5 and p.loc['t00','y0']<100

def test_unrelated_environment_missingness_does_not_drop_taxa():
    m=module();d=observations();d['chelsa_bio19']=np.nan
    assert len(m.aggregate_pair(d,'angle','chelsa_bio01'))==25

def test_photo_qc_is_not_ignored():
    m=module();d=observations();d.loc[d.accepted_key.eq('t00'),'orientation_image_vertical_angle__eligible']=False
    assert len(m.aggregate_pair(d,'angle','chelsa_bio01'))==24

def test_duplicate_id_refused():
    m=module();d=observations();d=pd.concat([d,d.iloc[[0]]])
    with pytest.raises(ValueError,match='IDENTITY'): m.validate_observations(d)

def test_missing_taxon_refused_as_input_identity_error():
    m=module();d=observations();d.loc[0,'accepted_key']=None
    with pytest.raises(ValueError,match='IDENTITY'):m.validate_observations(d)

def test_scalar_hc3_agrees_with_independent_statsmodels():
    m=module();p=m.aggregate_pair(observations(),'angle','chelsa_bio01')
    from statsmodels.api import OLS,add_constant
    ref=OLS(p.y0,add_constant(p.x)).fit(cov_type='HC3',use_t=True)
    out=m.fit_pair(p,'angle')
    assert out['beta_raw']==pytest.approx(ref.params.iloc[1],rel=1e-10)
    assert out['p_value']==pytest.approx(ref.pvalues.iloc[1],rel=1e-8)

def test_constant_predictor_is_not_significant_or_fake_zero():
    m=module();p=m.aggregate_pair(observations(),'angle','chelsa_bio01');p['x']=1.
    out=m.fit_pair(p,'angle');assert out['status']=='constant_predictor' and out['p_value'] is None

def test_constant_response_not_called_zero_effect():
    m=module();p=m.aggregate_pair(observations(),'angle','chelsa_bio01');p['y0']=3.
    out=m.fit_pair(p,'angle');assert out['status']=='constant_response' and out['p_value'] is None

def test_small_panel_not_fit():
    m=module();p=m.aggregate_pair(observations(),'angle','chelsa_bio01').iloc[:10]
    assert m.fit_pair(p,'angle')['status']=='insufficient_taxa'

def test_joint_hue_rotation_invariant_wald():
    m=module();rng=np.random.default_rng(3);n=80;x=np.linspace(-2,2,n)
    th=.4*x+rng.normal(0,.2,n)
    p=pd.DataFrame(dict(taxon=[f't{i}' for i in range(n)],x=x,y0=np.sin(th),y1=np.cos(th),n_observations=6))
    a=m.fit_pair(p,'hue');q=p.copy();v=p[['y0','y1']].to_numpy();rot=np.array([[np.cos(.8),-np.sin(.8)],[np.sin(.8),np.cos(.8)]])
    q[['y0','y1']]=v@rot;b=m.fit_pair(q,'hue')
    assert a['p_value']==pytest.approx(b['p_value'],rel=1e-7)
    assert a['beta_raw'] is None and a['joint_dimensions']==2

def test_hue_mean_wraparound_and_ambiguous_mean():
    m=module();d=observations();a=np.tile([179,-179,179,-179,179,-179],25)
    d['corolla_hue_sin']=np.sin(np.deg2rad(a));d['corolla_hue_cos']=np.cos(np.deg2rad(a))
    d['corolla_hue_sin__eligible']=True;d['corolla_hue_cos__eligible']=True
    p=m.aggregate_pair(d,'hue','chelsa_bio01');assert (p.y1<-.99).all()
    mask=d.accepted_key.eq('t00');d.loc[mask,'corolla_hue_sin']=[0]*6;d.loc[mask,'corolla_hue_cos']=[1,-1,1,-1,1,-1]
    assert len(m.aggregate_pair(d,'hue','chelsa_bio01'))==24

def test_covariance_alignment_and_gls_match_statsmodels():
    m=module();p=m.aggregate_pair(observations(),'angle','chelsa_bio01');n=len(p)
    C=.2*np.ones((n,n))+.8*np.eye(n)
    obj=m.covariance_from_payload(dict(kind='phylogeny_and_space_total_residual',taxa=p.taxon.tolist()[::-1],matrix=C[::-1,::-1].tolist(),source_hashes={'phylogeny':'a'*64,'space':'b'*64},parameter_scope='externally_fixed'))
    V=obj.subset(p.taxon.tolist());out=m.fit_pair(p,'angle',V)
    from statsmodels.api import GLS,add_constant
    ref=GLS(p.y0,add_constant(p.x),sigma=C).fit()
    assert out['beta_raw']==pytest.approx(ref.params.iloc[1],rel=1e-10)
    assert out['p_value']==pytest.approx(ref.pvalues.iloc[1],rel=1e-7)

def test_covariance_without_phylogeny_or_space_not_accepted():
    m=module()
    with pytest.raises(ValueError,match='COVARIANCE'):
        m.covariance_from_payload(dict(kind='identity',taxa=['a','b'],matrix=np.eye(2).tolist()))

def test_nonpositive_covariance_is_not_silently_repaired():
    m=module()
    with pytest.raises(ValueError,match='COVARIANCE'):
        m.covariance_from_payload(dict(kind='phylogeny_and_space_total_residual',taxa=['a','b'],matrix=[[1,2],[2,1]],source_hashes={'phylogeny':'a'*64,'space':'b'*64},parameter_scope='externally_fixed'))

def test_fdr_retains_failed_tests_in_denominator():
    m=module();r=[{'p_value':.01},{'p_value':None},{'p_value':.03},{'p_value':None}]
    m.adjust_family(r);assert r[0]['q_bh']==pytest.approx(.04)
    assert r[1]['q_bh'] is None and r[0]['multiplicity_family_size']==4

def test_all_184_rows_in_each_panel_even_when_missing_inputs():
    m=module();out=m.scan(observations());assert len(out['results'])==3*184
    assert sum(r['panel']=='marginal' and r['status']=='estimated' for r in out['results'])==2
    assert all(r['p_value'] is None for r in out['results'] if r['panel']!='marginal')
    text=json.dumps(out,allow_nan=False)
    for private in ('"taxon"','"obs_id"','o0_0','t00','"latitude"','"longitude"'):assert private not in text

def test_covariance_subset_has_matched_marginal_and_explicit_smaller_n():
    m=module();d=observations();labels=sorted(d.accepted_key.unique())[:22]
    co=m.covariance_from_payload(dict(kind='phylogeny_and_space_total_residual',taxa=labels,matrix=np.eye(22).tolist(),source_hashes={'phylogeny':'a'*64,'space':'b'*64},parameter_scope='externally_fixed'))
    rows=[r for r in m.scan(d,covariance=co)['results'] if r['construct']=='angle' and r['predictor']=='chelsa_bio01']
    assert [r['n_taxa'] for r in rows]==[25,22,22]

def test_cluster_input_does_not_access_trait_values():
    m=module();d=observations();a=m.environment_clusters(d)
    d['orientation_image_vertical_angle']=np.arange(len(d))*100
    assert a==m.environment_clusters(d)
    assert a['cluster_by_predictor']['chelsa_bio01']==a['cluster_by_predictor']['chelsa_bio12']

def test_run_requires_actual_checksum_and_keeps_execution_status(tmp_path):
    m=module();source=tmp_path/'table.csv';observations().to_csv(source,index=False)
    receipt=tmp_path/'receipt.json';meta={'schema':'azami_authorized_joined_s1_v1','table_sha256':m.sha256(source),'source_row_count':150,'environment_values_in_physical_units':True,'source_lineage_sha256':'a'*64}
    receipt.write_text(json.dumps(meta));output=tmp_path/'out.json'
    out=m.run(source,receipt,output,expected_n=150)
    assert out['source_row_count']==150 and out['independent_confirmation'] is False
    meta['table_sha256']='0'*64;receipt.write_text(json.dumps(meta))
    with pytest.raises(ValueError,match='SOURCE'):m.run(source,receipt,output,expected_n=150)

def test_input_taxon_count_is_separate_from_pairwise_eligible_count():
    m=module();d=observations();d.loc[d.accepted_key.eq('t00'),'orientation_image_vertical_angle__eligible']=False
    out=m.scan(d)
    assert out['n_input_taxa']==25
    r=next(r for r in out['results'] if r['construct']=='angle' and r['predictor']=='chelsa_bio01' and r['panel']=='marginal')
    assert r['n_taxa']==24

def test_all244_single_cell_taxa_can_enter_without_rf_geographic_filter():
    m=module();rng=np.random.default_rng(8);d=observations().iloc[:0].copy();records=[]
    for i in range(244):
        for j in range(5):
            records.append({'obs_id':f'i{i}_{j}','accepted_key':f't{i}',
                'chelsa_bio01':float(i%21),'orientation_image_vertical_angle':float(30+.1*i+rng.normal()),
                'orientation_image_vertical_angle__eligible':True,'latitude':35.,'longitude':140.})
    out=m.aggregate_pair(pd.DataFrame(records),'angle','chelsa_bio01')
    assert len(out)==244

def test_gls_joint_hue_rotation_invariant():
    m=module();rng=np.random.default_rng(123);n=40;x=rng.normal(size=n);th=.4*x+rng.normal(0,.8,n)
    p=pd.DataFrame(dict(taxon=[str(i) for i in range(n)],x=x,y0=np.sin(th),y1=np.cos(th),n_observations=5))
    C=.1*np.ones((n,n))+.9*np.eye(n)
    a=m.fit_pair(p,'hue',C);r=p.copy();r['y0']=p.y1;r['y1']=-p.y0
    b=m.fit_pair(r,'hue',C)
    assert a['p_value']==pytest.approx(b['p_value'],rel=1e-9)

def test_row_order_invariance():
    m=module();d=observations()
    assert m.scan(d)==m.scan(d.sample(frac=1,random_state=9))

def test_support_counts_survive_insufficient_taxon_fit():
    m=module();p=m.aggregate_pair(observations(),'angle','chelsa_bio01').iloc[:10]
    r=m.fit_pair(p,'angle')
    assert r['n_photo_contributions']==60
    assert r['observations_per_taxon_min_median_max']==[6,6.,6]
    assert r['status']=='insufficient_taxa' and r['p_value'] is None


def test_empty_pair_reports_zero_not_unknown_photo_support():
    m=module();p=m.aggregate_pair(observations(),'chroma','chelsa_bio01')
    r=m.fit_pair(p,'chroma')
    assert r['n_photo_contributions']==0
    assert r['observations_per_taxon_min_median_max'] is None
