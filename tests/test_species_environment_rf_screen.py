"""All source observations below are synthetic; no protected S1 used."""
import json
import numpy as np
import pandas as pd
import pytest
from analysis.v4 import species_environment_scan as atlas
from analysis.v4 import species_environment_rf_screen as rf


def fake(n_taxa=65, n_photos=7, *, signal=True, geo=False, missing=False):
    rng = np.random.default_rng(20261107)
    rows=[]
    for t in range(n_taxa):
        rain=rng.normal()
        temp=rng.normal()
        for photo in range(n_photos):
            x=rain+rng.normal(scale=.04)
            y=rain*12+rng.normal(scale=1.) if signal else rng.normal(scale=12.)
            record={'obs_id':f's{t}_{photo}','accepted_key':f'secret_taxon_{t}',
                    'orientation_image_vertical_angle':y,
                    'orientation_image_vertical_angle__eligible': True,
                    'latitude': -50+t%3*42+photo*.01,'longitude':[-120,15,115][t%3]+photo*.01}
            for name in atlas.PREDICTORS:
                record[name]=x if name=='chelsa_bio12' else (x+.0001*rng.normal() if name=='chelsa_bio18' else temp+rng.normal())
            if missing and t%5==0:record['chelsa_bio18']=np.nan
            for trait,cols in atlas.FEATURES.items():
                if trait=='angle':continue
                for c in cols:
                    record[c]=(.7+0.1*rng.normal()) if trait!='hue' else (0.3 if c.endswith('_sin') else 0.95)
                    record[c+'__eligible']=True
            rows.append(record)
    return pd.DataFrame(rows)


def test_single_cell_species_are_retained_and_all_eight_traits_summarized():
    f=fake(n_taxa=65)
    d=rf.aggregate_species(f,'angle')
    assert d['taxon'].nunique()==65
    assert int(d.n_source_photos.sum())==455
    assert 'chelsa_bio12' in d.columns
    assert len(rf.TRAITS)==8
    for c in rf.TRAITS:
        assert len(rf.aggregate_species(f,c))==65


def test_response_cohort_is_not_removed_by_missing_one_environment_axis():
    f=fake(missing=True)
    a=rf.aggregate_species(f,'angle')
    assert len(a)==65
    assert a['chelsa_bio18'].isna().sum()==13
    assert a['chelsa_bio12'].notna().all()


def test_train_imputation_does_not_consult_test_taxa():
    xtr=np.array([[1.,np.nan],[3.,2.],[np.nan,4.]])
    xte=np.array([[999.,np.nan]])
    a,b,m=rf.impute_train_only(xtr,xte)
    assert m.tolist()==[2.,3.]
    assert np.allclose(a,[[1,3],[3,2],[2,4]])
    assert np.allclose(b,[[999,3]])


def test_planted_rain_signal_is_importance_candidate_and_crossval_works():
    f=fake(n_taxa=70)
    cells=rf.aggregate_species(f,'angle')
    result=rf.evaluate_trait(cells,'angle',trees=35,repeats=1,permutation_repeats=1,n_splits=5,seed=42,geography=False)
    assert result['status']=='exploratory_fitted_synthetic_or_authorized_inputs'
    assert result['n_taxa'] == 70
    assert result['pooled_skill_vs_training_mean'] > .5
    imp={row['name']:row for row in result['variable_importance']}
    assert 'chelsa_bio12' in imp
    assert imp['chelsa_bio12']['n_folds'] == 5
    fs={r['family']:r for r in result['family_importance']}
    assert fs['precipitation']['mean_skill_loss'] > 0


def test_output_does_not_contain_taxon_ids_coordinates_or_predictions():
    f=fake(n_taxa=60)
    result=rf.evaluate_trait(rf.aggregate_species(f,'angle'),'angle',trees=20,repeats=1,permutation_repeats=1,seed=42,geography=False)
    text=json.dumps(result)
    for word in ('secret_taxon','"taxon"','"obs_id"','"latitude"','"longitude"','predicted_values'):
        assert word not in text


def test_missing_geographic_source_returns_explicit_status():
    f=fake(n_taxa=65).drop(columns=['latitude','longitude'])
    r=rf.evaluate_trait(rf.aggregate_species(f,'angle'),'angle',trees=20,repeats=1,permutation_repeats=1,seed=42,geography=True)
    assert r['geographic_holdout']['status']=='unavailable_no_source_locations'


def test_taxon_hue_is_unit_circle_and_not_scored_as_scalar_r():
    f=fake(n_taxa=60)
    x=rf.aggregate_species(f,'hue')
    assert np.allclose(np.hypot(x.y0,x.y1),1)
    z=rf.evaluate_trait(x,'hue',trees=15,repeats=1,permutation_repeats=1,seed=42,geography=False)
    assert z['response_dimensions']==2
    assert z['loss']=='joint_unit_hue_chord_squared'
    assert z['status'] in ('exploratory_fitted_synthetic_or_authorized_inputs','insufficient_response_variation')


def test_inadequate_covariate_support_is_not_zero_importance():
    f=fake(n_taxa=60)
    f.loc[f.accepted_key.str.contains(r'secret_taxon_(?:[0-4]?\d)$'), 'chelsa_bio04']=np.nan
    x=rf.aggregate_species(f,'angle')
    r=rf.evaluate_trait(x,'angle',trees=15,repeats=1,permutation_repeats=1,seed=42,geography=False)
    statuses={z['name']:z['status'] for z in r['predictor_coverage']}
    assert statuses['chelsa_bio04']!='used'
    assert 'chelsa_bio04' not in {z['name'] for z in r['variable_importance']}


def test_all_outputs_reproducible_after_row_order_shuffle():
    f=fake(n_taxa=60)
    a=rf.evaluate_trait(rf.aggregate_species(f,'angle'),'angle',trees=20,repeats=1,permutation_repeats=1,seed=42,geography=False)
    b=rf.evaluate_trait(rf.aggregate_species(f.sample(frac=1,random_state=33),'angle'),'angle',trees=20,repeats=1,permutation_repeats=1,seed=42,geography=False)
    assert a == b


def test_source_with_missing_required_trait_columns_reports_status():
    f=fake(n_taxa=45).drop(columns=['capitulum_outline_aspect_ratio'])
    a=rf.aggregate_species(f,'aspect')
    assert a.empty
    r=rf.evaluate_trait(a,'aspect',geography=False)
    assert r['status']=='source_trait_columns_unavailable'


def test_variate_observer_qc_policy_not_claimed_by_species_rf():
    f=fake(n_taxa=42)
    r=rf.run_screen(f,trees=10,repeats=1,permutation_repeats=1,geography=False)
    assert len(r['results'])==8
    assert r['new_species_mean_measurement_reliability_estimated'] is False
    assert r['causal_or_confirmatory_discovery'] is False
    assert len(r['registered_predictors']) == 23
    assert all(len(x['predictor_coverage']) == 23 for x in r['results'] if x['status'].startswith('exploratory'))


def test_geographic_holdout_only_reports_broad_three_sector_statistics():
    f=fake(n_taxa=66)
    d=rf.aggregate_species(f,'angle')
    r=rf.evaluate_trait(d,'angle',trees=15,repeats=1,permutation_repeats=1,geography=True,seed=17)
    assert r['geographic_holdout']['status']=='estimated_exploratory_three_sectors'
    assert len(r['geographic_holdout']['sectors'])==3
    assert all(z['heldout_taxa']>=10 for z in r['geographic_holdout']['sectors'])
    assert 'secret_taxon' not in json.dumps(r)


def test_existing_authorized_source_adapter_reports_all_eight_without_private_rows(tmp_path):
    from test_species_environment_existing_inputs import files
    from analysis.v4.species_environment_scan import sha256
    p=files(tmp_path);output=tmp_path/'rf_aggregate.json'
    r=rf.run_existing_rf(**p,source_receipt_sha256=sha256(p['source_receipt']),
                         output=output,trees=10,repeats=1,permutation_repeats=1,
                         expected_n=150,geography=False)
    assert output.is_file() and len(r['results'])==8
    assert all(x['status'] in ('insufficient_species_support','source_trait_columns_unavailable') for x in r['results'])
    s=output.read_text()
    for private in ('private_0','taxon_00','"accepted_key"','"obs_id"','"longitude"'):
        assert private not in s


def test_same_observation_id_cannot_be_counted_twice():
    f=fake(n_taxa=60)
    f.loc[len(f)] = f.iloc[0]
    with pytest.raises(ValueError,match='RF_SOURCE_OBSERVATION_IDS_INVALID'):
        rf.aggregate_species(f,'angle')


def test_correlated_bio12_and_bio18_grouped_together_in_training_only():
    f=fake(n_taxa=60)
    a=rf.aggregate_species(f,'angle')
    d=a[['chelsa_bio12','chelsa_bio18','chelsa_bio01']].to_numpy(float)
    cs=rf._corr_groups(d,['chelsa_bio12','chelsa_bio18','chelsa_bio01'])
    assert [0,1] in cs
    assert [2] in cs


def test_constant_trait_reports_not_evaluated_not_missing_geography():
    f=fake(n_taxa=55,geo=True)
    out=rf.evaluate_trait(rf.aggregate_species(f,'hue'),'hue',trees=10,repeats=1,permutation_repeats=1,geography=True)
    assert out['status']=='insufficient_response_variation'
    assert out['geographic_holdout']['status']=='not_evaluated_no_rf_fit'
