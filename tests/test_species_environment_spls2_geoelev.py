"""Synthetic, non-protected test for original V4 S1 climate-axis/geo/elevation model."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4.species_environment_spls2 import CLIMATE, YCOLS, SparsePLS2
from analysis.v4.species_environment_spls2_geoelev import (
    aggregate_joint_elevation, project_scores, compare_cv, run, _make_design,
    _standardized_y, FORBIDDEN,
)


def make_source(n_taxa=48, n_photos=7, seed=111):
    rng=np.random.default_rng(seed)
    rows=[]
    for s in range(n_taxa):
        sector=s%3
        sectorlon=[-90.,20.,135.][sector]
        z=rng.normal(size=22)
        e=700+300*z[0]+rng.normal(0,90)
        for j in range(n_photos):
            v=z+rng.normal(0,.10,size=22)
            yi=np.array([75+4*z[11]+rng.normal(),65-3*z[7]+rng.normal(),
                         30-2*z[15]+rng.normal(),np.sin(.4*z[7]),np.cos(.4*z[7]),
                         1.2+.15*z[0]+rng.normal(0,.05),.7+.05*z[10]+rng.normal(0,.03),
                         .8+.03*z[11]+rng.normal(0,.02),.15+.01*z[0]+rng.normal(0,.007)])
            row={'obs_id':f'obs-{s:03d}-{j:02d}','accepted_key':f'taxon_{s:03d}',
                 'latitude':-35+2*sector+rng.normal(0,1),
                 'longitude':sectorlon+rng.normal(0,3),
                 'topo_elevation':e+rng.normal(0,20)}
            row.update(dict(zip(CLIMATE,v)))
            row.update(dict(zip(YCOLS,yi)))
            row.update({q+'__eligible':True for q in YCOLS})
            rows.append(row)
    return pd.DataFrame(rows)


def test_joint_source_and_elevation_cohort():
    source=make_source()
    taxa,receipt=aggregate_joint_elevation(source)
    assert len(taxa)==48
    assert receipt['n_source_species_before_elevation_gate']==48
    assert receipt['n_source_species_matching_elevation_gate']==48
    assert receipt['matched_joint_source_photo_contributions']==len(source)
    assert len(taxa.columns)>22
    assert taxa.groupby('sector').size().to_dict()=={0:16,1:16,2:16}


def test_source_missing_elevation_gate_uses_exact_eligible_cohort():
    source=make_source()
    source.loc[source.accepted_key=='taxon_005','topo_elevation']=np.nan
    taxa,meta=aggregate_joint_elevation(source)
    assert len(taxa)==47
    assert meta['n_source_species_before_elevation_gate']==48
    assert meta['n_source_species_matching_elevation_gate']==47


def test_source_excludes_ineligible_photo_elevation():
    source=make_source()
    mask=(source.accepted_key=='taxon_001')&(source.obs_id.str.endswith('00'))
    source.loc[mask,'corolla_lab_chroma__eligible']=False
    source.loc[mask,'topo_elevation']=100000.
    taxa,meta=aggregate_joint_elevation(source)
    assert meta['matched_joint_source_photo_contributions']==len(source)-1
    got=taxa.loc[taxa._taxon_private=='taxon_001','topo_elevation'].iloc[0]
    expected=source.loc[(source.accepted_key=='taxon_001')&(~mask),'topo_elevation'].mean()
    assert abs(got-expected)<1e-10


def test_same_taxon_support_for_all_models():
    source=make_source()
    taxa,_=aggregate_joint_elevation(source)
    x=compare_cv(taxa,repeats=1)
    assert x['n_partitions']==5
    assert x['geographic_holdout']['status']=='estimated'
    assert x['geographic_holdout']['n_partitions']==3
    assert len(x['model_skill_against_training_mean'])==5
    assert len(x['model_skill_by_construct'])==8
    assert all(len(v)==5 for v in x['model_skill_by_construct'].values())
    assert sum(x['selected_sparsity_frequency'].values())==5
    assert all(k in x['incremental_skill_relative_to_train_mean']
               for k in ('climate_after_geography','elevation_after_climate_geography'))


def test_projection_honours_pls2_deflation():
    rng=np.random.default_rng(10);X=rng.normal(size=(80,22));Y=rng.normal(size=(80,9))
    Y[:,0]=X[:,0]*2+Y[:,0]
    m=SparsePLS2(16,2).fit(X,Y)
    projection=project_scores(m,X)
    standardized=(X-m.xmean)/m.xstd
    expected=standardized@m.weights@np.linalg.pinv(m.xloadings.T@m.weights)
    assert projection.shape==(80,2)
    assert np.allclose(projection,expected)
    assert np.linalg.norm(projection-standardized@m.weights)>1e-6


def test_training_only_projection_does_not_depend_on_test_response():
    rng=np.random.default_rng(12)
    X=rng.normal(size=(80,22));Y=rng.normal(size=(80,9))
    model=SparsePLS2(11,2).fit(X[:60],Y[:60])
    score=project_scores(model,X[60:])
    Y[60:]=rng.normal(0,1000,size=(20,9))
    assert np.array_equal(score,project_scores(model,X[60:]))


def test_geo_and_elevation_same_train_dimensions():
    rng=np.random.default_rng(14)
    tr={'climate':rng.normal(size=(30,2)),'geography':rng.normal(size=(30,8)),'elevation':rng.normal(size=(30,1))}
    te={'climate':rng.normal(size=(20,2)),'geography':rng.normal(size=(20,8)),'elevation':rng.normal(size=(20,1))}
    expected={'climate_axes':2,'geography_only':8,'climate_axes_plus_geography':10,
              'geography_plus_elevation':9,'climate_axes_plus_geography_plus_elevation':11}
    from analysis.v4.species_environment_spls2_geoelev import MODELS
    for k,v in MODELS.items():
        x,y=_make_design(tr,te,v)
        assert x.shape==(30,expected[k])
        assert y.shape==(20,expected[k])


def test_hue_joint_variance_counted_one():
    rng=np.random.default_rng(15)
    y=rng.normal(size=(80,9))
    x,z=_standardized_y(y[:60],y[60:])
    assert np.isclose(x[:,3:5].var(axis=0,ddof=1).sum(),1,rtol=1e-12)
    assert z.shape==(20,9)


def test_run_emits_anonymous_aggregates():
    d=run(make_source(),repeats=1)
    assert d['original_input_observations']==48*7
    assert d['species_support']['n_source_species_matching_elevation_gate']==48
    assert d['heldout']['geographic_holdout']['status']=='estimated'
    assert d['no_confirmatory_p_values'] is True
    assert d['original_frozen_S1_S2_H3c_unchanged'] is True
    keys=set()
    def walk(o):
        if isinstance(o,dict):
            for k,v in o.items():keys.add(k);walk(v)
        if isinstance(o,list):
            for v in o:walk(v)
    walk(d)
    assert not keys.intersection(FORBIDDEN)
    json.dumps(d,allow_nan=False)


def test_no_leakage_from_y_postphotos():
    d=make_source()
    d.loc[d.accepted_key=='taxon_010','corolla_hue_sin__eligible']=False
    t,a=aggregate_joint_elevation(d)
    assert len(t)==47
    assert a['n_source_species_matching_elevation_gate']==47


def test_strict_geo_source_validation():
    d=make_source()
    d.loc[0,'latitude']=200.
    with pytest.raises(ValueError,match='COORDINATES'):
        aggregate_joint_elevation(d)


def test_repeats_validation():
    d,_=aggregate_joint_elevation(make_source())
    with pytest.raises(ValueError,match='REPEATS'):
        compare_cv(d,repeats=0)
