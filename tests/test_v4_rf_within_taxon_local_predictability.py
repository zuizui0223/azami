"""Synthetic tests: no original S1 records or remote files required."""
import numpy as np
import pandas as pd
import pytest
from analysis.v4.rf_worldclim_terrain import MODEL_SET
from analysis.v4.rf_within_taxon_local_predictability import heldout_spatial_folds, within_taxon_heldout


def make_cells(n_taxa=12, blocks=10, n_per_block=3, noise=0.25):
    rng=np.random.default_rng(419)
    predictors=MODEL_SET['chelsa_vpd_wind23_plus_worldclim_terrain28']
    out=[]
    for s in range(n_taxa):
        for b in range(blocks):
            for c in range(n_per_block):
                env=rng.normal(size=len(predictors))
                # Correlated synthetic environmental effect with taxon intercept
                trait=10*s + 3*env[0] - 2*env[3] + rng.normal(scale=noise)
                out.append({'taxon':str(s),'lat_cell':20+b*5+c,
                            'lon_cell':50+b*5+c,
                            'orientation_image_vertical_angle':trait,
                            **dict(zip(predictors,env))})
    return pd.DataFrame(out)


def test_no_shared_5_degree_block_in_spatial_holdout():
    cells=make_cells()
    for a,b in heldout_spatial_folds(cells):
        blocks=(cells.lat_cell//5)*72+cells.lon_cell//5
        assert set(blocks.iloc[a]).isdisjoint(blocks.iloc[b])


def test_models_predict_same_taxon_local_deviations_with_train_only_baselines():
    cells=make_cells()
    result=within_taxon_heldout(cells,'angle',trees=20)
    assert result['n_taxa']==12
    assert result['train_species_means_never_use_heldout_trait_values']
    assert len(result['folds_anonymous_counts'])==5
    assert all(v['test_taxa_with_training_baseline']>=5 for v in result['folds_anonymous_counts'])
    assert set(result['models'])==set(MODEL_SET)
    assert all(np.isfinite(v['skill_over_species_training_means']) for v in result['models'].values())
    assert result['models']['chelsa_vpd_wind23']['skill_over_species_training_means']>0
    assert 'taxon' not in result['models']


def test_no_species_lookup_uses_heldout_outcome():
    data=make_cells(n_taxa=12,blocks=10,n_per_block=2)
    train,test=heldout_spatial_folds(data)[0]
    result=within_taxon_heldout(data,'angle',trees=8)
    modified=data.copy()
    # Changing only test y can change scoring, not training fit; don't assert
    # same skill but source training taxon mean must not be obtained from test y.
    modified.loc[test,'orientation_image_vertical_angle'] += 100
    assert within_taxon_heldout(modified,'angle',trees=8)['models']['chelsa_vpd_wind23']['skill_over_species_training_means'] != result['models']['chelsa_vpd_wind23']['skill_over_species_training_means']


def test_invalid_scope_fails_closed():
    frame=make_cells()
    with pytest.raises(ValueError,match='UNREGISTERED'):
        within_taxon_heldout(frame,'hue',trees=10)
    with pytest.raises(ValueError,match='INSUFFICIENT'):
        heldout_spatial_folds(frame.iloc[:4],n_splits=5)