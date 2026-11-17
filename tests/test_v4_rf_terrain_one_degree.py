"""Synthetic validation only; never use original S1 identity or image data."""
import numpy as np
import pytest
from analysis.v4.rf_extended_environment import BIO9,TERRAIN5
from analysis.v4.rf_archived_atmosphere import ATMOSPHERE23
from analysis.v4.rf_terrain_one_degree import score_new_geographic_grain
from test_v4_rf_multitrait_predictability import tiny_synthetic


def synthetic_wider_predictors():
    d=tiny_synthetic(n_taxa=30,n_cells=8,reps=6)
    rng=np.random.default_rng(14)
    for name in BIO9:
        d[name]=rng.normal(size=len(d))
    for name in ("chelsa_vpd_mean","chelsa_sfcwind_mean"):
        d[name]=rng.normal(size=len(d))
    for name in TERRAIN5:
        d[name]=rng.normal(size=len(d))
    assert all(n in d for n in ATMOSPHERE23+TERRAIN5)
    return d


def test_1degree_models_holdout_taxa_and_sector_without_predictor_name_leak():
    d=synthetic_wider_predictors()
    report=score_new_geographic_grain(d,constructs=("angle",),trees=8,repeats=1)
    assert len(report)==2
    assert {r['generalisation'] for r in report}=={'unseen_taxa_balanced','unseen_longitude_sectors'}
    for r in report:
        assert r['n_taxa']==30
        assert r['n_cells']==240
        assert list(r['models'])==[
            'chelsa_vpd_wind23','chelsa_vpd_wind23_plus_worldclim_terrain28']
        assert np.isfinite(r['paired_skill_gain_upper_minus_lower'])
    assert all(5 <= x['test_taxa'] <= 6 for x in
         report[0]['models']['chelsa_vpd_wind23']['fold_diagnostics_aggregate_only'])


def test_only_predeclared_1degree_min3_protocol_accepted():
    d=synthetic_wider_predictors()
    with pytest.raises(ValueError,match='UNDECLARED_1_DEGREE'):
        score_new_geographic_grain(d,degrees=2,constructs=('angle',),trees=8,repeats=1)
    with pytest.raises(ValueError,match='UNDECLARED_1_DEGREE'):
        score_new_geographic_grain(d,min_photos_per_cell=1,constructs=('angle',),trees=8,repeats=1)