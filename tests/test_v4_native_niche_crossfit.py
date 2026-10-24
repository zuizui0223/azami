"""Synthetic-only tests of spatially cross-fitted native climate envelopes.

True results are calculated on authenticated original 665115-source Actions.
These tests never invent a real photometric/evolutionary association.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_niche_crossfit import (
    FOLDS, SCHEME_VERSION, cross_fitted_geographic_novelty,
    geographic_cell_profile, _weighted_limits
)


def native_introduced_toy(*, native_cells=45, introduced_cells=26,
                          per_cell=5, rain_shift=0.0, rsds_shift=0.0):
    def make(ncells, offset, rain, sun):
        data=[]
        rng=np.random.default_rng(20261021+offset)
        for cell in range(ncells):
            # Each region samples somewhat varied native climate cells.
            season=np.sin(cell*.34)
            for j in range(per_cell):
                row={
                    "_cell":cell+offset,
                    "obs_id":f"synthetic-{offset}-{cell}-{j}",
                    "user_id":f"synthetic-human-{offset}-{cell}-{j%3}",
                    "orientation_angle_degrees_median":float(45+30*rng.random()),
                    "chelsa_bio01":float(10 + 5*season+rng.normal(scale=.4)),
                    "chelsa_bio04":float(500+170*np.cos(cell*.45)+rng.normal(scale=10)),
                    "chelsa_bio12":float(600+95*season+rain+rng.normal(scale=11)),
                    "chelsa_bio15":float(33+10*np.sin(cell*.19)+rng.normal(scale=1)),
                    "chelsa_rsds_mean":float(12+1.5*np.cos(cell*.33)+sun+rng.normal(scale=.3))
                }
                data.append(row)
        return pd.DataFrame(data)
    return (
        make(native_cells,0,0.,0.),
        make(introduced_cells,1000,rain_shift,rsds_shift)
    )


def test_crossfit_detects_strong_climate_shift_without_looking_at_angle_values():
    n,i=native_introduced_toy(rain_shift=350.,rsds_shift=4.)
    a=cross_fitted_geographic_novelty(
        n,i,source_name="Cirsium synthetic",repeats=8,folds=5
    )
    assert a["status"].startswith("EXECUTED_SOURCE_CLIMATE")
    assert a["schema_version"]==SCHEME_VERSION
    assert (a["original_native_equal_sin_5degree_cells"],
            a["original_introduced_equal_sin_5degree_cells"])==(45,26)
    assert a["native_reference_trained_on_nonheldout_cells"]
    assert a["introduced_and_native_evaluated_without_angle_values"]
    assert a["equal_total_weight_per_original_cell"]
    assert a["folds"]==FOLDS
    expected=set(("annual_precipitation_only","four_without_rsds",
                  "five_original_climate_predictors",
                  "rsds_only_technical_extrapolation"))
    assert set(a["results"])==expected
    for key,item in a["results"].items():
        introduced=item["introduced_novel_fraction_partition_q025_q50_q975"]
        native=item["genuinely_heldout_native_novel_fraction_partition_q025_q50_q975"]
        diffs=item["introduced_minus_heldout_native_fraction_partition_q025_q50_q975"]
        assert all(np.isfinite(x) for x in introduced+native+diffs)
        assert 0<=introduced[0]<=introduced[1]<=introduced[2]<=1
        assert 0<=native[0]<=native[1]<=native[2]<=1
        assert diffs[0]<=diffs[1]<=diffs[2]
        assert item["not_a_sampling_confidence_interval"]
        assert item["head_angle_never_entered_climate_reference_or_fold_selection"]
    assert a["results"]["annual_precipitation_only"][
        "introduced_minus_heldout_native_fraction_partition_q025_q50_q975"][1] > .2
    assert a["results"]["five_original_climate_predictors"][
        "fraction_of_random_partitions_with_positive_excess"] == 1.

    # Changing every orientation without changing any climate MUST have no effect
    # on climate novelty or native reference partition.
    m=n.copy();t=i.copy()
    m["orientation_angle_degrees_median"]=1.
    t["orientation_angle_degrees_median"]=179.
    b=cross_fitted_geographic_novelty(
        m,t,source_name="Cirsium synthetic",repeats=8,folds=5
    )
    assert a==b


def test_geo_holdout_is_deterministic_under_row_order_shuffle():
    n,i=native_introduced_toy()
    a=cross_fitted_geographic_novelty(
        n,i,source_name="Cirsium synthetic",repeats=6,folds=5
    )
    b=cross_fitted_geographic_novelty(
        n.sample(frac=1,random_state=21),i.sample(frac=1,random_state=31),
        source_name="Cirsium synthetic",repeats=6,folds=5
    )
    assert a["results"].keys()==b["results"].keys()
    for key in a["results"]:
        for field in (
            "introduced_novel_fraction_partition_q025_q50_q975",
            "genuinely_heldout_native_novel_fraction_partition_q025_q50_q975",
            "introduced_minus_heldout_native_fraction_partition_q025_q50_q975"
        ):
            assert a["results"][key][field]==pytest.approx(b["results"][key][field],abs=1e-12)


def test_not_many_regions_does_not_turn_gap_into_zero_evolution():
    n,i=native_introduced_toy(native_cells=12,introduced_cells=10)
    a=cross_fitted_geographic_novelty(n,i,source_name="test",repeats=6)
    assert a["status"]=="STOP_INSUFFICIENT_SOURCE_GEOGRAPHIC_CELLS"
    assert a["not_equivalent_to_no_niche_shift"] is True
    assert "results" not in a


def test_private_original_identifiers_are_not_emitted():
    n,i=native_introduced_toy()
    a=cross_fitted_geographic_novelty(n,i,source_name="Cirsium synthetic",
                                      repeats=3,folds=3)
    payload=json.dumps(a)
    assert "synthetic-human" not in payload
    assert "synthetic-0-0" not in payload
    assert '"native_mask"' not in payload
    assert "orientation_angle_degrees_median" not in payload
    assert "not a 95%" in payload or "not a 95%" in str(a).lower() or any(
        "not a 95%" in x.lower() for x in a["interpretation_limitations"]
    )


def test_source_validation_and_cell_equal_quantiles():
    n,i=native_introduced_toy()
    unique,inverse,weights=geographic_cell_profile(n)
    assert len(unique)==45
    assert np.allclose(np.bincount(inverse,weights=weights),np.ones(45))
    assert len(_weighted_limits(n["chelsa_bio12"].to_numpy(float),weights))==2
    g=n.copy()
    g.loc[g.index[0],"chelsa_bio12"]=np.nan
    with pytest.raises(ValueError,match="Non-finite|Missing climate"):
        cross_fitted_geographic_novelty(g,i,source_name="test",repeats=3)
    x=n.copy()
    x.loc[x.index[0],"obs_id"]=x.iloc[1]["obs_id"]
    with pytest.raises(ValueError,match="not one-to-one"):
        cross_fitted_geographic_novelty(x,i,source_name="test",repeats=3)
    with pytest.raises(ValueError,match="Exact preselected source taxon"):
        cross_fitted_geographic_novelty(n,i,source_name="",repeats=3)
