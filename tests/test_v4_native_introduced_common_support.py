"""S2 source-verified climate common-support and falsification tests.

Synthetic only. No true introduced-range fitness, photo orientation error,
plasticity, phylogeny or causal evolution results are invented.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_introduced_climate_transfer import (
    dated_and_balanced, standardization
)
from analysis.v4.native_introduced_common_support import (
    REFERENCE_RUN, REFERENCE_RESULTS_SHA256,
    overlap_rectangle, report_support, compare_on_support,
    continent_sensitivity,
)


def toy_climates(n=400, nonoverlap_radiation=False):
    rng=np.random.default_rng(20261021)
    out={}
    for state in ("native","introduced"):
        rows=[]
        for i in range(n):
            cell=i%40
            features={p:float(rng.normal()) for p in VARIABLES}
            features["chelsa_bio12"]=float(800+200*rng.normal())
            features["chelsa_rsds_mean"]=(
                float(rng.normal()+8) if state=="introduced" and nonoverlap_radiation
                else float(rng.normal())
            )
            value=float(80+.015*features["chelsa_bio12"]
                        +.7*features["chelsa_bio01"]+rng.normal(scale=12))
            row={
                "obs_id":state+"-"+str(i),
                "user_id":"human-"+state+"-"+str(i%167),
                "taxon_name":"Cirsium arvense",
                "range_state":state,
                "tdwg1":1 if state=="native" else 7,
                "latitude":-55.+(cell%20)*5.5,
                "longitude":-160.+(cell//20)*90.,
                "observed_on":f"202{1+i%4}-{(i%12+1):02d}-01",
                "orientation_angle_degrees_median":value,
            }
            row.update(features)
            rows.append(row)
        out[state]=dated_and_balanced(pd.DataFrame(rows))
    return out["native"],out["introduced"]


def test_BIO12_overlap_cannot_excuse_four_unmatched_covariates():
    n,i=toy_climates(nonoverlap_radiation=True)
    rain=overlap_rectangle(n,i,("chelsa_bio12",))
    five=overlap_rectangle(n,i,tuple(VARIABLES))
    assert rain["climate_intersections_nonempty"]
    assert rain["native_rows_in_overlap"]>50
    assert rain["introduced_rows_in_overlap"]>50
    assert not five["climate_intersections_nonempty"]
    assert five["native_rows_in_overlap"]==five["introduced_rows_in_overlap"]==0
    assert five["limits"]["chelsa_rsds_mean"]["intersection"] is None
    assert not any("mask" in x for x in report_support(five))


def test_full_five_overlap_fit_is_conditional_and_not_evolutionary():
    n,i=toy_climates()
    rectangle=overlap_rectangle(n,i,tuple(VARIABLES))
    assert rectangle["climate_intersections_nonempty"]
    f=pd.concat([n,i],ignore_index=True)
    scales=standardization(f)
    report=compare_on_support(n,i,scales,rectangle)
    assert report["status"]=="ESTIMATED_CONDITIONAL_MARGINAL_COMMON_SUPPORT_SENSITIVITY"
    for state in ("native","introduced"):
        assert report[state+"_support"]["unique_original_5degree_cells"]>=15
        assert report["refitted_models"][state]["status"]=="estimated_five_exposure_cell_cluster_WLS"
    assert report["slope_difference"]["not_equivalence_testing"]
    assert report["native_trained_transfer"]["climatic_overlap_gated"] is False
    assert report["support_does_not_remove_camera_technical_error"] is True


def test_no_result_is_misread_as_true_zero_after_overlap_screening():
    n,i=toy_climates(nonoverlap_radiation=True)
    scales=standardization(pd.concat([n,i],ignore_index=True))
    rect=overlap_rectangle(n,i,tuple(VARIABLES))
    result=compare_on_support(n,i,scales,rect)
    assert result["status"]=="STOP_INSUFFICIENT_CLIMATE_COMMON_SUPPORT"
    assert result["slope_difference"] is None
    assert result["native_trained_transfer"] is None
    assert result["introduced_support"]["n_original_observer_cell_balanced"]==0


def test_region_sensitivity_reports_failure_instead_of_changing_WCVP_native_status():
    n,i=toy_climates(nonoverlap_radiation=True)
    scales=standardization(pd.concat([n,i],ignore_index=True))
    out=continent_sensitivity(n,i,scales)
    assert set(out)=={"North America"}
    assert out["North America"]["status"]=="STOP_INSUFFICIENT_CLIMATE_COMMON_SUPPORT"
    assert "introduced_region_observer_cell_balanced_rows" in out["North America"]


def test_reject_malformed_climate_columns_ambiguity_and_missingness():
    n,i=toy_climates(n=80)
    with pytest.raises(ValueError,match="unique non-empty"):
        overlap_rectangle(n,i,("chelsa_bio12","chelsa_bio12"))
    with pytest.raises(ValueError,match="Frozen climate predictors"):
        overlap_rectangle(n,i,("rain_unknown",))
    faulty=n.copy(); faulty.loc[faulty.index[0],"chelsa_bio12"]=np.nan
    with pytest.raises(ValueError,match="missing climate"):
        overlap_rectangle(faulty,i,("chelsa_bio12",))


def test_frozen_prior_executed_numerical_source_is_not_same_as_new_common_support():
    assert REFERENCE_RUN==39539209624
    assert REFERENCE_RESULTS_SHA256==(
        "a95e0eaa085cc4cbd55a8302f5d527e26e3800761f7bb043ae9e65bfecdbdb26"
    )
