"""Direct S2 annual BIO12 / warm BIO18 vs photographed head angle.

This completes the rainfall-angle question omitted from the first combined
Part I+II comparison. Both fits use the exact same old V5 source photos,
same 1-degree cells and 20-degree region blocks; predictors differ only in
BIO01/BIO12 (annual) versus BIO05/BIO18 (warm). Both retain BIO04, BIO15,
and RSDS.

Per species two designs are given: pooled across regions (descriptive
geographic association), and geographically within-block demeaned (the
principal sensitivity). Cluster SEs use the 20-degree region as sampling
unit, not photographs or pairwise pseudo-replicates. A range x rain slope
contrast is tested rather than treating "significant in only one range" as
evidence for a difference.

Observational PHOTO-frame angles; not gravity-calibrated, not causal
rain protection, not population genetics. No climate rarity.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v5.annual_warm_s2 import MODELS, _photo_cells

AXES = {
    "annual_plus_seasonality": (
        "chelsa_bio12", "chelsa_bio01", "chelsa_rsds_mean",
        "chelsa_bio04", "chelsa_bio15"),
    "warm_plus_seasonality": (
        "chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean",
        "chelsa_bio04", "chelsa_bio15"),
}
ANGLE = "orientation_angle_degrees_median"
MIN_CELLS_PER_REGION = 8
MIN_REGIONS_PER_RANGE = 3
CLUSTER_CONDITION_LIMIT = 1.e8


def _supported_cells(cells: pd.DataFrame, species: str) -> pd.DataFrame:
    """Only original one-degree cells in 20-degree regions with >=8 cells."""
    d = cells.loc[cells.taxon_name.eq(species)].copy()
    region_count = d.groupby(["range_state","block"]).block.transform("size")
    d = d.loc[region_count >= MIN_CELLS_PER_REGION].copy()
    d["cluster"] = d.range_state.astype(str)+"_"+d.block.astype(str)
    support=d.groupby("range_state").cluster.nunique().to_dict()
    if (len(d)<120 or any(support.get(state, 0)<MIN_REGIONS_PER_RANGE
                       for state in ("native","introduced"))):
        raise ValueError("S2_RAIN_ANGLE_REGIONAL_SUPPORT_TOO_SMALL_"+species)
    if d[["range_state","block","i","j"]].duplicated().any():
        raise ValueError("S2_REGION_CELL_NOT_UNIQUE")
    return d


def _fit(d: pd.DataFrame, axes, demeaning: bool):
    """A single two-range interaction model, cluster inference by 20deg block."""
    present=d[list(axes)+[ANGLE]].to_numpy(float)
    if not np.isfinite(present).all():
        raise ValueError("S2_RAIN_ANGLE_REQUIRED_VALUES_MISSING")
    x=present[:,:len(axes)]
    y=present[:,-1]
    # Species-wide, all-range 1-degree-cell equally weighted SD.
    mu=x.mean(axis=0);sd=x.std(axis=0)
    yy=y.std()
    if (sd<=1e-10).any() or yy<=1e-10:
        raise ValueError("S2_RAIN_ANGLE_ZERO_PREDICTOR_OR_RESPONSE_SD")
    xx=(x-mu)/sd
    y=(y-y.mean())/yy
    is_intro=d.range_state.eq("introduced").to_numpy(int)
    groups=d.cluster.to_numpy()
    if demeaning:
        # Control all 20degree region-specific intercepts, including range.
        for g in np.unique(groups):
            mask=groups==g
            xx[mask]-=xx[mask].mean(axis=0)
            y[mask]-=y[mask].mean()
        X=np.column_stack((xx,xx*is_intro[:,None]))
    else:
        # Native baseline intercept + five climate partial slopes;
        # introduced intercept and five climate-by-range interactions.
        X=np.column_stack((
            np.ones(len(d)),xx,is_intro,xx*is_intro[:,None]
        ))
    nterms=len(axes)
    offset=0 if demeaning else 1
    interaction=offset+nterms+(0 if demeaning else 1)
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("S2_RAIN_ANGLE_MODEL_NOT_FULL_RANK")
    if np.linalg.cond(X)>=CLUSTER_CONDITION_LIMIT:
        raise ValueError("S2_RAIN_ANGLE_NEAR_COLLINEAR")
    fit=sm.OLS(y,X).fit(
        cov_type="cluster",
        cov_kwds={"groups":groups},
        use_t=True)
    def test(l):
        r=fit.t_test(l)
        lo,hi=np.asarray(r.conf_int(),float).reshape(-1)
        return {
            "beta":float(np.asarray(r.effect).item()),
            "ci95_pointwise":[float(lo),float(hi)],
            "p_pointwise":float(np.asarray(r.pvalue).item()),
        }
    tests={}
    for idx,name in enumerate(("rain","temperature","radiation",
                               "temperature_seasonality",
                               "precipitation_seasonality")):
        contrast=np.zeros(X.shape[1])
        contrast[offset+idx]=1
        native=test(contrast)
        contrast[interaction+idx]=1
        introduced=test(contrast)
        diff=np.zeros(X.shape[1]);diff[interaction+idx]=1
        tests[name]={
            "native":native,
            "introduced":introduced,
            "introduced_minus_native":test(diff),
        }
    return {
        "model":"within_20degree_region" if demeaning else "pooled_spatial_regions",
        "n_geographic_region_clusters":int(d.cluster.nunique()),
        "n_1degree_cells":int(len(d)),
        "n_native_cells":int((is_intro==0).sum()),
        "n_introduced_cells":int((is_intro==1).sum()),
        "n_native_regions":int(d.loc[d.range_state.eq("native"),"cluster"].nunique()),
        "n_introduced_regions":int(d.loc[d.range_state.eq("introduced"),"cluster"].nunique()),
        "predictor_standardisation":"equal_weight_per_eligible_1degree_cell_both_ranges",
        "response_standardisation":"photographic_angle_cell_equal_both_ranges",
        "clustered_se_for_geographic_20degree_region":True,
        "condition_number":float(np.linalg.cond(X)),
        "model_comparison_same_cells":True,
        "rain_partial_slope":tests["rain"],
        "all_axis_slopes":tests,
        "no_temporal_photo_month_control":True,
        "no_causal_rain_exposure_or_heritable_adaptation_inference":True,
    }


def run(photo_wild:pd.DataFrame, measured:pd.DataFrame) -> dict:
    _,cells=_photo_cells(photo_wild,measured)
    output={
        "status":"SOURCE_VERIFIED_DIRECT_S2_RAIN_ANGLE_FULL_YEAR_VS_WARM",
        "description":"original_V5_same_photo_eligible_1degree_cells",
        "annual_vs_warm_models":{k:list(v) for k,v in AXES.items()},
        "photographic_frame_head_angle_not_true_vertical":True,
        "effects_not_mechanistic_or_evolutionary":True,
        "comparisons_posthoc_and_exploratory":True,
        "results":[],
    }
    for species in ("Cirsium arvense","Cirsium vulgare"):
        d=_supported_cells(cells,species)
        for name,axes in AXES.items():
            if set(axes)!=set(MODELS[name]):
                raise ValueError("S2_CLIMATE_AXES_NOT_SAME_AS_FULL_NICHE")
            for within in (False,True):
                fit=_fit(d,axes,within)
                output["results"].append({
                    "species":species,
                    "climate_model":name,
                    "rain_variable":"BIO12" if "annual" in name else "BIO18",
                    **fit,
                })
    # One family per geographic specification: 2 species × 2 climate models,
    # testing the *range difference* in rain-vs-angle partial slopes.
    for which in ("within_20degree_region","pooled_spatial_regions"):
        rows=[v for v in output["results"] if v["model"]==which]
        if len(rows)!=4:raise ValueError("S2_DIRECT_RAIN_PRIMARY_FAMILY_NOT_FOUR")
        adjusted=multipletests([
            v["rain_partial_slope"]["introduced_minus_native"]["p_pointwise"]
            for v in rows],method="holm")[1]
        for v,p in zip(rows,adjusted):
            v["holm_4_species_x_climate_for_range_interaction"]=float(p)
    return output
