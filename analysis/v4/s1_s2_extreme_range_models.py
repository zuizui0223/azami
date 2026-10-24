"""V4 S2 REAL native-vs-introduced climate gradients on original measured traits.

Within the two original introduced-range species, compare annual BIO12 and
warmest-quarter BIO18 on EXACTLY matched WCVP wild, assessor-eligible
photographs. Baseline controls include temperature/radiation; elevation is a
separately declared sensitivity, not silently treated as independent climate.
Each original observer x 5-degree cell contributes at most one photo record.
Cell-cluster errors are observational sensitivity, not independent evolution.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.api as sm

from analysis.v4.native_spatial import cells
from analysis.v4.s1_s2_extreme_lego import (
    ALL_ENV, FEATURES, diagnostic_matrix, require)

S2_MODELS = {
    "annual_BIO12": ["chelsa_bio12","chelsa_bio01","chelsa_rsds_mean"],
    "warmest_quarter_BIO18": ["chelsa_bio18","chelsa_bio01","chelsa_rsds_mean"],
    "annual_BIO12_plus_elevation": ["chelsa_bio12","chelsa_bio01","chelsa_rsds_mean","topo_elevation"],
    "warmest_quarter_BIO18_plus_elevation": ["chelsa_bio18","chelsa_bio01","chelsa_rsds_mean","topo_elevation"],
    "dry_quarter_BIO17_plus_elevation": ["chelsa_bio17","chelsa_bio01","chelsa_rsds_mean","topo_elevation"],
    "warmest_month_BIO05": ["chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","topo_elevation"],
    "coldest_month_BIO06": ["chelsa_bio18","chelsa_bio06","chelsa_rsds_mean","topo_elevation"],
}
MIN_INDEPENDENT_CELLS=15
MIN_PHOTOGRAPHERS=10
MIN_PHOTOS=100


def _get_dated_balanced(frame):
    if frame.obs_id.duplicated().any() or frame.user_id.isna().any():
        raise ValueError("S2_ORIGINAL_PHOTO_OR_OBSERVER_ID_INVALID")
    d=frame.copy()
    d["_cell"]=cells(d.latitude.to_numpy(float),d.longitude.to_numpy(float))
    dt=pd.to_datetime(d.observed_on,errors="coerce",format="mixed")
    d=d.loc[dt.notna()].copy()
    if len(d)==0:
        return d
    dt=dt.loc[d.index]
    d["_month"]=dt.dt.month.to_numpy()
    d["_year"]=dt.dt.year.to_numpy()
    d=d.sort_values("obs_id",kind="stable").drop_duplicates(
        ["user_id","_cell"],keep="first").copy()
    return d


def _formula(frame,columns):
    x=frame[columns].to_numpy(float)
    if not np.isfinite(x).all():
        raise ValueError("S2_MISSING_CLIMATE_AFTER_MATCHED_SAMPLE")
    means=x.mean(axis=0)
    std=x.std(axis=0)
    if np.any(std<1e-12):
        raise ValueError("S2_ENVIRONMENT_CONSTANT")
    z=(x-means)/std
    month=frame["_month"].to_numpy(float)
    year=frame["_year"].to_numpy(float)
    ystd=year.std()
    nuisance=np.column_stack([
        np.sin(2*np.pi*month/12),
        np.cos(2*np.pi*month/12),
        (year-year.mean())/ystd if ystd>0 else np.zeros(len(year))])
    design=np.column_stack([np.ones(len(x)),z,nuisance])
    if np.linalg.matrix_rank(design)<design.shape[1]:
        raise ValueError("S2_DESIGN_RANK_DEFICIENT")
    return design,means,std


def model_one(frame,outcome,variables):
    if len(frame)<MIN_PHOTOS or frame._cell.nunique()<MIN_INDEPENDENT_CELLS or \
       frame.user_id.nunique()<MIN_PHOTOGRAPHERS:
        return {"status":"INSUFFICIENT_RANGE_GEOGRAPHIC_SUPPORT",
                "n":len(frame),"cells":int(frame._cell.nunique()),
                "observers":int(frame.user_id.nunique())}
    coldiagnostic=diagnostic_matrix(frame[variables],variables)
    if coldiagnostic["status"]!="OK" or \
       coldiagnostic["max_abs_pair_r"]>.8 or \
       max(coldiagnostic["vif"].values())>5:
        return {"status":"UNRESOLVED_S2_CLIMATE_COLLINEARITY",
                "n":len(frame),"diagnostic":coldiagnostic}
    try:
        design,means,std=_formula(frame,variables)
    except (ValueError,np.linalg.LinAlgError) as exc:
        return {"status":"S2_EXPOSURE_OR_DATE_DESIGN_NOT_IDENTIFIABLE",
                "reason":str(exc)[:160],"n":len(frame),
                "cells":int(frame._cell.nunique())}
    y=frame[outcome].to_numpy(float)
    if not np.isfinite(y).all() or y.std()<1e-12:
        return {"status":"INVALID_OR_CONSTANT_S2_TRAIT_VALUE","n":len(frame)}
    cell=frame._cell.to_numpy()
    ncell=frame.groupby("_cell").size()
    weights=frame._cell.map(1/ncell).to_numpy(float)
    fit=sm.WLS(y,design,weights=weights).fit(
        cov_type="cluster",cov_kwds={"groups":cell},use_t=False)
    return {"status":"EXPLORATORY_OBSERVED_S2_CLIMATE_GRADIENT",
            "n":len(frame),"cells":int(len(ncell)),
            "observers":int(frame.user_id.nunique()),
            "same_source_cohort_before_models":True,
            "trait":outcome,"predictors":variables,
            "hydric_predictor":variables[0],
            "hydric_coefficient_per_sample_SD":float(fit.params[1]),
            "hydric_coefficient_per_unit":float(fit.params[1]/std[0]),
            "hydric_source_SD":float(std[0]),
            "hydric_cluster_se_per_SD":float(fit.bse[1]),
            "hydric_cluster_p_uncorrected":float(fit.pvalues[1]),
            "AIC_same_matched_range_rows":float(fit.aic),
            "seasonal_month_sin_cos_and_year_adjusted":True,
            "5degree_cell_cluster_robust_not_spatial_causal":True}


def analyse_s2_range_models(focal):
    required={"obs_id","user_id","observed_on","latitude","longitude","taxon_name",
              "range_state",*ALL_ENV}
    require(required<=set(focal),"S2_ORIGINAL_SOURCE_OR_CLIMATE_LAYERS_MISSING")
    # Exactly seven original climate layers + five baseline ones are required,
    # not models on separately missing or selectively imputed records.
    table={}
    for (species,region),g in focal.groupby(["taxon_name","range_state"]):
        full=np.isfinite(g[ALL_ENV].to_numpy(float)).all(axis=1)
        basis=g.loc[full].copy()
        selection={}
        for key,names in S2_MODELS.items():
            selection[key]={
               "all_source_predictor_diagnostic":diagnostic_matrix(basis[names],names),
               "predictors":names}
        outcomes={}
        for block,features in FEATURES.items():
            if not all(v in basis for v in features):
                outcomes[block]={"status":"ORIGINAL_S2_MEASUREMENT_NOT_AVAILABLE",
                                 "source_columns_required":features}
                continue
            trait_relevant=basis.dropna(subset=features).copy()
            if block=="angle":
                trait_relevant=trait_relevant.loc[
                    trait_relevant[features[0]].between(0,180)].copy()
            # Hue is a 2D circular phenotype, so don't interpret sine and
            # cosine as two independent biological trait discoveries.
            cohort=_get_dated_balanced(trait_relevant)
            scenario_results={}
            for key,names in S2_MODELS.items():
                components={f:model_one(cohort,f,names) for f in features}
                scenario_results[key]={
                    "source_phenotype_unit":block,
                    "joint_circular_hue_needs_multivariate_joint_test":block=="hue",
                    "component_diagnostics_not_independent_endpoints":components}
            outcomes[block]={
                "source_measured_n":len(trait_relevant),
                "date_photographer_5deg_balanced_n":len(cohort),
                "same_cohort_used_for_all_7_s2_climate_scenarios":True,
                "range_models":scenario_results}
        table[species+"/"+region]={
            "original_WCVP_range_source_observations":len(g),
            "complete_all12_source_climate":len(basis),
            "source_only_collinearity":selection,
            "original_visible_trait_environment_models":outcomes}
    return {"status":"S2_NEW_EXPLORATORY_CLIMATE_MODELS_SOURCE_RANGE_STRATIFIED",
            "strata":table,
            "source_strata_are_observation_rows_not_independent_colonizations":True,
            "alternative_annual_vs_warm_precip_no_multiple_test_claim":True,
            "no_physical_plant_genetic_LEGO_or_adaptation_inference":True}
