"""V4 S1 annual rain vs image angle: original-source regional/spatial confounding.

One audited S1 source (170,066 native orientation images / 244 accepted keys)
and the SAME five preselected CHELSA variables BIO1/BIO4/BIO12/BIO15/RSDS.
Test strict 237 noncollided taxon keys and 191 tree-covered keys. Use real
native observation latitude/longitude; no inferred biogeographic-region tags.
Geographic controls are response-blind fixed longitude sectors, or spherical
degree-2 terms. Compare original photo-weighted and geographically equal 5°
taxon-cell source measurements; both analysed with equal taxon weights.

This is spatial *confounding sensitivity* and not spatially blocked causal
inference, ecological fitness, or gravity-calibrated nodding of capitula.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

from analysis.v4.joint_s1_t2_phylogeny import (
    ENDPOINT, EXPOSURES, MODEL_PATH, MAPPING_PATH,
)
from analysis.v4.native_integration_audit import eligible

ORIGINAL_SHA="342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
SCALES=np.array([1.,100.,100.,10.,1.],dtype=float)
SECTOR_NAMES=("W180_W030","W030_E060","E060_E180")
GEO_COLS=("g_x","g_y","g_z","g_xy","g_xz","g_yz",
          "g_x2_minus_y2","g_3z2_minus1")
PREDICTOR_NAMES=("chelsa_bio12","chelsa_bio15")
GEO_ADJUSTMENTS=("none","three_longitude_sector_proportions",
                 "spherical_second_order")
CROP_COHORTS=("clean_native","clean_tree","clean_tree_min10")


def geo_features(lat_deg,lon_deg):
    """Fixed a priori 8-term spherical location basis, no response fitting."""
    lat=np.asarray(lat_deg,dtype=float)
    lon=np.asarray(lon_deg,dtype=float)
    if lat.shape!=lon.shape or not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("Original source positions missing/nonfinite")
    if np.any(abs(lat)>90) or np.any(abs(lon)>180):
        raise ValueError("Original source coordinates outside valid degrees")
    la=np.deg2rad(lat)
    lo=np.deg2rad(lon)
    x=np.cos(la)*np.cos(lo)
    y=np.cos(la)*np.sin(lo)
    z=np.sin(la)
    return np.column_stack((x,y,z,x*y,x*z,y*z,x*x-y*y,3*z*z-1))


def longitude_sector(lon_deg):
    lon=np.asarray(lon_deg,dtype=float)
    if np.any(~np.isfinite(lon)) or np.any(np.abs(lon)>180):
        raise ValueError("Longitude invalid for fixed geographic strata")
    # Three a priori longitudinal sectors, not political/biogeographic realms.
    return np.select((lon< -30,lon<60), (0,1),default=2).astype(np.int8)


def authenticated_source(table):
    original=Path(table)
    info=json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    if (hashlib.sha256(original.read_bytes()).hexdigest()!=ORIGINAL_SHA
        or info["input_sha256"]!=ORIGINAL_SHA
        or info["selection"]["selected"]!=list(EXPOSURES)):
        raise ValueError("Cannot change source S1 cohort or reselect five CHELSA predictors")
    df=pd.read_csv(original,low_memory=False)
    required={"obs_id","accepted_key","taxon_name","latitude","longitude",
             ENDPOINT,ENDPOINT+"__eligible",*EXPOSURES}
    if not required<=set(df) or not df.obs_id.is_unique:
        raise ValueError("Missing authentic native S1 coordinates, phenotype or identity")
    mask=eligible(df,[ENDPOINT]) & np.isfinite(df[list(EXPOSURES)]).all(axis=1)
    if int(mask.sum())!=170066 or int(df.loc[mask,"accepted_key"].nunique())!=244:
        raise ValueError("Frozen S1 eligible population 170066/244 drift")
    work=df.loc[mask,["accepted_key","taxon_name","latitude","longitude",
                       ENDPOINT,*EXPOSURES]].copy()
    work["accepted_key"]=work["accepted_key"].astype(str)
    work["latitude"]=pd.to_numeric(work.latitude,errors="coerce")
    work["longitude"]=pd.to_numeric(work.longitude,errors="coerce")
    valid=(np.isfinite(work.latitude)&np.isfinite(work.longitude)
          &work.latitude.between(-90,90)&work.longitude.between(-180,180))
    n_invalid=int((~valid).sum())
    if n_invalid>len(work)*.1:
        raise ValueError("Too many source S1 image locations unavailable for spatial audit")
    work=work.loc[valid].copy()
    values=geo_features(work.latitude.to_numpy(),work.longitude.to_numpy())
    for j,c in enumerate(GEO_COLS):
        work[c]=values[:,j]
    sectors=longitude_sector(work.longitude.to_numpy())
    for s,n in enumerate(SECTOR_NAMES):
        work[f"sector_{n}"]=(sectors==s).astype(float)
    work["_cell_lat"]=(np.floor((work.latitude+90)/5)).astype(int)
    work["_cell_lon"]=(np.floor((work.longitude+180)/5)).astype(int)
    initial={"n_S1_five_climate_eligible":170066,"n_S1_accepted_keys":244,
             "n_with_valid_source_coordinates":len(work),
             "n_source_invalid_coordinates_excluded":n_invalid,
             "n_taxa_remaining_with_coordinates":int(work.accepted_key.nunique()),
             "geography":"observed native S1 photo locations; no inferred national boundaries",
             "five_selected_exposures":list(EXPOSURES)}
    return work,initial


def taxon_rosters(work,mode):
    if mode not in ("photo_weighted","equal_5deg_spatial_cell"):
        raise ValueError("Unknown source sampling-robustness weighting")
    if mode=="photo_weighted":
        sample=work
        item=work.groupby("accepted_key",sort=True)
        a=item.agg(image_mean=(ENDPOINT,"mean"),image_median=(ENDPOINT,"median"),
                   n_source_images=(ENDPOINT,"size"))
        means=item[list(EXPOSURES)+list(GEO_COLS)+
                   [f"sector_{n}" for n in SECTOR_NAMES]].mean()
        n_cells=work.groupby("accepted_key")[["_cell_lat","_cell_lon"]].apply(
            lambda x:len(x.drop_duplicates()))
        a["n_geo_cells"]=n_cells
    else:
        grouped=work.groupby(["accepted_key","_cell_lat","_cell_lon"],
                             as_index=False,sort=True)
        cell=grouped.agg(**{ENDPOINT:(ENDPOINT,"mean"),
                         **{c:(c,"mean") for c in list(EXPOSURES)+list(GEO_COLS)+
                            [f"sector_{n}" for n in SECTOR_NAMES]}})
        item=cell.groupby("accepted_key",sort=True)
        a=item.agg(image_mean=(ENDPOINT,"mean"),image_median=(ENDPOINT,"median"),
                   n_geo_cells=(ENDPOINT,"size"))
        a["n_source_images"]=work.groupby("accepted_key").size()
        means=item[list(EXPOSURES)+list(GEO_COLS)+
                   [f"sector_{n}" for n in SECTOR_NAMES]].mean()
    names=work.groupby("accepted_key").taxon_name.agg(
       lambda s:sorted(set(s.dropna().astype(str).str.strip())))
    output=a.join(means).join(names.rename("source_names")).reset_index()
    output["source_concept_collision"]=output.source_names.map(lambda x:len(x)!=1)
    props=output[[f"sector_{s}" for s in SECTOR_NAMES]].to_numpy(float)
    output["dominant_sector_id"]=np.argmax(props,axis=1).astype(int)
    output["dominant_sector_fraction"]=props.max(axis=1)
    if mode=="photo_weighted" and output.n_source_images.sum()!=len(work):
        raise ValueError("S1 per-taxon observation sum missing")
    return output


def tree_subset_keys():
    mapping=json.loads(MAPPING_PATH.read_text(encoding="utf-8"))
    matched={str(x["accepted_key"]) for x in mapping if x["status"]=="matched"
             and len(x["source_names"])==1 and len(x["tree_matches"])==1}
    if len(matched)!=191:
        raise ValueError("Frozen S1 strict191 source/tree matching drift")
    return matched


def select_cohorts(frame,mapped):
    clean=frame.loc[~frame.source_concept_collision].copy()
    tree=clean.loc[clean.accepted_key.isin(mapped)].copy()
    sampled=tree.loc[tree.n_source_images>=10].copy()
    if len(tree)!=191:
        raise ValueError("Original 191 concept-clean tree matched S1 taxa unavailable")
    if len(sampled)!=143:
        raise ValueError("Original 143 >=10-photo S1 taxon source unavailable")
    if len(clean)<220:
        raise ValueError("Expected >220 noncollided native S1 taxa")
    return {
      "clean_native":clean,
      "clean_tree":tree,
      "clean_tree_min10":sampled,
    }


def design(cohort,geo_control):
    rows=cohort.sort_values("accepted_key")
    env=rows[list(EXPOSURES)].to_numpy(float)
    mean=env.mean(axis=0)
    fixed=(env-mean)/SCALES
    spatial=[]
    if geo_control=="none":
        pass
    elif geo_control=="three_longitude_sector_proportions":
        spatial=rows[[f"sector_{name}" for name in SECTOR_NAMES[:2]]].to_numpy(float)
    elif geo_control=="spherical_second_order":
        spatial=rows[list(GEO_COLS)].to_numpy(float)
    else:
        raise ValueError("Unregistered spatial control")
    X=np.column_stack([np.ones(len(rows)),fixed])
    if len(spatial):
        sig=np.std(spatial,axis=0)
        if np.any(sig<=1e-12):
            raise ValueError("One or more spatial predictors constant in cohort")
        X=np.column_stack((X,(spatial-spatial.mean(axis=0))/sig))
    if np.linalg.matrix_rank(X)!=X.shape[1] or np.linalg.cond(X)>1e7:
        raise ValueError("Model design collinear for selected taxon geography")
    return X,rows


def fit_hc3(X,y):
    n,p=X.shape
    if n<=p+8:
        raise ValueError("Too few species in geographic comparison")
    solve=np.linalg.inv(X.T@X)
    beta=solve@(X.T@y)
    err=y-X@beta
    h=np.einsum("ij,jk,ik->i",X,solve,X)
    if np.any(h>=1-1e-8):
        raise ValueError("Unstable taxon-level leverage after geographic controls")
    score=X*(err/(1-h))[:,None]
    v=solve@(score.T@score)@solve
    se=np.sqrt(np.maximum(0,np.diag(v)))
    q=float(t.ppf(.975,n-p))
    interval=np.column_stack((beta-q*se,beta+q*se))
    return beta,interval,err,h


def model_report(group,trait,geography):
    try:
        X,rows=design(group,geography)
        outcome=rows[trait].to_numpy(float)
        beta,interval,err,h=fit_hc3(X,outcome)
    except (ValueError,np.linalg.LinAlgError) as e:
        return {"status":"not_identifiable","reason":str(e),
                "n_source_taxa":len(group),"spatial_control":geography}
    results={}
    for axis in PREDICTOR_NAMES:
        idx=list(EXPOSURES).index(axis)+1
        results[axis]={
           "estimate_image_degree_per_registered_increment":float(beta[idx]),
           "HC3_taxon_interval95":interval[idx].tolist(),
           "increment":float(SCALES[idx-1]),
           "positive_pointwise_HC3_interval":bool(interval[idx,0]>0),
           "not_spatially_cluster_bootstrap_or_causal_effect":True,
        }
    return {"status":"estimated","n_source_taxa":len(rows),
            "degrees_of_freedom":int(len(rows)-X.shape[1]),
            "predictor_rank":int(X.shape[1]),
            "spatial_control":geography,
            "n_source_name_collisions":int(rows.source_concept_collision.sum()),
            "source_study_area":"taxon means across observed native photographic records",
            "max_taxon_leverage":float(h.max()),
            "response":trait,
            "exposures":results}


def directional_sector_holdouts(group,trait):
    rows=group.sort_values("accepted_key")
    base=rows.dominant_sector_id.to_numpy(int)
    counts={SECTOR_NAMES[j]:int(np.sum(base==j)) for j in range(3)}
    result={}
    for j,sector in enumerate(SECTOR_NAMES):
        if counts[sector]<8:
            result[sector]={"status":"not_estimable_too_few_source_taxa",
                            "taxa_held_out":counts[sector]}
            continue
        kept=rows.loc[base!=j]
        if len(kept)<60:
            result[sector]={"status":"not_estimable_remaining_cohort_small",
                            "taxa_held_out":counts[sector]}
            continue
        report=model_report(kept,trait,"none")
        report["n_heldout_source_taxa"]=counts[sector]
        report["holdout_rule"]="fixed_dominant_longitude_sector_preselected_before_y_fit"
        result[sector]=report
    return {"all_original_source_taxa_sector_counts":counts,
            "fixed_predefined_sector_holdouts":result,
            "taxa_with_dominant_sector_share_below_0_7":
                int(sum(rows.dominant_sector_fraction < .7)),
            "not_previously_defined_terrestrial_biogeographic_realm":True}



def regional_within_sector_models(cohort,trait):
    """Direct within-sector S1 photo-angle climate contrasts and data support.

    These fixed longitude-sector aggregates represent native photo-geographic
    composition, not political continents nor independent within-species tests.
    Species are assigned by their dominant observed longitude sector.
    """
    rows=cohort.sort_values("accepted_key")
    geography=rows.dominant_sector_id.to_numpy(int)
    summary={}
    for number,label in enumerate(SECTOR_NAMES):
        sub=rows.loc[geography==number]
        if len(sub)<15:
            summary[label]={"status":"too_few_taxa",
                            "n_taxa":int(len(sub))}
            continue
        reg=model_report(sub,trait,"none")
        annual=sub["chelsa_bio12"].to_numpy(float)
        angle=sub[trait].to_numpy(float)
        summary[label]={
            "n_species":int(len(sub)),
            "n_source_images":int(sub.n_source_images.sum()),
            "n_species_primary_source_sector_concentration_gt70pct":
                 int(sum(sub.dominant_sector_fraction>=.7)),
            "BIO12_taxon_mean_yearly_mm_q10_q50_q90":
                np.quantile(annual,[.1,.5,.9]).tolist(),
            "image_angle_taxon_mean_degree_q10_q50_q90":
                np.quantile(angle,[.1,.5,.9]).tolist(),
            "taxon_equal_mean_original_image_angle_deg":float(angle.mean()),
            "taxon_equal_mean_BIO12_mm_year":float(annual.mean()),
            "within_sector_5_CHELSA_HC3_result":reg,
            "geographic_qualification":
                "fixed longitude sector only, not a biological realm, and no species-level evolutionary origin conclusion",
        }
    return {"sector_selection_predeclared_without_response":True,
            "regional_results":summary}


def execute(table,out):
    work,audit=authenticated_source(table)
    mapped=tree_subset_keys()
    models={}
    for weighting in ("photo_weighted","equal_5deg_spatial_cell"):
        roster=taxon_rosters(work,weighting)
        groups=select_cohorts(roster,mapped)
        variant={}
        for cohort_key,cohort in groups.items():
            tests={}
            for phenotype in ("image_mean","image_median"):
                fits={model:model_report(cohort,phenotype,model)
                      for model in GEO_ADJUSTMENTS}
                holds=directional_sector_holdouts(cohort,phenotype)
                within_region=regional_within_sector_models(cohort,phenotype)
                tests[phenotype]={"model_comparisons":fits,
                                   "sector_holdouts":holds,
                                   "within_sector_and_composition":within_region}
            variant[cohort_key]={
                "n_species":int(len(cohort)),
                "n_geographical_5degree_cells_total":int(cohort.n_geo_cells.sum()),
                "n_images":int(cohort.n_source_images.sum()),
                "n_taxa_with_majority_geographic_sector_lt70pct":
                    int(sum(cohort.dominant_sector_fraction<.7)),
                "response_and_geography_models":tests}
        models[weighting]=variant
    receipt={
       "schema":"azami_v4_S1_original_spatial_rainfall_photo_angles_five_CHELSA_cohort_20261022",
       "status":"EMPIRICAL_ORIGINAL_S1_GEO_SPATIAL_RAINFALL_SENSITIVITY",
       "original_S1_input_sha256":ORIGINAL_SHA,
       "original_S1_fixed_five_climate":list(EXPOSURES),
       "S1_source_audit":audit,
       "geography_base_lon_sectors":list(SECTOR_NAMES),
       "spherical_second_order_basis":list(GEO_COLS),
       "analysis_groups":models,
       "scope":{
         "same_continuous_image_axis_as_original_S1":True,
         "native_only_original_S1_source":True,
         "photographic_vertical_not_gravity":True,
         "equal_taxon_weighting_but_photo_or_spatial_cell_within_taxon":True,
         "climate_seasonality_is_BIO15_not_dry_season_extremes":True,
         "spatial_spherical_and_longitude_controls_are_sensitivity_not_causation":True,
         "original_S1_between_within_and_V3_frozen_not_refitted":True,
         "no_independent_evolutionary_nodding_gain_count_or_wetting_fitness_inferred":True,
         "never_publish_exact_observation_coordinates":True,
       },
    }
    dst=Path(out);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(receipt,indent=2,sort_keys=True,
                               ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    wanted={}
    for sample in ("photo_weighted","equal_5deg_spatial_cell"):
      for cohort in ("clean_native","clean_tree","clean_tree_min10"):
        c=models[sample][cohort]
        wanted[f"{sample}|{cohort}"]={
          "n_species":c["n_species"],
          "image_mean":c["response_and_geography_models"]["image_mean"]}
    print("S1_RAIN_SPATIAL_FIVE_COVAR_SUMMARY="+json.dumps({
       "status":receipt["status"],"source":audit,"panels":wanted,
       "scope":receipt["scope"]},sort_keys=True,ensure_ascii=False),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table",required=True)
    p.add_argument("--out",required=True)
    options=p.parse_args()
    execute(options.table,options.out)


if __name__=="__main__":
    main()
