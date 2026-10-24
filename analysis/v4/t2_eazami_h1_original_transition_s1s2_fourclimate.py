"""Exact EAzami U->D CTMC/Brownian transition workflow on NEW S1/S2 four axes.

Reuse EAzami 2026-09-04 frozen 6 real Comp1061 nuclear AU topologies,
original FDT4 authentic Japanese/Taiwanese occurrence data, historical
flora-derived botanical U/D taxon-concept state crosswalk, exact 792 and
1716 count-preserving state-map null. New exposure data from the SAME public
CHELSA BIO18/BIO05/RSDS and WorldClim 30" elevation rasters as azami S1/S2.
The EAzami H1 helper is SHA-pin verified and imported rather than replaced
with PIC or PGLS. A separate frozen original BIO1/BIO15 n10 numeric
reproduction must pass before interpreting the four NEW axes.

New S1-selected directional prediction is U->D vs HIGHER BIO18; the other
three axes are two-sided/omnibus diagnostics, not manufactured hypotheses.
This is post-result exploratory cross-study translation, NOT a new originally
preregistered EAzami H1 or independent recoding of photographic S1 angles.
The original taxon-authority U/D status is NEVER an original photograph
measurement, and branchwise CTMC weights are not counts of independent gains.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from analysis.v4 import s1_s2_real_bio18_bio05_rsds_elevation as four
from analysis.v4.native_cached_climate import verified_cache, sample_cached_pixels

EA_HEAD_COMMIT="c9870e91604ed1514344ddc8f765eae5237dbf27"
EA_H1_SOURCE_BLOB="09731b563907b3ec2778ea04086d729771376f7e"
EA_SIX_TREE_BLOB="d7f8c0b96ee9c257c815db06a5a0aa9fdae4841e"
EA_UD_BLOB="f655544e77f1659323c1e915567dacbb51b01cc3"
EA_CONTRACT_BLOB="59acc74cbcace6da78226c62fc976c0e543f0225"
EA_COVERAGE_BLOB="47e97cb35d241c0dcba8cb9c844038b568db03d5"
EA_LEGACY_BLOB="5cb259ee373b33c1f16032b3ff9bbfd699913826"
EA_REFERENCE_BLOB="fb7992e683301ddf307be3ab7be35a2f246dddde"
EA_JP_ARTIFACT_ZIP_SHA="cb5717c90cbc769398f9ffce38787704a02e34253c1ff16cf8f368ba55aeda96"
EA_TW_ARTIFACT_ZIP_SHA="819936ceeb6008f356c580622fb72972617700402afaa3d3ab125bff2e6d9106"
RSDS_ORIGINAL_SHA="cc1ad14a70a1a940be3d01d918d4c3f7e35db95b86b12af9be5487e1bbc1fc7d"
AXES=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
PANELS={
 "n5_primary":("5",12,7,5,792),
 "n3_sensitivity":("3",13,7,6,1716),
 "n10_strict_sensitivity":("10",9,5,4,126),
}
MAX_CENSUS=748
MIN_COVERAGE=0.90
NEGATIVE_NICHE_CAUSALITY=False

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def checked_blob(file,expected,name):
    data=Path(file).read_bytes()
    if git_blob(data)!=expected:
        raise ValueError(f"Frozen {name} EAzami source Git blob changed")
    return data

def load_original_python(path):
    checked_blob(path,EA_H1_SOURCE_BLOB,"EAzami original CTMC/Brownian implementation")
    spec=importlib.util.spec_from_file_location("eazami_original_transition_H1",path)
    if spec is None or spec.loader is None:
        raise ValueError("EAzami source module unavailable")
    h1=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h1)
    for name in ("read_trees","zscore","fit_symmetric_q",
                 "ctmc_likelihood_and_edge_joint","brownian_internal_values",
                 "tree_structure","panel_state_map","prepare_topology_assets",
                 "observed_topology_stats","method_check","exact_panel_test",
                 "normalize_tip","prune_to_panel"):
        if not hasattr(h1,name):
            raise ValueError("Frozen EAzami original H1 implementation missing "+name)
    return h1

def frozen(h1,files):
    checked_blob(files["trees"],EA_SIX_TREE_BLOB,"6 Comp1061 nuclear AU trees")
    checked_blob(files["orientation"],EA_UD_BLOB,"taxonomic original U/D")
    checked_blob(files["contract"],EA_CONTRACT_BLOB,"original H1 hypothesis")
    checked_blob(files["coverage"],EA_COVERAGE_BLOB,"original species coverage")
    checked_blob(files["legacy"],EA_LEGACY_BLOB,"legacy method benchmark")
    checked_blob(files["reference"],EA_REFERENCE_BLOB,"original H1 saved results")
    trees=h1.read_trees(files["trees"],6)
    orientation=pd.read_csv(files["orientation"])
    if orientation.accepted_taxon.duplicated().any():
        raise ValueError("Original EAzami U/D crosswalk duplicate")
    coverage=json.loads(Path(files["coverage"]).read_text())
    contract=json.loads(Path(files["contract"]).read_text())
    legacy=json.loads(Path(files["legacy"]).read_text())
    reference=json.loads(Path(files["reference"]).read_text())
    if (coverage["version"]!="chapter2_orientation_occurrence_coverage_audit_result_v1"
        or contract["version"]!="chapter2_orientation_transition_regime_hypothesis_contract_v1"
        or legacy["contract_version"]!="fdt4_branchwise_niche_transition_concordance_v1"):
        raise ValueError("Frozen original EAzami H1 contracts no longer match")
    for name,(gate,n,nu,nd,nmaps) in PANELS.items():
        taxa=coverage["threshold_summaries"][gate]["taxa"]
        states=orientation.set_index("accepted_taxon").loc[taxa].analysis_state
        if (len(taxa),int(sum(states=="U")),int(sum(states=="D"))) != (n,nu,nd):
            raise ValueError("EAzami true U/D source-gated species census changed "+name)
        if math.comb(n,nd)!=nmaps:
            raise ValueError("EAzami count-preserving null size changed")
    return trees,orientation,coverage,contract,legacy,reference

def original_occurrences(jp,tw,h1,trees,orientation,coverage,contract,legacy,reference):
    jp=pd.read_csv(jp,low_memory=False)
    tw=pd.read_csv(tw,low_memory=False)
    expected={"Japan":693,"Taiwan":55}
    if len(jp)!=expected["Japan"] or len(tw)!=expected["Taiwan"]:
        raise ValueError("Frozen original EAzami environment-complete occurrences drifted")
    jp["_region"]="Japan"
    tw["_region"]="Taiwan"
    old=pd.concat([jp,tw],ignore_index=True)
    required={"scientific_name_query","scientificName","latitude","longitude",
        "chelsa_bio01","chelsa_bio15","environment_complete"}
    if not required<=set(old):
        raise ValueError("EAzami source occurrence or original climate fields changed")
    old["latitude"]=pd.to_numeric(old.latitude,errors="raise")
    old["longitude"]=pd.to_numeric(old.longitude,errors="raise")
    if (not np.isfinite(old[["latitude","longitude"]]).all().all()
        or not old.latitude.between(-90,90).all() or not old.longitude.between(-180,180).all()):
        raise ValueError("Original EAzami FDT4 coordinates invalid")
    if len(set(old.scientific_name_query))!=14:
        raise ValueError("EAzami original 14-taxonomic concept public occurrence source drifted")
    # Faithfully validate author's own H1 branch score from original BIO1/BIO15.
    strict=coverage["threshold_summaries"]["10"]["taxa"]
    counts,old_env=h1.build_panel_environment(old,strict)
    states=h1.panel_state_map(orientation,strict)
    assets=h1.prepare_topology_assets(trees,strict,old_env)
    result=h1.method_check(contract,legacy,h1.observed_topology_stats(assets,states))
    if result["status"]!="pass":
        raise ValueError("EAzami actual original 9 species numeric H1 reproduction FAILED")
    # Exact n5 observed old composite on six true AU topologies must remain.
    taxa=coverage["threshold_summaries"]["5"]["taxa"]
    _,env=h1.build_panel_environment(old,taxa)
    original_scores=h1.observed_topology_stats(
        h1.prepare_topology_assets(trees,taxa,env),
        h1.panel_state_map(orientation,taxa))
    original_median=float(np.median([x["composite"] for x in original_scores]))
    expected=float(reference["n5_primary"]["observed_composite_median"])
    if abs(original_median-expected)>0.001:
        raise ValueError("EAzami original BIO1/BIO15 H1 n5 composite cannot be reproduced on 6 true AU trees")
    return old,{"status":"ORIGINAL_EAZAMI_FROZEN_CTMC_BROWNIAN_REPRODUCED",
       "n_original_occurrences":len(old),"n_Japan":len(jp),"n_Taiwan":len(tw),
       "source_taxa":len(set(old.scientific_name_query)),
       "old_n10_numeric_method_check":result,
       "old_n5_six_AU_observed_composite_median":original_median,
       "old_n5_expected_frozen_composite_median":expected,
       "old_BIO15up_BIO1down_primary_exact_rank":"16/792",
       "old_method_different_hypothesis_than_new_four_axes":True}

def sample_new_climate(df,raster_dir,rsds_dir,rsds_meta):
    latitude=df.latitude.to_numpy(float)
    longitude=df.longitude.to_numpy(float)
    site=df.copy()
    evidence={}
    for axis,file in (("chelsa_bio18","CHELSA_bio18_1981-2010_V.2.1.tif"),
                      ("chelsa_bio05","CHELSA_bio05_1981-2010_V.2.1.tif"),
                      ("elevation_m","wc2.1_30s_elev.tif")):
        result,meta=four.raster_sample(Path(raster_dir)/file,latitude,longitude)
        site[axis]=result
        evidence[axis]=meta
    meta=json.loads(Path(rsds_meta).read_text())["chelsa_rsds_mean"]
    if meta["source_sha256"]!=RSDS_ORIGINAL_SHA:
        raise ValueError("Original S1/S2 raster RSDS provenance changed")
    location=Path(rsds_dir)/"chelsa_rsds_mean"
    grid=verified_cache(location,meta)
    readings=sample_cached_pixels(location,grid,latitude,longitude)
    site["chelsa_rsds_mean"]=readings*float(grid["scale"])+float(grid["offset"])
    evidence["chelsa_rsds_mean"]={
        "source_sha256":grid["source_sha256"],
        "source_mode":"original_S1_S2_316_SHA_verified_tiles_pixel_values",
        "n_matched":int(site["chelsa_rsds_mean"].notna().sum())}
    complete=site[list(AXES)].notna().all(axis=1)
    if float(complete.mean())<MIN_COVERAGE:
        raise ValueError("Original EAzami species environment lacks 90pct exact 4 climate matching")
    keep=site.loc[complete].copy()
    if len(keep)<680:
        raise ValueError("Below original native/GBIF environment record coverage threshold")
    for field,low,high in (("chelsa_bio18",0,10000),("chelsa_bio05",-70,70),
                           ("chelsa_rsds_mean",0,50),("elevation_m",-500,9000)):
        median=float(keep[field].median())
        if not low<median<high:
            raise ValueError("New physical climate unit suspicious: "+field+" median "+str(median))
    return keep,{
        "original_occurrences_before_four_axis":len(site),
        "n_observations_four_axis_complete":len(keep),
        "four_axis_complete_fraction":float(complete.mean()),
        "n_source_species_after_complete":int(keep.scientific_name_query.nunique()),
        "by_species_old_vs_new":{name:{
           "old_n":int(len(allrows)),
           "new_four_axis_n":int((keep.scientific_name_query==name).sum())}
             for name,allrows in site.groupby("scientific_name_query")},
        "raster_sampling_evidence":evidence,
        "source_identity":"EAzami_original_geo_thinned_Japan693_Taiwan55_FDT4_occurrences",
        "environment_source_same_as_S1_S2_BIO18_BIO05_RSDS_elevation":True,
    }

def edge_assets(h1,trees,taxa,env):
    """Identical Brownian reconstruction as EAzami; one Δ per actual nuclear edge."""
    tip_z={axis:{h1.normalize_tip(t):float(z) for t,z in
        zip(taxa,h1.zscore(env.loc[taxa,axis].to_numpy(float)))}
        for axis in AXES}
    result=[]
    for tree in trees:
        tr=h1.prune_to_panel(tree,taxa)
        source_edges=h1.tree_structure(tr)[2]
        keys=sorted(source_edges)
        vals=[]
        for axis in AXES:
            _,delta=h1.brownian_internal_values(tr,tip_z[axis])
            if set(delta)!=set(keys):
                raise ValueError("EAzami Brownian reconstructed topology edge drift")
            vals.append([float(delta[c]) for c in keys])
        result.append({"tree":tr,"edge_keys":keys,
                       "deltas":np.asarray(vals,float).T})
    return result

def transition_weighted_scores(h1,assets,state_map):
    """Calculate EXACT EAzami edge-joint U->D minus D->U weighted shifts."""
    output=[];changes=[]
    for item in assets:
        tr=item["tree"]
        q=h1.fit_symmetric_q(tr,state_map)
        Z,joint=h1.ctmc_likelihood_and_edge_joint(tr,state_map,q,need_joint=True)
        if not Z>0 or joint is None:
            raise ValueError("EAzami two-state CTMC likelihood invalid")
        weights=np.array([joint[e][0,1]-joint[e][1,0]
                          for e in item["edge_keys"]],float)
        transitions=np.array([joint[e][0,1]+joint[e][1,0]
                               for e in item["edge_keys"]],float)
        if not np.isfinite(weights).all() or transitions.sum()<1e-12:
            raise ValueError("Original CTMC expected transition mass invalid")
        output.append(weights@item["deltas"]/transitions.sum())
        changes.append(float(transitions.sum()))
    return np.asarray(output,float),changes

def tails_null(scores,obs,alternative):
    """Exact finite-map ranks, not independent replicate P values."""
    null=np.asarray(scores,float)
    return {
      "higher_tail_count":int(np.count_nonzero(null>=obs-1e-12)),
      "lower_tail_count":int(np.count_nonzero(null<=obs+1e-12)),
      "n_null_state_maps":len(null),
      "higher_tail_fraction":float(np.mean(null>=obs-1e-12)),
      "lower_tail_fraction":float(np.mean(null<=obs+1e-12)),
      "two_sided_finite_map_fraction":float(min(1,2*min(
         np.mean(null>=obs-1e-12),np.mean(null<=obs+1e-12)))),
      "prespecified_tail":alternative,
      "one_sided_directional_fraction":float(np.mean(null>=obs-1e-12)) if alternative=="positive"
          else None,
      "finite_permutation_exact_rank_not_independent_biological_p_value":True}

def map_rank(h1,trees,taxa,states,env,*,region_by_taxon=None,all_state_maps=True):
    if not all_state_maps and not region_by_taxon:
        raise ValueError("Region-conditioned null needs original regional names")
    assets=edge_assets(h1,trees,taxa,env)
    ids=[h1.normalize_tip(t) for t in taxa]
    observed,expchanges=transition_weighted_scores(h1,assets,states)
    obsmed=np.median(observed,axis=0)
    d_count=int(sum(states.values()))
    if all_state_maps:
        assignments=itertools.combinations(range(len(taxa)),d_count)
        expected=math.comb(len(taxa),d_count)
        label="all_tip_U_D_count_preserving_exact"
    else:
        regions={r:[i for i,t in enumerate(taxa) if region_by_taxon[t]==r]
                 for r in sorted(set(region_by_taxon[t] for t in taxa))}
        choices=[]
        for region,index in regions.items():
            need=sum(states[ids[i]] for i in index)
            choices.append([c for c in itertools.combinations(index,need)])
        assignments=(tuple(i for chunk in group for i in chunk)
                     for group in itertools.product(*choices))
        expected=int(np.prod([len(q) for q in choices]))
        label="source_region_specific_U_D_count_preserving_exact"
    null=np.empty((expected,len(AXES)),dtype=float)
    same=0
    for it,assignment in enumerate(assignments):
        subset=set(assignment)
        proxy={n:int(i in subset) for i,n in enumerate(ids)}
        if proxy==states:same+=1
        stat,_=transition_weighted_scores(h1,assets,proxy)
        null[it]=np.median(stat,axis=0)
    if len(null)!=expected or same!=1:
        raise ValueError("Original finite EAzami state map enumeration not exhaustive")
    signed_ref={k:float(v) for k,v in zip(AXES,obsmed)}
    summaries={}
    for i,axis in enumerate(AXES):
        more=tails_null(null[:,i],float(obsmed[i]),"positive" if axis=="chelsa_bio18" else "two_sided")
        more.update({
         "median_observed_signed_transition_weighted_branch_niche_shift":float(obsmed[i]),
         "six_nuclear_topology_observed_signed_effects":observed[:,i].tolist(),
         "n_positive_six_topologies":int(np.sum(observed[:,i]>0)),
         "null_median":float(np.median(null[:,i])),
         "null_95_interval":[float(v) for v in np.quantile(null[:,i],[.025,.975])],
        })
        summaries[axis]=more
    # Null calibration preserves correlation among all 4 climate axes.
    scales=null.std(axis=0,ddof=0)
    if np.any(scales<1e-12):
        raise ValueError("Source climate null axis has degenerate true nuclear contrast")
    score=np.max(np.abs((null-null.mean(0))/scales),axis=1)
    orig=float(np.max(np.abs((obsmed-null.mean(0))/scales)))
    maxT=float(np.mean(score>=orig-1e-12))
    for i,axis in enumerate(AXES):
        point=float(abs((obsmed[i]-null[:,i].mean())/scales[i]))
        summaries[axis]["four_axis_two_sided_maxT_FWER_fraction"]=float(
            np.mean(score>=point-1e-12))
    omnibus2=np.sum(((null-null.mean(0))/scales)**2,axis=1)
    v_obs=float(np.sum(((obsmed-null.mean(0))/scales)**2))
    return {
      "status":"EAzami_IDENTICAL_SYMMETRIC_CTMC_BROWNIAN_EDGE_POSTERIOR_MAP_TEST_4AXES",
      "n_original_taxa":len(taxa),
      "observed_U":len(taxa)-d_count,"observed_D":d_count,
      "n_accepted_comp1061_nuclear_topologies":len(trees),
      "n_exact_state_maps":expected,
      "source_taxa":list(taxa),
      "null_condition":label,
      "number_observed_state_maps":same,
      "original_EAzami_expected_CTMC_change_count_per_topology":expchanges,
      "axis_effects":summaries,
      "BIO18_positive_one_sided_observed":bool(observed[:,0].min()>0),
      "BIO18_positive_one_sided_tail_fraction":summaries["chelsa_bio18"]["higher_tail_fraction"],
      "four_axis_global_maxT_abs_observed":orig,
      "four_axis_global_maxT_two_sided_FWER_fraction":maxT,
      "four_axis_omnibus_squared_abs_z_score":v_obs,
      "four_axis_omnibus_finite_map_fraction":float(np.mean(omnibus2>=v_obs-1e-12)),
      "alternative_exploratory_not_frozen_H1":True,
      "not_convergence_or_adaptive_nodding_evolution":True,
      "no_true_transition_time_climate_observed":True,
    }

def original_H1_matched_new_four_sites(h1,trees,coverage,orientation,environment,reference):
    """The same AUTHOR CTMC/Brownian BIO15↑-BIO1↓ exact rank on exact NEW four-complete sites.

    Preserves 12 original taxonomic states / 6 original nuclear topologies /
    792 state maps, but source environmental niche centroids are computed from
    the same 689 occurrence sites retained by the four new climate layers.
    This is an attrition sensitivity, not a new independent validation.
    """
    taxa=coverage["threshold_summaries"]["5"]["taxa"]
    counts,env=h1.build_panel_environment(environment,taxa)
    if len(taxa)!=12 or any(int(counts.get(t,0))<5 for t in taxa):
        raise ValueError("New four-climate complete observations cannot support original EAzami n5 panel")
    states=h1.panel_state_map(orientation,taxa)
    assets=h1.prepare_topology_assets(trees,taxa,env)
    original_new,full_maps=h1.exact_panel_test(
       "matched_EAzami_H1_BIO15up_BIO1down_12_original_taxa",
       taxa,states,assets,792)
    observed=float(original_new["observed"]["composite_median"])
    old=float(reference["n5_primary"]["observed_composite_median"])
    new_rank=original_new["exact_primary_rank"]
    if new_rank["n_maps"]!=792 or len(full_maps)!=792:
        raise ValueError("EAzami original matched climate exact state-map rank incomplete")
    old_rank=reference["n5_primary"]["exact_primary_rank"]
    return {
       "status":"EAzami_ORIGINAL_BIO15up_BIO1down_H1_RECOMPUTED_ON_EXACT_4CLIMATE_COMPLETE_SITES",
       "original_taxa":12,
       "original_n_Japan_Taiwan":int(sum(counts[t] for t in taxa)),
       "n_original_4climate_complete_sites_entire_14_taxon_universe":len(environment),
       "six_frozen_nuclear_AU_trees":len(trees),
       "state_preserving_null_12_choose_5":792,
       "original_frozen_748site_BIO15up_minus_BIO1down_composite_median":old,
       "original_frozen_748site_rank":{"n_maps":int(old_rank["n_maps"]),"count_at_least_observed":int(old_rank["count_at_least_observed"])},
       "matched_4climate_site_BIO15up_minus_BIO1down_composite_median":observed,
       "matched_fourclimate_exact_rank":new_rank,
       "matched_fourclimate_all_six_topology_scores":original_new["observed"]["topology_composite"],
       "matched_old_positive_all_six":all(float(t)>0 for t in original_new["observed"]["topology_composite"]),
       "same_site_old_H1_original_threshold_passes":
          bool(all(float(t)>0 for t in original_new["observed"]["topology_composite"]) and
               new_rank["exact_fraction"]<=.05),
       "old_data_are_biological_states_not_camera_head_angle":True,
       "new_source_site_attrition_can_explain_part_of_cross_variable_difference":True,
       "not_independent_cross_study_confirmation":True,
    }


def sensitivity(h1,trees,original,orientation,environment,region):
    # n5 maps exactly original source, n3 and n10 coverage sensitivity.
    output={}
    n5=original["threshold_summaries"]["5"]["taxa"]
    for panel,(gate,n,nu,nd,maps) in PANELS.items():
        taxa=original["threshold_summaries"][gate]["taxa"]
        count=environment.scientific_name_query.value_counts()
        below={t:int(count.get(t,0)) for t in taxa if count.get(t,0)<int(gate)}
        if below:raise ValueError("Original EAzami taxon inclusion requires original occurrence support: "+str(below))
        state=h1.panel_state_map(orientation,taxa)
        env=environment.groupby("scientific_name_query")[list(AXES)].mean().loc[taxa]
        main=map_rank(h1,trees,taxa,state,env)
        if main["n_exact_state_maps"]!=maps:
            raise ValueError("Finite H1 map census changed after new axes")
        rawcounts={t:int(count[t]) for t in taxa}
        main["new_four_axis_source_counts_by_frozen_taxon"]=rawcounts
        main["old_EAzami_n5_primary_or_sensitivity"]=panel
        output[panel]=main
    # The EAzami regional null is a fixed separate stress, not a new discovery.
    taxa=n5
    state=h1.panel_state_map(orientation,taxa)
    env=environment.groupby("scientific_name_query")[list(AXES)].mean().loc[taxa]
    region_map={t:region[t] for t in taxa}
    if list(region_map.values()).count("Japan")!=8 or list(region_map.values()).count("Taiwan")!=4:
        raise ValueError("Original EAzami n5 Japan8/Taiwan4 stratification changed")
    check=map_rank(h1,trees,taxa,state,env,region_by_taxon=region_map,all_state_maps=False)
    if check["n_exact_state_maps"]!=336:
        raise ValueError("Expected source-preserving 336 Japan/Taiwan state placements")
    output["n5_region_frequency_preserving_stress"]=check
    return output

def execute(h1_path,trees_path,orientation_path,coverage_path,contract_path,
             legacy_path,reference_path,jp_path,tw_path,
             rasters,rsds,rsds_meta,outdir):
    h1=load_original_python(h1_path)
    paths={"trees":trees_path,"orientation":orientation_path,"coverage":coverage_path,
      "contract":contract_path,"legacy":legacy_path,"reference":reference_path}
    trees,orientation,coverage,contract,legacy,reference=frozen(h1,paths)
    original,method=original_occurrences(jp_path,tw_path,h1,trees,orientation,
                                         coverage,contract,legacy,reference)
    sampled,qc=sample_new_climate(original,rasters,rsds,rsds_meta)
    matched_old=original_H1_matched_new_four_sites(h1,trees,coverage,orientation,sampled,reference)
    region={t:group["_region"].iloc[0] for t,group in
        original.groupby("scientific_name_query")}
    if any(group._region.nunique()!=1 for _,group in original.groupby("scientific_name_query")):
        raise ValueError("Original Japan/Taiwan source region disjointness changed")
    result=sensitivity(h1,trees,coverage,orientation,sampled,region)
    r={
      "schema":"azami_v4_EAzami_exact_U_to_D_CTMC_Brownian_four_original_S1_S2_exposures_20261024_v1",
      "status":"COMPLETED_EAzami_FROZEN_BINARY_UD_CTMC_WITH_REAL_NEW_CLIMATES_NOT_S1_PHOTO_ANGLE",
      "original_EAzami_frozen_H1_method":method,
      "original_H1_BIO15up_BIO1down_recomputed_on_same_689_sites":matched_old,
      "original_EAzami_method_Git_blob":EA_H1_SOURCE_BLOB,
      "original_EAzami_source_commit":EA_HEAD_COMMIT,
      "original_six_accepted_nuclear_AU_topologies_from_EAzami_blob":EA_SIX_TREE_BLOB,
      "original_frozen_taxon_authority_U_D_from_EAzami_blob":EA_UD_BLOB,
      "original_continuous_photo_angle_not_used_as_U_D":True,
      "new_S1_S2_four_climate_exposures":list(AXES),
      "new_EAzami_original_photo_independent_GBIF_occurrence_pixel_QC":qc,
      "frozen_panels_exact_maps":result,
      "post_hoc_primary_direction_based_on_S1_between_species":
          "U_to_D_weighted_edge_climate_change toward greater BIO18 alone; NOT originally preregistered EAzami H1.",
      "other_BIO5_RSDS_elevation_axes_exploratory_two_sided":True,
      "all_four_joint_permutation_maxT_corrects_family_of_four":True,
      "cannot_combine_older_frozen_original_H1_significance_with_new_four_axes_as_two_independent_experiments":True,
      "no_S1_photo_angle_effect_estimated_here":True,
      "no_Historical_transition_time_environment_or_causal_rain_mechanism":True,
      "no_true_gravity_referenced_photographs_inferred":True,
      "S1_S2_original_4axis_source_reused_at_exact_EAzami_GBIF_survey_sites":True,
      "source_observation_coordinates_NOT_in_public_outputs":True,
    }
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    (outdir/"EAzami_S1S2_BIO18_BIO05_RSDS_ELEV_CTMC_BRANCHWISE_RESULTS.json").write_text(
      json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n")
    # Leave a taxon-level climate sample-count receipt only, never source sites.
    cols=["original_taxon","source_region","n_original_EAzami_sites",
          "n_new_four_axis_complete_EAzami_sites",*AXES,"original_U_D_authority"]
    taxon=[]
    state_by=dict(zip(orientation.accepted_taxon,orientation.analysis_state))
    counts=qc["by_species_old_vs_new"]
    for name in sorted(set(original.scientific_name_query)):
        rows=sampled.loc[sampled.scientific_name_query==name]
        taxon.append({
          "original_taxon":name,"source_region":region[name],
          "n_original_EAzami_sites":counts[name]["old_n"],
          "n_new_four_axis_complete_EAzami_sites":counts[name]["new_four_axis_n"],
          **{k:float(rows[k].mean()) if len(rows) else "" for k in AXES},
          "original_U_D_authority":state_by.get(name,"unresolved")
        })
    with (outdir/"EAzami_ORIGINAL_TAXON_CENTROIDS_S1S2_4AXIS.csv").open("w",newline="",encoding="utf-8") as fd:
        writer=csv.DictWriter(fd,fieldnames=cols)
        writer.writeheader();writer.writerows(taxon)
    print("EAzami_FROZEN_TRANSITION_NEW_S1_S2_FOUR_AXES="+json.dumps({
       "original_method_numeric_check":method["status"],
       "original_H1_on_four_complete_matched_sites_rank":matched_old["matched_fourclimate_exact_rank"],
       "original_H1_on_four_complete_matched_sites_supported":matched_old["same_site_old_H1_original_threshold_passes"],
       "original_12_taxa_BIO18_positive_tail":result["n5_primary"]["BIO18_positive_one_sided_tail_fraction"],
       "original_12_taxa_BIO18_six_topology_scores":result["n5_primary"]["axis_effects"]["chelsa_bio18"]["six_nuclear_topology_observed_signed_effects"],
       "original_12_taxa_all_four_omnibus_rank":result["n5_primary"]["four_axis_omnibus_finite_map_fraction"],
       "original_12_taxa_region_frequency_BIO18_tail":result["n5_region_frequency_preserving_stress"]["BIO18_positive_one_sided_tail_fraction"],
       "n5_exact_state_maps":result["n5_primary"]["n_exact_state_maps"],
       "n3_exact_state_maps":result["n3_sensitivity"]["n_exact_state_maps"]
    },sort_keys=True),flush=True)
    return r

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("h1","trees","orientation","coverage","contract","legacy","reference",
                 "jp","tw","rasters","rsds","rsds-meta","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    execute(a.h1,a.trees,a.orientation,a.coverage,a.contract,a.legacy,a.reference,
            a.jp,a.tw,a.rasters,a.rsds,a.rsds_meta,a.out)
if __name__=="__main__":main()
