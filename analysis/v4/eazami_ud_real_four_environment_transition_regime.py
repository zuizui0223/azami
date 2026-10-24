"""Re-execute EAzami's EXACT U->D CTMC transition-regime method on S1/S2's 4 environments.

THIS IS THE EAzami METHOD, NOT PGLS OR PIC:
 - Fixed East-Asian taxon concepts, flora U/D labels, EAzami true multilocus
   Comp1061 first six AU-non-rejected nuclear phylogenies.
 - Original Japan/Taiwan exact frozen GBIF occurrence rows restored untouched.
 - First replicate the ORIGINAL EAzami BIO15↑ and BIO1↓ branch statistics on
   strict n>=10 panel with original author's functions and original legacy bounds.
 - Newly sample *these same occurrence sites*, not azami's photographic sites,
   for CHELSA v2.1 BIO18, BIO05, original CHELSA RSDS cache, WorldClim DEM.
 - Compute on each of 6 nuclear trees actual symmetric 2-state CTMC ML q,
   exact-edge P(U->D)-P(D->U), and Brownian squared-change edge increments.
 - Enumerate every U/D count-preserving configuration for the original fixed
   n>=5 (12, choose(12,5)=792) PRIMARY, n>=10 (9, choose(9,4)=126)
   and n>=3 (13, choose(13,6)=1716) SENSITIVITY cohorts.
 - BIO18-positive S1-derived directional hypothesis and direction-free 4-axis
   omnibus conditional state-map null are both declared before reading new
   transition effects. BIO5, RSDS and elevation secondary two-sided checks.

Important: Unlike S1's own measured continuous image-frame angle, this EAzami
layer is literature-characterized botanical U/D. It MUST NOT be relabelled as
a PHOTO measurement or evidence of adaptation, dated environmental exposure,
or genuinely independent repeated transitions. This is a post-S1 exploratory
extension of an EXISTING EAzami transition workflow, not a prospective test.

Raw point observations, images, donor metadata, user IDs, lat/lon, and state-map
taxon strings are NOT published. Only public source-linked taxon aggregates and
summary statistics, without source occurrence IDs, are exported.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
from analysis.v4.native_cached_climate import verified_cache, sample_cached_pixels
from analysis.v4.s1_s2_real_bio18_bio05_rsds_elevation import raster_sample

EA_REPO_SHA="c9870e91604ed1514344ddc8f765eae5237dbf27"
EA_CORE_GIT_BLOB="09731b563907b3ec2778ea04086d729771376f7e"
EA_TREES_BLOB="d7f8c0b96ee9c257c815db06a5a0aa9fdae4841e"
EA_STATES_BLOB="f655544e77f1659323c1e915567dacbb51b01cc3"
EA_COVERAGE_BLOB="47e97cb35d241c0dcba8cb9c844038b568db03d5"
EA_CONTRACT_BLOB="59acc74cbcace6da78226c62fc976c0e543f0225"
EA_JP_CSV_SHA="d115635f57e11d8a9f5c417bf28992a3d2a00f2a4d02ebb710c96d892ad0a9b2"
EA_TW_CSV_SHA="63b43b670f55a16b5557e414e37297a9ee18942f06714094fc1b5788ec07c70c"
B18_TIF_SHA="68b0b1cb905cdf6972364e947e137f714a68d287f38863b1516161a2a4bbd7ec"
B05_TIF_SHA="90078d99eeb0b08dd1a115a5ffa60b7c19e5849686854db950283cbac9d3817d"
ELEV_TIF_SHA="1d57ce9fd3584cea6a227396bc75ef88522fd5346e8487c17370885c0499716e"
RSDS_SOURCE_SHA="cc1ad14a70a1a940be3d01d918d4c3f7e35db95b86b12af9be5487e1bbc1fc7d"
AXES=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
PANELS={
 "n5_primary":("5",12,7,5,792),
 "n10_strict_sensitivity":("10",9,5,4,126),
 "n3_resolution_sensitivity":("3",13,7,6,1716),
}
MISSING_STOP_RATE=0.03
MIN_VALID_PER_ORIGINAL_TAXON=0.75
EPS=1e-12

def sha(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()

def blob(raw:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def checked(path,expected_gitblob):
    raw=Path(path).read_bytes()
    if blob(raw)!=expected_gitblob:
        raise ValueError("Frozen EAzami source Git blob altered: "+str(path))
    return raw

def original_implementation(path):
    checked(path,EA_CORE_GIT_BLOB)
    spec=importlib.util.spec_from_file_location("exact_eazami_frozen_UD_CTMC",path)
    lib=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lib)
    required=["read_trees","prune_to_panel","brownian_internal_values",
              "panel_state_map","method_check","observed_topology_stats",
              "fit_symmetric_q","ctmc_likelihood_and_edge_joint","zscore",
              "tree_structure","normalize_tip","prepare_topology_assets"]
    if not all(callable(getattr(lib,name,None)) for name in required):
        raise ValueError("Missing function in EAzami author implementation")
    return lib

def original_occurrence_cohort(japan,taiwan,states,coverage):
    jp=Path(japan).read_bytes();tw=Path(taiwan).read_bytes()
    if sha(jp)!=EA_JP_CSV_SHA or sha(tw)!=EA_TW_CSV_SHA:
        raise ValueError("Original Japan/Taiwan frozen GBIF occurrence bytes altered")
    a=pd.read_csv(japan,low_memory=False)
    b=pd.read_csv(taiwan,low_memory=False)
    if (len(a),len(b))!=(693,55):
        raise ValueError("Original EAzami GBIF original 748 occurrence cohort changed")
    source=pd.concat([a,b],ignore_index=True)
    need={"scientific_name_query","scientificName","latitude","longitude",
          "chelsa_bio01","chelsa_bio15"}
    if not need<=set(source):
        raise ValueError("EAzami source positions or original climate axes missing")
    if source.scientific_name_query.nunique()!=14:
        raise ValueError("Original EAzami taxon concept census changed")
    for name in coverage["threshold_summaries"]["3"]["taxa"]:
        if name not in set(source.scientific_name_query):
            raise ValueError("EAzami n3 panel taxon source absent")
    if len(states)!=20 or not {"accepted_taxon","analysis_state"}<=set(states):
        raise ValueError("Full original floristic U/D source crosswalk lost")
    return source,{
      "original_japan_occurrence_rows":len(a),
      "original_taiwan_occurrence_rows":len(b),
      "original_source_taxa":source.scientific_name_query.nunique(),
      "source_japan_CSV_SHA256":sha(jp),
      "source_taiwan_CSV_SHA256":sha(tw),
      "original_literature_U_D_source_not_S1_photo":True}

def sample_four_axes(source,root,cache_root,metadata):
    lat=source.latitude.to_numpy(float)
    lon=source.longitude.to_numpy(float)
    if not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("Focal original GBIF positions missing")
    out=source.copy()
    source_receipts={}
    for name,member,wanted in (
       ("chelsa_bio18","CHELSA_bio18_1981-2010_V.2.1.tif",B18_TIF_SHA),
       ("chelsa_bio05","CHELSA_bio05_1981-2010_V.2.1.tif",B05_TIF_SHA),
       ("elevation_m","wc2.1_30s_elev.tif",ELEV_TIF_SHA)):
        path=Path(root)/member
        if sha(path.read_bytes())!=wanted:
            raise ValueError("New source raster bytes differ from actual S1/S2")
        values,audit=raster_sample(path,lat,lon)
        out[name]=values
        source_receipts[name]=audit
    with Path(metadata).open() as handle:
        fixed=json.load(handle)["chelsa_rsds_mean"]
    cache=Path(cache_root)/"chelsa_rsds_mean"
    grid=verified_cache(cache,fixed)
    if grid["source_sha256"]!=RSDS_SOURCE_SHA:
        raise ValueError("Original S1/S2 RSDS physical source altered")
    vals=sample_cached_pixels(cache,grid,lat,lon)
    out["chelsa_rsds_mean"]=vals*float(grid["scale"])+float(grid["offset"])
    source_receipts["chelsa_rsds_mean"]={
      "source_raster_SHA256":grid["source_sha256"],
      "original_316_tile_S1_S2_archive_authenticated":True,
      "n_EAzami_source_sites_with_original_RSDS":int(np.isfinite(vals).sum()),
      "n_source_sites_total":len(out),
      "n_sites_without_original_S1_S2_RSDS_tile":int(np.sum(~np.isfinite(vals)))}
    if any(float(np.isfinite(out[name]).mean()) < MIN_VALID_PER_ORIGINAL_TAXON for name in AXES):
        raise ValueError("New four-axis EAzami source site coverage insufficient")
    filtered=out.dropna(subset=list(AXES)).copy()
    rec={
        "original_GBIF_source_observations":len(source),
        "complete_original_GBIF_new_four_axis_observations":len(filtered),
        "n_original_taxa_with_new_four_axes":int(filtered.scientific_name_query.nunique()),
        "original_input_counts_by_species":{str(k):int(v) for k,v in source.groupby("scientific_name_query").size().items()},
        "new_complete_counts_by_species":{str(k):int(v) for k,v in filtered.groupby("scientific_name_query").size().items()},
        "new_original_pixel_source_provenance":source_receipts,
        "no_observation_coordinates_published":True,
        "not_an_independent_new_GBIF_occurrence_sample":True}
    return filtered,rec

def original_method_reproduction(core,crosswalk,coverage,contract,legacy,rawtrees,occ):
    taxa=list(coverage["threshold_summaries"]["10"]["taxa"])
    env=occ.groupby("scientific_name_query")[["chelsa_bio15","chelsa_bio01"]].mean().loc[taxa].copy()
    env["z15"]=core.zscore(env["chelsa_bio15"].to_numpy(float))
    env["z1"]=core.zscore(env["chelsa_bio01"].to_numpy(float))
    state=core.panel_state_map(crosswalk,taxa)
    assets=core.prepare_topology_assets(rawtrees,taxa,env)
    scores=core.observed_topology_stats(assets,state)
    reproduction=core.method_check(contract,legacy,scores)
    if reproduction["status"]!="pass":
        raise ValueError("Original EAzami original BIO1/BIO15 CTMC validation failed")
    return {
        "original_EAzami_BIO15_BIO1_method_check":reproduction,
        "original_frozen_9taxon_n_U":5,
        "original_frozen_9taxon_n_D":4,
        "frozen_original_branchwise_method_exactly_reused":True}

def prepare_environment(core,taxa,occ):
    eligible=occ.loc[occ.scientific_name_query.isin(taxa)]
    d=eligible.groupby("scientific_name_query")[list(AXES)].mean().loc[taxa]
    z={name:core.zscore(d[name].to_numpy(float)) for name in AXES}
    return {core.normalize_tip(taxa[i]):{x:float(z[x][i]) for x in AXES}
          for i in range(len(taxa))}

def create_assets(core,rawtrees,taxa,zvalues):
    assets=[]
    for original in rawtrees:
        tree=core.prune_to_panel(original,taxa)
        edge={}
        for label in AXES:
            mapping={tip:values[label] for tip,values in zvalues.items()}
            _,shift=core.brownian_internal_values(tree,mapping)
            edge[label]=shift
        assets.append({"tree":tree,"edge":edge})
    if len(assets)!=6:
        raise ValueError("EAzami six original accepted nuclear topologies not preserved")
    return assets

def exact_transition_vector(core,tree,state,edges):
    q=core.fit_symmetric_q(tree,state)
    probability,joint=core.ctmc_likelihood_and_edge_joint(tree,state,q,need_joint=True)
    if joint is None or not probability>0:
        raise ValueError("Original EAzami CTMC edge likelihood invalid")
    numerator=np.zeros(len(AXES))
    denom=0.0
    signed=[]
    for ci,J in joint.items():
        delta=float(J[0,1]-J[1,0])
        changes=float(J[0,1]+J[1,0])
        denom+=changes
        signed.append(delta)
        for j,name in enumerate(AXES):
            numerator[j]+=delta*edges[name][ci]
    if denom<=EPS or not np.isfinite(numerator).all():
        raise ValueError("EAzami numerator/transition denominator cannot be estimated")
    return numerator/denom,denom,float(q)

def rank(x,observed,greater=True):
    scores=np.asarray(x,float)
    if greater:
        extreme=int(np.sum(scores>=observed-EPS))
    else:
        extreme=int(np.sum(scores<=observed+EPS))
    return {"extreme_state_maps":extreme,
            "all_fixed_count_state_maps":len(scores),
            "conditional_exact_state_assignment_fraction":float(extreme/len(scores)),
            "not_a_biological_replicate_or_frequentist_P_value":True}

def exact_panel(core,asset,counts,states,*,max_maps):
    taxa=list(counts)
    tipnames=[core.normalize_tip(t) for t in taxa]
    d=sum(int(states[t]) for t in tipnames)
    combinations=itertools.combinations(range(len(taxa)),d)
    total=math.comb(len(taxa),d)
    if total!=max_maps:
        raise ValueError("Original EAzami count-preserving null universe changed")
    all_stats=np.empty((total,len(AXES)))
    original_map=np.empty(len(AXES))
    observed_once=0
    original_topo=None
    q_observed=[]
    transitions_observed=[]
    for n,comb in enumerate(combinations):
        choose=set(comb)
        testmap={tipnames[i]:(1 if i in choose else 0) for i in range(len(taxa))}
        values=[];expected=[];qvalues=[]
        for row in asset:
            z,events,q=exact_transition_vector(core,row["tree"],testmap,row["edge"])
            values.append(z)
            expected.append(events)
            qvalues.append(q)
        arr=np.asarray(values)
        median=np.median(arr,axis=0)
        all_stats[n]=median
        if all(testmap[t]==states[t] for t in tipnames):
            original_map[:]=median
            original_topo=arr.tolist()
            q_observed=qvalues
            transitions_observed=expected
            observed_once+=1
    if observed_once!=1 or original_topo is None:
        raise ValueError("EAzami original published state map absent/duplicated")
    j=AXES.index("chelsa_bio18")
    projected=np.sqrt(np.sum(all_stats**2,axis=1))
    obs4=float(np.linalg.norm(original_map))
    metrics={}
    for i,field in enumerate(AXES):
        obs=float(original_map[i])
        null=all_stats[:,i]
        metrics[field]={
            "original_EAzami_six_topology_median_signed_U_to_D_minus_D_to_U_environment_shift":obs,
            "by_topology": [float(x[i]) for x in original_topo],
            "direction_consistent_on_six_topologies":
                len(set(np.sign(np.array(original_topo)[:,i])))==1 and obs!=0,
            "positive_fraction":rank(null,obs,True),
            "negative_fraction":rank(null,obs,False),
            "two_sided_absolute_fraction":rank(np.abs(null),abs(obs),True),
            "null_map_quantiles_q025_q50_q975":[float(t) for t in np.quantile(null,[.025,.5,.975])],
            "one_phenotype_tip_map_drawn_before_environment_sampling":True,
        }
    return {
       "original_EAzami_frozen_tip_count":len(taxa),
       "original_EAzami_frozen_n_U":len(taxa)-d,
       "original_EAzami_frozen_n_D":d,
       "state_maps_exhaustively_enumerated":total,
       "actual_original_published_state_assignment_included_once":True,
       "same_EAzami_6_accepted_nuclear_topologies":True,
       "fixed_count_null_state_map_is_not_photographic_phenotype":True,
       "original_state_map_fit_rates_six_topologies":list(map(float,q_observed)),
       "original_state_map_expected_changes_six_topologies":list(map(float,transitions_observed)),
       "new_four_axes":metrics,
       "S1_anchored_positive_BIO18_hypothesis":{
           "scope":"posthoc_axis_chosen_from_original_S1_positive_BIO18_species_association",
           "predeclared_direction":"higher_BIO18_on_U_to_D_edges",
           "original_signed_transition_BIO18":float(original_map[j]),
           "expected_direction_match_all_six_topologies":all(z[j]>0 for z in original_topo),
           "exact_positive_state_map_rank":metrics["chelsa_bio18"]["positive_fraction"],
           "not_frozen_original_EAzami_BIO15_BIO1_hypothesis":True,
        },
       "four_axis_direction_free_omnibus":{
           "statistic":"Euclidean_norm_of_four_six_topology_median_signed_standardized_environment_edge_shifts",
           "observed_norm":obs4,
           "conditional_frozen_count_state_map_rank":rank(projected,obs4,True),
           "preserves_between_axis_covariance_through_joint_state_assignment_null":True,
           "no_posthoc_direction_search":True,
        },
    }

def execute(args):
    core=original_implementation(args.original_eazami_algorithm)
    raw_trees=checked(args.eazami_six_trees,EA_TREES_BLOB)
    checked(args.eazami_crosswalk,EA_STATES_BLOB)
    checked(args.eazami_coverage,EA_COVERAGE_BLOB)
    checked(args.eazami_contract,EA_CONTRACT_BLOB)
    from io import StringIO
    trees=core.read_trees(args.eazami_six_trees,6)
    cross=pd.read_csv(args.eazami_crosswalk)
    coverage=json.loads(args.eazami_coverage.read_text())
    contract=json.loads(args.eazami_contract.read_text())
    legacy=json.loads(args.eazami_legacy.read_text())
    all_occ,original_proof=original_occurrence_cohort(args.japan_occurrences,args.taiwan_occurrences,cross,coverage)
    method=original_method_reproduction(core,cross,coverage,contract,legacy,trees,all_occ)
    measured,source=sample_four_axes(all_occ,args.raster_root,args.rsds_cache_root,args.rsds_metadata)
    original_counts=all_occ.groupby("scientific_name_query").size().to_dict()
    current=measured.groupby("scientific_name_query").size().to_dict()
    results={}
    for label,(key,n,nu,nd,nmaps) in PANELS.items():
        taxa=list(coverage["threshold_summaries"][key]["taxa"])
        if len(taxa)!=n:raise ValueError("Original EAzami panel membership not fixed")
        if any(original_counts.get(t,0)<int(key) for t in taxa):
            raise ValueError("Original EAzami specimen occurrence threshold not met")
        insufficient={t:int(current.get(t,0)) for t in taxa if current.get(t,0)<int(key)}
        if insufficient:
            results[label]={"status":"STOP_NEW_FOUR_CLIMATE_OCCURRENCE_SOURCE_TOO_FEW",
                "n_original_frozen_panel_taxa":n,
                "n_failed_species":len(insufficient),
                "no_values_imputed":True}
            continue
        mapst=core.panel_state_map(cross,taxa)
        if (sum(mapst.values()),len(mapst))!=(nd,n):
            raise ValueError("EAzami frozen U D concept changed")
        env=prepare_environment(core,taxa,measured)
        assets=create_assets(core,trees,taxa,env)
        model=exact_panel(core,assets,taxa,mapst,max_maps=nmaps)
        model["status"]="EAzami_AUTHENTIC_BRANCHWISE_UD_CTMC_FOUR_CLIMATE_EXACT_MAP_COMPLETED"
        model["n_original_source_occurrence_rows_before_4_climate_extension"]=sum(int(original_counts[t]) for t in taxa)
        model["n_new_4environment_complete_original_occurrence_rows"]=sum(int(current[t]) for t in taxa)
        model["post_result_source_driven_axis_extension_not_independent_confirmation"]=True
        results[label]=model
        print("EAzami_REAL_FOUR_ENV_UD_TRANSITION="+json.dumps({
           "panel":label,
           "n":n,
           "maps":nmaps,
           "BIO18_signed":model["S1_anchored_positive_BIO18_hypothesis"]["original_signed_transition_BIO18"],
           "BIO18_topology_consistent":model["S1_anchored_positive_BIO18_hypothesis"]["expected_direction_match_all_six_topologies"],
           "BIO18_map_rank":model["S1_anchored_positive_BIO18_hypothesis"]["exact_positive_state_map_rank"],
           "four_axis_map_rank":model["four_axis_direction_free_omnibus"]["conditional_frozen_count_state_map_rank"],
        },sort_keys=True),flush=True)
    report={
       "schema":"azami_v4_20261024_EAzami_ORIGINAL_exact_ud_CTMC_BIO18_BIO05_RSDS_elev_S1S2_T2_v1",
       "status":"ORIGINAL_EAzami_TRANSITION_WORKFLOW_REEXECUTED_WITH_S1S2_FOUR_ENVIRONMENT_SITE_VALUES",
       "scientific_target":"Does source-authenticated U->D branch probability preferentially accompany larger warmest-quarter rain BIO18, and/or a non-directional 4-axis niche shift, in the ORIGINAL frozen EAzami East-Asian panel?",
       "scientific_role":"posthoc_S1_BIO18_motivated_transition_workflow_extension_not_confirmation_of_original_EAzami_BIO1_BIO15_and_not_a_S1_photo_angle_measurement",
       "same_4_axis_environment_definitions_as_S1_S2":["chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m"],
       "source_climatic_epochs":{"CHELSA_v2_1":"1981-2010","elevation":"WorldClim_v2_1_static_30arcsec"},
       "exact_frozen_EAzami_algorithm_git_blob":EA_CORE_GIT_BLOB,
       "exact_EAzami_accepted_nuclear_topologies_git_blob":EA_TREES_BLOB,
       "fixed_20tip_botanical_U_D_source_git_blob":EA_STATES_BLOB,
       "original_sample_proof":original_proof,
       "original_BIO1_BIO15_algorithm_numeric_reproduction":method,
       "same_original_GBIF_sites_new_four_axis_coverage":source,
       "panels":results,
       "primary_panel":"n5_primary",
       "side_panels_are_sensitivity_not_independent_replicates":True,
       "no_trait_photo_measurement_or_true_gravity_conflation":True,
       "no_independent_molecular_UD_origins_counted":True,
       "no_transition_time_real_climate_or_organismal_fitness_evidence":True,
       "null_is_exact_finite_state_map_rank_not_independent_replicate_P_value":True,
       "two_trees_and_S1_S2_overlap_not_independent_evolutionary_replications":True,
       "conditional_on_original_taxon_states_six_trees_Brownian_internal_reconstruction":True,
       "no_model_selection_from_BIO18_BIO05_RSDS_elevation_outcomes":True,
       "no_raw_observation_location_user_or_image_data_output":True,
    }
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    Path(args.out).write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    return report

def main():
    p=argparse.ArgumentParser()
    for name in ("original-eazami-algorithm","eazami-six-trees","eazami-crosswalk",
                 "eazami-coverage","eazami-contract","eazami-legacy",
                 "japan-occurrences","taiwan-occurrences","raster-root",
                 "rsds-cache-root","rsds-metadata","out"):
        p.add_argument("--"+name.replace(" ",""),type=Path,required=True)
    a=p.parse_args()
    execute(a)
if __name__=="__main__":
    main()
