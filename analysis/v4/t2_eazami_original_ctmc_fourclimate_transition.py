"""Reapply the EXACT EAzami 2026-09-04 U->D CTMC/Brownian workflow.

This is NOT S1's EXIF-relative continuous image angle or a PIC/PGLS substitute.
The original nine (5 U/4 D), twelve (7 U/5 D) and thirteen (7 U/6 D)
genuine Japan/Taiwan source taxa, six original accepted AU Comp1061 multilocus
NUCLEAR trees, source-specific real occurrences and exact tip-label null are
preserved. EAzami's ORIGINAL BIO15/BIO1 effect is reproduced before touching
the new exposures. Four NEW descriptors BIO18 (warmest-quarter rainfall),
BIO05 (max temperature warmest month), original CHELSA RSDS and WorldClim
DEM are sampled at precisely the old source occurrences.

Primary user-directed directional test: EAzami U->D aligns with increasing
BIO18, suggested by S1 *between-taxon* image-angle result. Posthoc omnibus
four-variable and sign-free individual axes are exploratory. Six nuclear-tree
statistics are medians across fixed topologies (NOT six replicates).
Original status, read-only; do not edit or reclassify original EAzami states.
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
import zipfile

import numpy as np
import pandas as pd
from scipy import linalg

from analysis.v4.s1_s2_real_bio18_bio05_rsds_elevation import raster_sample
from analysis.v4.native_cached_climate import verified_cache, sample_cached_pixels

AXES=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
TREE_BLOB="d7f8c0b96ee9c257c815db06a5a0aa9fdae4841e"
SCRIPT_BLOB="09731b563907b3ec2778ea04086d729771376f7e"
CROSSWALK_BLOB="f655544e77f1659323c1e915567dacbb51b01cc3"
N_ORIG_OCC={"Japan":693,"Taiwan":55}
STATE_PANELS={"10":(9,5,4,126),"5":(12,7,5,792),"3":(13,7,6,1716)}
LAYER_RASTERS={
    "chelsa_bio18":"CHELSA_bio18_1981-2010_V.2.1.tif",
    "chelsa_bio05":"CHELSA_bio05_1981-2010_V.2.1.tif",
    "elevation_m":"wc2.1_30s_elev.tif"}
EXPECTED_RASTER_SHA={
 "chelsa_bio18":"68b0b1cb905cdf6972364e947e137f714a68d287f38863b1516161a2a4bbd7ec",
 "chelsa_bio05":"90078d99eeb0b08dd1a115a5ffa60b7c19e5849686854db950283cbac9d3817d",
 "elevation_m":"1d57ce9fd3584cea6a227396bc75ef88522fd5346e8487c17370885c0499716e"}
EXPECTED_RSDS_SHA="cc1ad14a70a1a940be3d01d918d4c3f7e35db95b86b12af9be5487e1bbc1fc7d"
LEGACY_BIO15=[.267907,.269417]
LEGACY_BIO1=[-.199306,-.199151]

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def reference_module(root):
    p=Path(root)/"analysis"/"run_chapter2_orientation_transition_regime_hypothesis_v1.py"
    raw=p.read_bytes()
    if git_blob(raw)!=SCRIPT_BLOB:
        raise ValueError("EAzami reference-original transition source Git blob changed")
    spec=importlib.util.spec_from_file_location("eazami_original_frozen_ctmc",p)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def original_assets(root):
    base=Path(root)/"data"/"evidence"
    tr=base/"comp1061_six_au_nonrejected_trees_v1.nwk"
    cross=base/"orientation_comp1061_20tip_source_crosswalk_v1.csv"
    if git_blob(tr.read_bytes())!=TREE_BLOB:
        raise ValueError("Original 6 accepted multilocus nuclear trees changed")
    if git_blob(cross.read_bytes())!=CROSSWALK_BLOB:
        raise ValueError("Original EAzami U/D authority crosswalk changed")
    cov=json.loads((base/"chapter2_orientation_occurrence_coverage_audit_result_v1.json").read_text())
    contract=json.loads((base/"chapter2_orientation_transition_regime_hypothesis_contract_v1.json").read_text())
    legacy=json.loads((base/"fdt4_branchwise_niche_transition_concordance_v1.json").read_text())
    if contract["version"]!="chapter2_orientation_transition_regime_hypothesis_contract_v1":
        raise ValueError("Unexpected original biological hypothesis")
    return tr,pd.read_csv(cross),cov,contract,legacy

def occurrences_from_zip(path,expected):
    raw=Path(path).read_bytes()
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            member="focal_occurrences_environment_complete.csv"
            if member not in archive.namelist():
                raise ValueError("Missing original EAzami exact occurrence source")
            frame=pd.read_csv(archive.open(member),low_memory=False)
    else:
        # GitHub Actions 'gh run download --name' extracts original ZIP into
        # its authenticated CSV member. These are identical source rows.
        frame=pd.read_csv(path,low_memory=False)
    if len(frame)!=N_ORIG_OCC[expected] or not {"scientific_name_query",
       "latitude","longitude","chelsa_bio01","chelsa_bio15",
       "chelsa_bio12","chelsa_bio04"}<=set(frame):
        raise ValueError("Changed original EAzami Japan/Taiwan niche source cohort")
    return frame,sha(raw)

def original_reproduction(ref,orig,trees,crosswalk,coverage,contract,legacy):
    taxa=coverage["threshold_summaries"]["10"]["taxa"]
    n,u,d,maps=STATE_PANELS["10"]
    if len(taxa)!=n or len(taxa)!=len(set(taxa)):
        raise ValueError("Frozen EAzami n10 panel changed")
    states=ref.panel_state_map(crosswalk,taxa)
    if (sum(1 for v in states.values() if v==0),sum(states.values()))!=(u,d):
        raise ValueError("Original trait state baseline changed")
    counts,env=ref.build_panel_environment(orig,taxa)
    if not all(counts.loc[t]>=10 for t in taxa):
        raise ValueError("Original source does not reproduce n10 occurrence support")
    assets=ref.prepare_topology_assets(trees,taxa,env)
    stats=ref.observed_topology_stats(assets,states)
    receipt=ref.method_check(contract,legacy,stats)
    # Source-validated original U->D low BIO1 and high BIO15 before testing
    # ANY user-requested new climate or new photos.
    if receipt["status"]!="pass":
        raise ValueError("Frozen original EAzami CTMC numerical method failed")
    return receipt

def append_new_original_site_environment(occ,source_root,rsds_root,rsds_metadata):
    raw=occ.copy()
    lat=raw.latitude.to_numpy(float)
    lon=raw.longitude.to_numpy(float)
    audit={}
    for field,raster in LAYER_RASTERS.items():
        x,metadata=raster_sample(Path(source_root)/raster,lat,lon)
        if metadata["source_raster_sha256"]!=EXPECTED_RASTER_SHA[field]:
            raise ValueError("New physical CHELSA/WorldClim environmental raster SHA differs from S1/S2")
        raw[field]=x
        audit[field]=metadata
    contract=json.loads(Path(rsds_metadata).read_text())["chelsa_rsds_mean"]
    grid=verified_cache(Path(rsds_root)/"chelsa_rsds_mean",contract)
    if grid["source_sha256"]!=EXPECTED_RSDS_SHA:
        raise ValueError("Frozen original S1/S2 RSDS grid changed")
    x=sample_cached_pixels(Path(rsds_root)/"chelsa_rsds_mean",grid,lat,lon)
    raw["chelsa_rsds_mean"]=x*grid["scale"]+grid["offset"]
    audit["chelsa_rsds_mean"]={
       "same_original_316_tile_CHELSA_S1_S2_source_sha256":grid["source_sha256"],
       "n_complete":int(np.isfinite(raw["chelsa_rsds_mean"]).sum()),
       "n_total_original_occurrences":len(raw)}
    for v in AXES:
        raw[v]=pd.to_numeric(raw[v],errors="coerce")
    return raw,audit

def panel_data(ref,frame,taxa,threshold):
    required=set(AXES)
    valid=frame.loc[np.isfinite(frame[list(AXES)]).all(axis=1)].copy()
    counts=valid.groupby("scientific_name_query").size()
    missing={t:int(counts.get(t,0)) for t in taxa if int(counts.get(t,0))<threshold}
    if missing:
        return None,{"status":"STOP_INCOMPLETE_EAZAMI_ORIGINAL_TAXON_4_ENVIRONMENT_SITE_COVERAGE",
                     "minimum_source_observations":threshold,"n_taxa_original":len(taxa),
                     "missing":missing}
    env=valid.groupby("scientific_name_query")[list(AXES)].mean().loc[taxa].copy()
    for field in AXES:
        if not np.isfinite(env[field]).all() or env[field].std(ddof=1)<=0:
            raise ValueError("Real new CHELSA four-axis source environment variance absent")
        env["z_"+field]=ref.zscore(env[field].to_numpy(float))
    return env,{"status":"ORIGINAL_EAZAMI_TAXON_COHORT_FOUR_SOURCE_FIELDS_COMPLETE",
                "n_original_source_taxa":len(taxa),
                "n_original_source_occurrences":int(sum(int(counts.loc[t]) for t in taxa)),
                "by_taxon_four_climate_complete_counts":{t:int(counts.loc[t]) for t in taxa},
                "taxon_centroid_means_physical_units":{t:{a:float(env.loc[t,a]) for a in AXES}
                    for t in taxa}}

def tree_environment_asset(ref,tree,taxa,env):
    tr=ref.prune_to_panel(tree,taxa)
    vals={a:{ref.normalize_tip(t):float(env.loc[t,"z_"+a]) for t in taxa} for a in AXES}
    deltas={}
    for a in AXES:
        _,edge=ref.brownian_internal_values(tr,vals[a])
        deltas[a]=edge
    edges=sorted(deltas[AXES[0]])
    if any(set(deltas[a])!=set(edges) for a in AXES):
        raise ValueError("Brownian environment branch indexing drifted")
    matrix=np.column_stack([[deltas[a][edge] for edge in edges] for a in AXES])
    return tr,edges,matrix

def one_state_map(ref,assets,states):
    scores=[]
    expected=[]
    for tree,edges,climate in assets:
        q=ref.fit_symmetric_q(tree,states)
        likelihood,J=ref.ctmc_likelihood_and_edge_joint(tree,states,q,need_joint=True)
        if J is None or not np.isfinite(likelihood) or likelihood<=0:
            raise ValueError("True CTMC edge joint not identifiable")
        change=float(sum(J[c][0,1]+J[c][1,0] for c in edges))
        if change<=1e-12:raise ValueError("No posterior transition mass")
        signed=np.array([float(J[c][0,1]-J[c][1,0]) for c in edges])
        scores.append((signed@climate)/change)
        expected.append(change)
    arr=np.asarray(scores,float)
    return np.median(arr,axis=0),arr,np.array(expected)

def all_state_maps(ref,assets,taxa,observed_states,expect_maps):
    names=[ref.normalize_tip(t) for t in taxa]
    nD=sum(observed_states.values())
    combinations=list(itertools.combinations(range(len(names)),nD))
    if len(combinations)!=expect_maps:
        raise ValueError("Exact EAzami n5/n3 conditional state-map universe changed")
    medians=np.empty((len(combinations),len(AXES)),dtype=float)
    observed_index=None
    topology_obs=None
    changes_obs=None
    for i,combo in enumerate(combinations):
        dset=set(combo)
        state={name:int(j in dset) for j,name in enumerate(names)}
        med,topo,changes=one_state_map(ref,assets,state)
        medians[i]=med
        if state==observed_states:
            if observed_index is not None:raise ValueError("Observed EAzami states duplicated")
            observed_index=i
            topology_obs=topo
            changes_obs=changes
    if observed_index is None:
        raise ValueError("Frozen EAzami state map absent from conditional permutations")
    assert np.allclose(medians[observed_index],np.median(topology_obs,axis=0))
    observed=medians[observed_index]
    m=len(combinations)
    axis_result={}
    for j,a in enumerate(AXES):
        row=medians[:,j]
        median=float(np.median(row))
        two_sided=float(np.sum(np.abs(row-median)>=abs(observed[j]-median)-1e-12)/m)
        plus=float(np.sum(row>=observed[j]-1e-12)/m)
        minus=float(np.sum(row<=observed[j]+1e-12)/m)
        axis_result[a]={
          "six_topology_observed_scores":[float(t) for t in topology_obs[:,j]],
          "six_topology_median_observed_score":float(observed[j]),
          "null_median":median,
          "null_q025_q975":[float(q) for q in np.quantile(row,[.025,.975])],
          "exact_fraction_at_least_observed_in_U_to_D_increase_direction":plus,
          "exact_fraction_at_least_observed_in_U_to_D_decrease_direction":minus,
          "exact_two_sided_fraction_around_null_median":two_sided,
          "all_source_climate_branch_changes_Brownian_not_observed_historical_climate":True}
    # The four variables' null covariance accounts for strong collinearity.
    # Center the null and use a fixed pseudoinverse, deterministic exhaustive
    # rank; no post hoc direction search across environmental variables.
    mean=medians.mean(0)
    cov=np.cov(medians,rowvar=False)
    eig,rotation=np.linalg.eigh(cov)
    positive=eig>max(float(eig.max())*1e-8,1e-12)
    if int(positive.sum())<2:
        raise ValueError("Observed four-axis state-null covariance nonidentifiable")
    matrix=(rotation[:,positive]/eig[positive])@rotation[:,positive].T
    gaps=medians-mean
    sq=np.einsum("ni,ij,nj->n",gaps,matrix,gaps)
    stat=float(sq[observed_index])
    rank=float(np.mean(sq>=stat-1e-12))
    return {
      "n_original_taxa":len(taxa),
      "n_orig_U":len(taxa)-nD,"n_orig_D":nD,
      "original_taxa":list(taxa),
      "number_count_preserving_tip_state_maps":m,
      "fixed_original_genuine_accepted_nuclear_topologies":len(assets),
      "reestimated_ctmc_symmetric_q_per_topology_and_null_map":True,
      "signed_U_to_D_minus_D_to_U_edge_joint_weight":True,
      "axis_results":axis_result,
      "BIO18_focal_one_sided_increase_exact_fraction":
          axis_result["chelsa_bio18"][
            "exact_fraction_at_least_observed_in_U_to_D_increase_direction"],
      "all_four_posthoc_omnibus_state_map_rank":{
        "statistic":"squared_correlated_null_whitened_norm_of_four_environment_changes",
        "null_effective_rank":int(positive.sum()),
        "squared_distance_observed":stat,
        "exact_fraction_all_maps_as_or_more_exceptional":rank,
        "not_pre_registered_in_original_EAzami":True},
      "observed_branchwise_changes_each_of_six_trees":[float(v) for v in changes_obs],
      "phylogenetic_transition_time_climate_or_fitness_not_observed":True,
      "finite_label_map_rank_not_biological_replicate_probability":True,
    }

def execute(reference,japan,taiwan,source_root,rsds_root,rsds_meta,out):
    ref=reference_module(reference)
    treefile,crosswalk,coverage,contract,legacy=original_assets(reference)
    jp,jsha=occurrences_from_zip(japan,"Japan")
    tw,tsha=occurrences_from_zip(taiwan,"Taiwan")
    frame=pd.concat([jp,tw],ignore_index=True)
    six=ref.read_trees(treefile,6)
    if len(six)!=6:raise ValueError("EAzami original six non-AU-rejected topologies not recovered")
    reproduced=original_reproduction(ref,frame,six,crosswalk,coverage,contract,legacy)
    enriched,source_meta=append_new_original_site_environment(frame,source_root,rsds_root,rsds_meta)
    report={}
    for key in ("5","3","10"):
        n,U,D,maps=STATE_PANELS[key]
        taxa=coverage["threshold_summaries"][key]["taxa"]
        if len(taxa)!=n:raise ValueError("Original historical EAzami panel no longer same n")
        state=ref.panel_state_map(crosswalk,taxa)
        if (sum(s==0 for s in state.values()),sum(state.values()))!=(U,D):
            raise ValueError("Original EAzami CTMC states changed")
        env,diagnosis=panel_data(ref,enriched,taxa,int(key))
        name={"5":"n5_original_primary_12","3":"n3_original_sensitivity_13",
              "10":"n10_original_method_9"}[key]
        if env is None:
            report[name]={"status":"ORIGINAL_REAL_FOUR_ENVIRONMENT_INCOMPLETE",
                          "source_coverage_diagnosis":diagnosis}
            continue
        assets=[tree_environment_asset(ref,tr,taxa,env) for tr in six]
        scores=all_state_maps(ref,assets,taxa,state,maps)
        scores["source_coverage_diagnosis"]=diagnosis
        scores["status"]="ORIGINAL_EAzami_SIX_NUCLEAR_CTMC_BROWNIAN_STATE_MAPS_4CLIMATE_COMPLETE"
        report[name]=scores
        print("EAzami_TRUE_U_TO_D_FOUR_AXIS_PANEL="+json.dumps({
           "panel":name,"n":n,
           "BIO18_stat":scores["axis_results"]["chelsa_bio18"]["six_topology_median_observed_score"],
           "BIO18_expected_increase_exact_fraction":scores["BIO18_focal_one_sided_increase_exact_fraction"],
           "four_axes_omnibus_exact_fraction":
               scores["all_four_posthoc_omnibus_state_map_rank"]["exact_fraction_all_maps_as_or_more_exceptional"]
          },sort_keys=True),flush=True)
    doc={
      "schema":"EAzami_original_JP_TW_nuclear_UD_CTMC_Brownian_REAL_BIO18_BIO05_RSDS_elevation_20261024_v1",
      "status":"ORIGINAL_EAZAMI_HISTORICAL_BINARY_TRANSITION_WORKFLOW_REPLAYED_WITH_4_REAL_SITE_CLIMATES",
      "science_role":"Complement S1/S2 real four-climate continuous image-angle result with independent EAzami literature U/D lineage-transition evidence",
      "source_names":"original EAzami 20-tip authority crosswalk and unchanged occurrence-coverage n5/n3/n10 panels",
      "original_old_BIO1_BIO15_method_reproduced_before_new_axes":reproduced,
      "same_original_six_AU_nonrejected_nuclear_topologies":True,
      "original_EAzami_trait":"botanical-literature U(0) / D(1), NOT S1 photo angle",
      "four_real_site_climates":list(AXES),
      "new_source_provenance":source_meta,
      "n_original_Japan_occurrences":len(jp),"n_original_Taiwan_occurrences":len(tw),
      "source_original_archive_sha256":{"japan":jsha,"taiwan":tsha},
      "original_EAzami_reference_git_blob":{
         "hypothesis_method":SCRIPT_BLOB,"six_AU_NWK":TREE_BLOB,
         "botanical_U_D_authority_crosswalk":CROSSWALK_BLOB},
      "original_BIO1_BIO15_frozen_result_unmodified":True,
      "new_axes_user_requested_posthoc":True,
      "focal_BIO18_direction_selected_from_S1_among_Cirsium_taxa_not_from_EAzami_transition_outputs":True,
      "per_taxon_coordinates_and_observer_data_never_exported":True,
      "all_axes_exploratory_exact_ranks_not_p_values_of_independent_biological_events":True,
      "no_transition_time_temperature_precipitation_radiation_or_elevation_observed":True,
      "EAzami_binary_botanical_tip_states_not_interchangeable_with_S1_0_to_180_photo_angles":True,
      "not_a_T3_external_tribe_S1_photo_angle_replication":True,
      "panels":report}
    root=Path(out);root.mkdir(parents=True,exist_ok=True)
    (root/"EAzami_UD_original_method_4_REAL_CLIMATES_CTMC_BRANCH_RESULTS.json").write_text(
       json.dumps(doc,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    return doc

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for s in ("reference","japan-zip","taiwan-zip","source-root","rsds-root","rsds-metadata","outdir"):
        p.add_argument("--"+s,required=True,type=Path)
    a=p.parse_args()
    execute(a.reference,a.japan_zip,a.taiwan_zip,a.source_root,a.rsds_root,
            a.rsds_metadata,a.outdir)
if __name__=="__main__":
    main()
