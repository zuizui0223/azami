"""Independent North American EAzami exact U->D CTMC transition × S1/S2 climates.

Apply the ORIGINAL EAzami pin-verified symmetric CTMC edge posterior,
Brownian squared-change ancestral environmental values, and exhaustive
count-preserving tip U/D map null to Siniscalchi 2023 author's genuinely
multilocus 63-tip nuclear species tree (strict original S1 42-tip subset).
Author's independently described posture is a distinct botanical U/D trait,
NEVER obtained by binarizing S1 camera-photo angle. No cpDNA, grafted
placement, single-locus ITS or failed candidate nuclear gene tree is used.
BIO18/BIO5/RSDS/elevation come from actual S1/S2 original site pixels
fixed BEFORE the NA author posture data are used. No historical rainfall
selection, gravity angle or independent gain count is inferred.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,itertools,json,math
from pathlib import Path
from Bio import Phylo
import numpy as np
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as dna
from analysis.v4 import eazami_ud_real_four_environment_transition_regime as ea
from analysis.v4 import t2_s1_fourclimate_real_nuclear_pgls_pic as four
from analysis.v4.t2_pgls_placement_screen import read_author_states

NA_TRAITS="traits/traits.csv"
NA_TRAIT_BLOB="fa127f0587b676e6e004ea121424d26069586bcb"
FOUR=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
SOURCE_S1_SHA=dna.S1_SHA
SOURCE_NEW4_SHA=four.NEW_AGGREGATE_SHA

def read251(path,sha):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:raise ValueError("Unverified original S1 source or four-site raster means")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    by={x["original_source_name"]:x for x in rows}
    if len(rows)!=251 or len(by)!=251:raise ValueError("251 original S1 names drifted")
    return by

def inputs(original,new4,tree_raw,trait_raw):
    orig=read251(original,SOURCE_S1_SHA)
    four=read251(new4,SOURCE_NEW4_SHA)
    if set(orig)!=set(four):raise ValueError("Same-original-photo taxon matching changed")
    if dna.git_blob(tree_raw)!=dna.NA_GIT_BLOB:raise ValueError("True NA nuclear species tree Git blob drift")
    if dna.git_blob(trait_raw)!=NA_TRAIT_BLOB:raise ValueError("NA author botanical posture raw source drift")
    states=read_author_states(trait_raw) # independent full 63 U/D authority.
    tr=Phylo.read(io.StringIO(tree_raw.decode()),"newick")
    if len(tr.get_terminals())!=63 or {x.name for x in tr.get_terminals()}!=set(states):
        raise ValueError("NA published 63 original species/tree mismatch")
    names=sorted(n for n,x in orig.items()
                 if n in states and x["WCVP_key_collides_source_names"]=="False")
    if len(names)!=42 or sum(states[n] for n in names)!=3:
        raise ValueError("Genuine NA42 original D3 vs U39 roster changed")
    if sum(int(orig[n]["original_photo_count"]) for n in names)!=59917:
        raise ValueError("NA42 original 59917 camera angle photo records changed")
    for n in names:
        a,b=orig[n],four[n]
        if a["accepted_key"]!=b["accepted_key"]:
            raise ValueError("Photographic species taxon keys mismatched")
        if int(b["n_photo_complete_four_axis"])>int(a["original_photo_count"]):
            raise ValueError("New bio18 photos cannot exceed original")
        if not 0<=float(b["mean_image_angle_deg"])<=180:
            raise ValueError("S1 measured camera angles invalid")
        if not all(np.isfinite(float(b[z])) for z in FOUR):
            raise ValueError("Incomplete source observations for novel four climate")
    n10=[n for n in names if int(four[n]["n_photo_complete_four_axis"])>=10]
    if len(n10)!=40:raise ValueError("Original >=10 image count nuclear cohort NA40 drift")
    return orig,four,states,tr,names,n10

def environmental_edges(h1,tr,taxa,source):
    m=np.array([[float(source[n][f]) for f in FOUR] for n in taxa])
    if not np.isfinite(m).all() or min(np.std(m,axis=0))<=0:
        raise ValueError("Nonfinite or constant frozen source climates")
    return {axis:h1.brownian_internal_values(tr,dict(zip(taxa,h1.zscore(m[:,i]))))[1]
            for i,axis in enumerate(FOUR)}

def transition_score(h1,tr,state,env):
    q=h1.fit_symmetric_q(tr,state)
    L,joint=h1.ctmc_likelihood_and_edge_joint(tr,state,q,need_joint=True)
    if joint is None or L<=0:raise ValueError("Cannot fit original EAzami CTMC")
    den=sum(float(J[0,1]+J[1,0]) for J in joint.values())
    if den<=0:raise ValueError("Invalid posterior CTMC change mass")
    signed={c:float(J[0,1]-J[1,0]) for c,J in joint.items()}
    out=np.array([sum(signed[c]*env[z][c] for c in signed)/den for z in FOUR])
    if not np.isfinite(out).all():raise ValueError("Invalid original EAzami environmental shift")
    return out,float(q),float(den)

def exact_fit(h1,tree,by,states,names,tag):
    names=sorted(names)
    t=h1.prune_to_panel(tree,names)
    obs={n:int(states[n]) for n in names}
    d=sum(obs.values())
    count=math.comb(len(names),d)
    if (len(names),d,count)==(42,3,11480):pass
    elif len(names)==40 and d>=2 and count<=20000:pass
    else:raise ValueError("Unexpected true author nuclear finite state map count")
    env=environmental_edges(h1,t,names,by)
    original,q,expected=transition_score(h1,t,obs,env)
    score=np.empty((count,4),dtype=float)
    matches=0
    for i,group in enumerate(itertools.combinations(range(len(names)),d)):
        marked=set(group)
        state={name:int(j in marked) for j,name in enumerate(names)}
        score[i]=transition_score(h1,t,state,env)[0] # RE-FIT q per exact map
        if state==obs:
            matches+=1
            if not np.allclose(original,score[i],atol=1e-9,rtol=0):
                raise ValueError("Original map statistic not in exact state-map null")
        if (i+1)%2500==0:
            print(f"EAzami_NA_U_D_MAPS_{tag}={i+1}/{count}",flush=True)
    if matches!=1:raise ValueError("Observed original author U/D map not unique")
    mu=score.mean(axis=0);sd=score.std(axis=0)
    if min(sd)<=0:raise ValueError("Degenerate EAzami null on actual NA nuclear tree")
    z=(original-mu)/sd
    z0=(score-mu)/sd
    maxT=np.max(abs(z0),axis=1)
    axes={}
    for j,var in enumerate(FOUR):
        critical=abs(z[j])
        rank=int(np.count_nonzero(abs(z0[:,j])>=critical-1e-12))
        maxrank=int(np.count_nonzero(maxT>=critical-1e-12))
        axes[var]={
          "signed_branch_environment_shift_standardized":float(original[j]),
          "direction":"higher" if original[j]>0 else "lower",
          "two_sided_state_map_exact_count":rank,
          "two_sided_state_map_exact_fraction":rank/count,
          "four_axis_maxT_exact_count":maxrank,
          "four_axis_maxT_exact_fraction":maxrank/count}
    bio18=axes[FOUR[0]]
    onesided=int(np.count_nonzero(score[:,0]>=original[0]-1e-12))
    bio18["BIO18_positive_U_to_D_state_map_exact_count"]=onesided
    bio18["BIO18_positive_U_to_D_state_map_exact_fraction"]=onesided/count
    return {
      "panel":tag,"n_authentic_NA_nuclear_species":len(names),
      "n_original_author_biological_D":d,"n_original_author_biological_U":len(names)-d,
      "original_independent_author_nodding_taxa":sorted(n for n in names if states[n]==1),
      "n_original_S1_four_climate_measured_images":sum(int(by[n]["n_photo_complete_four_axis"]) for n in names),
      "n_combinatorial_count_preserving_exact_maps":count,
      "all_possible_count_preserving_maps_refitted_original_CTMC_q":True,
      "observed_map_included_once":True,
      "observed_author_original_posterior_change_mass":expected,
      "observed_author_original_CTMC_q":q,
      "four_environment":axes,
      "four_axis_omnibus_maxT_exact_fraction":float(
          np.count_nonzero(maxT>=max(abs(z))-1e-12)/count),
      "no_photographic_head_angles_dichotomised":True,
      "species_tree_is_independent_published_multilocus_NA63":True,
      "not_old_EAzami_12_taxon_six_AU_topologies":True,
      "state_map_fraction_not_independent_biological_event_p":True,
      "ancestral_climate_not_historically_observed":True}

def run(original,newfour,ea_path,outdir):
    h1=ea.original_implementation(ea_path)
    treebytes=dna.upstream(dna.NA_REPO,dna.NA_REF,dna.NA_TREE)
    traits=dna.upstream(dna.NA_REPO,dna.NA_REF,NA_TRAITS)
    source,climate,states,tree,full,n10=inputs(original,newfour,treebytes,traits)
    out={
      "schema":"azami_S1_4axis_NA_genuine_multilocus_original_U_D_EAzami_exact_transition_null_v1",
      "status":"TRUE_ORIGINAL_NORTH_AMERICA_AUTHOR_UD_4CLIMATE_S1_SOURCE_CTMC_ALL_FINITE_MAPS",
      "original_S1_251_photo_angles_used_as_an_environmental_cohort_not_for_U_D_classification":True,
      "four_source_climates":list(FOUR),
      "original_S1_251_sha256":SOURCE_S1_SHA,
      "same_S1_S2_site_fourclimate_sha256":SOURCE_NEW4_SHA,
      "published_NA_author_tree_git_blob":dna.NA_GIT_BLOB,
      "published_NA_author_U_D_trait_git_blob":NA_TRAIT_BLOB,
      "EAzami_original_U_D_CTMP_algorithm_git_blob":ea.EA_CORE_GIT_BLOB,
      "source_photo_angle_is_NOT_true_gravity_U_D":True,
      "one_independent_author_nuclear_tree_not_52_phylogenetic_grafts":True,
      "two_source_populations_not_two_independent_genetic_replicates":True,
      "NA42_strict_full":exact_fit(h1,tree,climate,states,full,"NA42"),
      "NA40_at_least10_S1_photos":exact_fit(h1,tree,climate,states,n10,"NA40"),
    }
    folder=Path(outdir);folder.mkdir(parents=True,exist_ok=True)
    (folder/"EAzami_NA42_NA40_INDEPENDENT_TRUE_NUCLEAR_CTMC_NEW_FOUR_CLIMATE.json").write_text(
         json.dumps(out,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    for label in ("NA42_strict_full","NA40_at_least10_S1_photos"):
        x=out[label]
        print("INDEPENDENT_NA_EAzami_CTMC_SOURCE4="+json.dumps({
          "panel":label,"n":x["n_authentic_NA_nuclear_species"],
          "D":x["n_original_author_biological_D"],
          "all_maps":x["n_combinatorial_count_preserving_exact_maps"],
          "bio18_one_sided":x["four_environment"]["chelsa_bio18"]["BIO18_positive_U_to_D_state_map_exact_fraction"],
          "bio18_maxT":x["four_environment"]["chelsa_bio18"]["four_axis_maxT_exact_fraction"],
          "omnibus":x["four_axis_omnibus_maxT_exact_fraction"]
        },sort_keys=True),flush=True)
    return out

def main():
    p=argparse.ArgumentParser()
    for arg in ("original251","newfour251","EAzami-original-method","outdir"):
        p.add_argument("--"+arg,required=True,type=Path)
    x=p.parse_args()
    run(x.original251,x.newfour251,x.EAzami_original_method,x.outdir)
if __name__=="__main__":main()
