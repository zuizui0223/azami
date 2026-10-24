"""EAzami ORIGINAL U->D CTMC true NA nuclear phenotype: all internal edges & every nodding-taxon LOO.

Strict reanalysis of an observed surprising NA42 D3 U39 positive BIO18
transition-climate alignment, using SAME actual published NA63 original
author U/D labels, NA63 bona fide multilocus nuclear species tree, S1/S2
original four site-climate environmental means, pinned original EAzami
symmetric CTMC and Brownian internal climate model.

Before claiming independently repeated environment-aligned U->D evolution,
verify the old EAzami internal-edge-only estimator (internal edge-only
normalization) and remove each of the three author-labeled D taxa in turn.
Exhaust every count-preserving tip state map (NA42 D3 = 11480, each LOO
NA41 D2 = 820). Refit symmetric CTMC q for EVERY assignment, not
unconditional or uncorrelated randomization; preserve axis maxT correction.
This cannot prove independent historical gains, climate at transition time,
fitness, or that S1 camera-angle itself is a botanical nodding state.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math
from pathlib import Path
import numpy as np
from analysis.v4 import t2_na_author_UD_eazami_true_nuclear_fourclimate_exact as na
from analysis.v4 import eazami_ud_real_four_environment_transition_regime as ea
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as dna

BASELINE_SHA="c97756b6806128aefdca7d09701474c7d21c0bd4ffd6afb9ce27f68369ecba03"
FOUR=na.FOUR

def internal_score(h1,tr,states,env):
    q=h1.fit_symmetric_q(tr,states)
    L,joint=h1.ctmc_likelihood_and_edge_joint(tr,states,q,need_joint=True)
    if joint is None or L<=0:raise ValueError("No original EAzami CTMC edge joint")
    _,children,_,_,_,_,_=h1.tree_structure(tr)
    indices=[i for i in joint if children[i]]
    den=sum(float(joint[i][0,1]+joint[i][1,0]) for i in indices)
    if den<=0 or not indices:
        raise ValueError("No internal EAzami phylogenetic branch change mass")
    out=np.array([sum((float(joint[i][0,1]-joint[i][1,0]))*env[v][i]
                    for i in indices)/den for v in FOUR])
    return out,float(den)

def finite_count_null(h1,original_tree,by,author,names,*,mode,tag):
    names=sorted(names)
    D=sum(int(author[n]) for n in names)
    size=math.comb(len(names),D)
    if (len(names),D,size) not in ((42,3,11480),(41,2,820)):
        raise ValueError("Published original author D-only deletion sample altered")
    t=h1.prune_to_panel(original_tree,names)
    env=na.environmental_edges(h1,t,names,by)
    observed_states={n:int(author[n]) for n in names}
    def score(st):
        if mode=="internal_only":
            return internal_score(h1,t,st,env)[0]
        elif mode=="all_edges":
            return na.transition_score(h1,t,st,env)[0]
        else:raise ValueError("Invalid EAzami branch scope")
    measured=score(observed_states)
    sims=np.empty((size,len(FOUR)))
    witness=0
    for i,comb in enumerate(itertools.combinations(range(len(names)),D)):
        mark=set(comb)
        st={n:int(j in mark) for j,n in enumerate(names)}
        sims[i]=score(st)
        if st==observed_states:
            witness+=1
            if not np.allclose(measured,sims[i],rtol=0,atol=1e-8):
                raise ValueError("Original state assignment not included in exact null")
    if witness!=1:raise ValueError("Original author U/D assignment is not unique")
    center=sims.mean(0);sds=sims.std(0)
    if (sds<=0).any():raise ValueError("Unidentifiable EAzami exact null")
    z=(measured-center)/sds;z0=(sims-center)/sds
    maxnull=np.max(abs(z0),axis=1)
    effects={}
    for axis,field in enumerate(FOUR):
        twoside=int(np.count_nonzero(abs(z0[:,axis])>=abs(z[axis])-1e-12))
        family=int(np.count_nonzero(maxnull>=abs(z[axis])-1e-12))
        effects[field]={"score":float(measured[axis]),
          "two_sided_count_preserved_rank":twoside,
          "two_sided_count_preserved_fraction":twoside/size,
          "four_axis_maxT_count_preserved_fraction":family/size}
    positive=int(np.count_nonzero(sims[:,0]>=measured[0]-1e-12))
    effects[FOUR[0]]["BIO18_positive_one_sided_exact_fraction"]=positive/size
    effects[FOUR[0]]["BIO18_positive_one_sided_exact_count"]=positive
    return {
      "tag":tag,"scope":mode,"species":len(names),"U":len(names)-D,"D":D,
      "n_exact_count_preserving_configurations":size,
      "observed_map_exhaustively_present_once":True,
      "refit_EAzami_q_for_each_state_map":True,
      "n_author_original_S1_fourclimate_photos":sum(
         int(by[n]["n_photo_complete_four_axis"]) for n in names),
      "scores_and_finite_map_ranks":effects,
      "four_axis_omnibus_maxT_fraction":
          float(np.count_nonzero(maxnull>=max(abs(z))-1e-12)/size),
      "S1_camera_angle_NOT_classified_as_U_D":True,
      "not_a_likelihood_test_over_independent_evolutionary_events":True,
    }

def run(original,newfour,ea_core,baseline,outfile):
    source=Path(baseline).read_bytes()
    if hashlib.sha256(source).hexdigest()!=BASELINE_SHA:
        raise ValueError("Actual completed NA42 D3 original finite CTMC baseline source drift")
    prev=json.loads(source)
    original42=prev["NA42_strict_full"]
    if (original42["n_combinatorial_count_preserving_exact_maps"]!=11480 or
        original42["n_original_author_biological_D"]!=3 or
        original42["four_environment"]["chelsa_bio18"]["BIO18_positive_U_to_D_state_map_exact_count"]!=63):
        raise ValueError("Original independently published source NA42 H1 score receipt changed")
    h1=ea.original_implementation(ea_core)
    treebytes=dna.upstream(dna.NA_REPO,dna.NA_REF,dna.NA_TREE)
    traitbytes=dna.upstream(dna.NA_REPO,dna.NA_REF,na.NA_TRAITS)
    orig,by,author,tree,NA42,NA40=na.inputs(original,newfour,treebytes,traitbytes)
    D=sorted(n for n in NA42 if author[n]==1)
    if D!=["Cirsium_conspicuum","Cirsium_funkiae","Cirsium_scopulorum"]:
        raise ValueError("Original independent NA nodding species source changed")
    tests={
      "NA42_only_internal_edges":
        finite_count_null(h1,tree,by,author,NA42,mode="internal_only",tag="NA42_internal_only")
    }
    for gone in D:
        new=[n for n in NA42 if n!=gone]
        tests["NA41_drop_"+gone]=finite_count_null(
            h1,tree,by,author,new,mode="all_edges",tag="NA41_drop_"+gone)
    result={
      "schema":"V4_original_S1_fourclimate_NA42_true_nuclear_original_EAzami_CTMC_internal_D_leaveoneout_v1",
      "status":"COMPLETED_INDEPENDENT_TRUE_NUCLEAR_ORIGINAL_AUTHOR_UD_CTMP_INTERNAL_AND_ALL_D_LOO",
      "original_NA42_full_measured_rank":
         original42["four_environment"]["chelsa_bio18"]["BIO18_positive_U_to_D_state_map_exact_fraction"],
      "original_NA42_full_measured_four_axis_maxT":
         original42["four_environment"]["chelsa_bio18"]["four_axis_maxT_exact_fraction"],
      "original_baseline_json_sha256":BASELINE_SHA,
      "original_multilocus_63_tip_NA_author_tree_SHA_Git":dna.NA_GIT_BLOB,
      "original_author_botanical_U_D_trait_git_blob":na.NA_TRAIT_BLOB,
      "original_EAzami_CTMC_Brownian_algorithm_pin":ea.EA_CORE_GIT_BLOB,
      "four_environment_alternative_hypotheses":list(FOUR),
      "all_original_D_species_exhaustively_one_by_one_deleted":D,
      "no_S1_photo_angle_dichotomization":True,
      "no_climate_at_ancestor_observed":True,
      "tests":tests,
    }
    path=Path(outfile);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    for tag,rec in tests.items():
        print("TRUE_NA_EAzami_CTMC_ROBUSTNESS="+json.dumps({
         "test":tag,"n_species":rec["species"],"D":rec["D"],
         "exact_maps":rec["n_exact_count_preserving_configurations"],
         "bio18_score":rec["scores_and_finite_map_ranks"]["chelsa_bio18"]["score"],
         "bio18_positive_one_side":rec["scores_and_finite_map_ranks"]["chelsa_bio18"]["BIO18_positive_one_sided_exact_fraction"],
         "BIO18_fourclimate_maxT":rec["scores_and_finite_map_ranks"]["chelsa_bio18"]["four_axis_maxT_count_preserved_fraction"]
        },sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser()
    for name in ("original251","newfour251","EAzami-method","NA42-baseline","out"):
        p.add_argument("--"+name,required=True,type=Path)
    z=p.parse_args()
    run(z.original251,z.newfour251,z.EAzami_method,z.NA42_baseline,z.out)
if __name__=="__main__":main()
