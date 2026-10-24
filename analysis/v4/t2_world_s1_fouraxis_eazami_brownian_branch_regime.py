"""World-S1 251 -> 191 exact placement tips -> nuclear ITS105: EAzami-style branch test.

The 2026-09-04 EAzami U->D test uses independent botanical binary states,
a 2-state CTMC and exhaustive count-fixed state maps. WORLD S1 contains only
continuous 0..180-degree photographic camera-axis angles, NOT validated
botanical U/D: do NOT binarize photographs or invent CTMC U->D gains.
Instead use EAzami's Brownian branch-length-weighted squared-change ancestor
model on both actual S1 original photo angle and actual CHELSA BIO18/BIO05,
RSDS, WorldClim elevation. Align signed inferred angle increments with
climatic increments, then permute all original continuous tip angles, keeping
one fixed phylogeny and all original climatic tips per run. Ensemble statistic
is the median across the 52 original, correlated GBOTB/V.PhyloMaker placement
SCENARIOS, not 52 independent evolutionary replicates.

Use 251 original photographic source names, 202 source-matching global tips,
191 excluding WCVP concept collisions, sensitivity n>=10, real 109 GenBank
ITS tips pruned to 105 exact alias-clean source names. ITS is a MULTICOPY
nuclear-rDNA gene tree, NOT a validated multilocus nuclear species phylogeny.
Neither global 52 placements nor ITS implies true inherited gravitropic
nodding, historical selection, phylogenetic causation or independent gains.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,copy
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from Bio import Phylo

ORIGINAL_SHA="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
FOUR_SHA="823ecb8442eb03598253cb33f3e16bac86538c6c2652f48ca2dc0886d0e19626"
ITS_SHA="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"
FOUR=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
RNG_SEED=20261024

def verified_csv(path,sha):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:
        raise ValueError("Original world S1 or measured four-site environment source SHA drift")
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))

def source_panels(original251,new4):
    source=verified_csv(original251,ORIGINAL_SHA)
    climate=verified_csv(new4,FOUR_SHA)
    names=[r["original_source_name"] for r in source]
    if len(names)!=251 or len(set(names))!=251 or sum(int(r["original_photo_count"]) for r in source)!=170066:
        raise ValueError("Original global 251/170066 native source photo universe changed")
    by={r["original_source_name"]:r for r in climate}
    if len(by)!=251 or set(by)!=set(names) or sum(int(r["n_photo_complete_four_axis"]) for r in climate)!=169020:
        raise ValueError("New true four-climate 169020 source photos / 251 names changed")
    matched=[];unambiguous=[]
    for item in source:
        name=item["original_source_name"];row=by[name]
        if item["accepted_key"]!=row["accepted_key"]:
            raise ValueError("World original photo and climate taxonomic concept drift")
        present=item["exact_source_name_has_frozen_placement_tip"]=="True"
        collide=item["WCVP_key_collides_source_names"]=="True"
        if item["exact_source_name_has_frozen_placement_tip"] not in ("True","False") or item["WCVP_key_collides_source_names"] not in ("True","False"):
            raise ValueError("Original source placement and WCVP alias status unknown")
        if int(row["n_photo_complete_four_axis"])>int(item["original_photo_count"]):
            raise ValueError("New four climate photos cannot exceed original eligible photos")
        if not 0<=float(row["mean_image_angle_deg"])<=180 or not all(np.isfinite(float(row[z])) for z in FOUR):
            raise ValueError("Source photographic continuous angle or climate invalid")
        if present:matched.append(name)
        if present and not collide:unambiguous.append(name)
    if len(matched)!=202 or len(unambiguous)!=191:
        raise ValueError("Original S1 source 202/191 placement coverage changed")
    return by,{"world_original202":sorted(matched),"world_collision_free191":sorted(unambiguous)},source

def pruned(tree,names,midpoint=False):
    t=copy.deepcopy(tree)
    available={x.name for x in t.get_terminals()}
    if not set(names)<=available:
        raise ValueError("Historical or ITS molecular tree lacks original S1 matched tip")
    for node in list(t.get_terminals()):
        if node.name not in names:t.prune(node)
    if midpoint:t.root_at_midpoint()
    if len(t.get_terminals())!=len(names):
        raise ValueError("Original world phylogeny has duplicate or collapsed source tips")
    pos=[float(node.branch_length) for node in t.find_clades()
         if node is not t.root and node.branch_length is not None and node.branch_length>0]
    if not pos:raise ValueError("Original tree has no valid positive branch lengths")
    eps=min(pos)*1e-7
    for node in t.find_clades():
        if node is t.root:continue
        if node.branch_length is None or not np.isfinite(node.branch_length) or node.branch_length<0:
            raise ValueError("Original molecular or placement tree negative/missing edge")
        if node.branch_length==0:node.branch_length=eps
    return t

def brownian_branch_operator(t,names):
    """Exactly the EAzami Luu inverse Brownian squared-change ancestral method.

    The conditional ancestor/terminal vector is a LINEAR function of original
    species angles. Compile the operator once per tree and apply the same
    999+ shuffled continuous source-tip vectors without any U/D threshold.
    """
    if len(names)<8 or len(set(names))!=len(names):
        raise ValueError("World 8+ biological source names are required")
    tips={node.name:node for node in t.get_terminals()}
    if set(tips)!=set(names):raise ValueError("True source-name taxon/tree mismatch")
    nodes=list(t.find_clades(order="preorder"))
    ix={node:j for j,node in enumerate(nodes)}
    n=len(names);k=len(nodes)
    L=np.zeros((k,k));edges=[]
    for parent in nodes:
        for child in parent.clades:
            a=ix[parent];b=ix[child];bl=max(float(child.branch_length or 0),1e-10)
            weight=1/bl
            L[a,a]+=weight;L[b,b]+=weight;L[a,b]-=weight;L[b,a]-=weight
            edges.append((a,b,child.is_terminal()))
    fixed=[ix[tips[name]] for name in names]
    free=[ix[node] for node in nodes if not node.is_terminal()]
    if len(fixed)+len(free)!=k:raise ValueError("Unknown non-tip nodes")
    F=np.zeros((k,n))
    F[fixed,np.arange(n)]=1
    F[free]=solve(L[np.ix_(free,free)],-L[np.ix_(free,fixed)],assume_a="pos")
    D=np.vstack([F[c]-F[p] for p,c,_ in edges])
    internal=np.array([not terminal for _,_,terminal in edges])
    if D.shape[1]!=n or internal.sum()<n//4 or not np.isfinite(D).all():
        raise ValueError("Original EAzami Brownian branch score incompatible with tree")
    if not np.allclose(D@np.ones(n),0,atol=1e-7):
        raise ValueError("Zero-angle-change constant ancestral trait was not preserved")
    return D,internal

def values(by,names):
    head=np.array([float(by[t]["mean_image_angle_deg"]) for t in names])
    env=np.array([[float(by[t][k]) for k in FOUR] for t in names])
    sd=env.std(axis=0)
    if not np.isfinite(env).all() or (sd<=0).any():
        raise ValueError("World original four-climate niches not all observed")
    return head,(env-env.mean(axis=0))/sd,sd

def branch_score(D,env,angles,mask):
    projection=D[mask]
    denv=projection@env
    dhead=angles@projection.T
    denom=np.abs(dhead).sum(axis=1)
    if np.any(denom<1e-7):raise ValueError("No inferred signed camera-angle branch change")
    # Same signed change x reconstructed climatic change structure as EAzami,
    # but unthresholded degrees rather than posterior binary state changes.
    return (dhead@denv)/denom[:,None]

def rank(obs,null):
    center=null.mean(axis=0);sd=null.std(axis=0)
    if obs.shape!=(4,) or null.ndim!=2 or null.shape[1]!=4 or (sd<=0).any():
        raise ValueError("Invalid four-axis continuous-tip permutation null")
    z=(obs-center)/sd;z0=(null-center)/sd
    total=len(null)+1
    p=[float((1+np.count_nonzero(np.abs(z0[:,j])>=abs(z[j])-1e-12))/total) for j in range(4)]
    maxnull=np.max(np.abs(z0),axis=1)
    maxcorr=[float((1+np.count_nonzero(maxnull>=abs(z[j])-1e-12))/total) for j in range(4)]
    return {
       "scores_standardised_climatic_shift_per_signed_degree":dict(zip(FOUR,obs.tolist())),
       "BIO18_positive_direction_one_sided_fraction":float((1+np.count_nonzero(null[:,0]>=obs[0]-1e-12))/total),
       "four_axis_two_sided_fraction":dict(zip(FOUR,p)),
       "four_axis_maxT_fraction":dict(zip(FOUR,maxcorr)),
       "four_axis_maxT_omnibus_fraction":float((1+np.count_nonzero(maxnull>=max(abs(z))-1e-12))/total),
       "finite_tip_permutation_not_count_preserving_U_D_CTMP":True,
       "conditional_null_not_species_tree_uncertainty":True,
    }

def evaluate_panel(by,names,trees,nperm,seed,edge_scope="all_edges"):
    if nperm<199:raise ValueError("Insufficient continuous-tip null resamples")
    names=sorted(names)
    y,env,sd=values(by,names)
    rng=np.random.default_rng(seed)
    shuffled=np.vstack([y]+[rng.permutation(y) for _ in range(nperm)])
    scenarios=[]
    audits=[]
    for name,tree in trees:
        current=pruned(tree,names,midpoint=name=="GenBank_ITS_midpoint")
        D,internal=brownian_branch_operator(current,names)
        mask=np.ones(len(D),dtype=bool) if edge_scope=="all_edges" else internal
        scores=branch_score(D,env,shuffled,mask)
        scenarios.append(scores)
        audits.append({"tree_name":name,"tree_edges":len(D),"internal_edges":int(internal.sum())})
    # All original 52 placement scenarios get the SAME tip permutation by
    # index; do not treat correlated placed topologies as 52 independent runs.
    combined=np.median(np.asarray(scenarios),axis=0)
    response=rank(combined[0],combined[1:])
    scores=np.array([x[0] for x in scenarios])
    response.update({
        "n_original_source_names":len(names),
        "n_original_fouraxis_head_photos":sum(int(by[x]["n_photo_complete_four_axis"]) for x in names),
        "n_original_tree_placement_scenarios":len(trees),
        "n_continuous_photo_tip_permutations":nperm,
        "per_scenario_BIO18_scores":[float(x[0]) for x in scores],
        "n_scenarios_BIO18_direction_positive":int(np.count_nonzero(scores[:,0]>0)),
        "BIO18_scenario_range":[float(np.min(scores[:,0])),float(np.max(scores[:,0]))],
        "tree_names":[x["tree_name"] for x in audits],
        "one_tree_audit":audits[0],
        "edge_scope":edge_scope,
        "within_panel_environment_sd":dict(zip(FOUR,sd.tolist())),
        "original_0_180_angle_not_gravity_nodding":True,
        "not_Mk_transition_model_or_botanical_U_D":True,
    })
    return response

def stable_seed(tag,seed=20261024):
    return seed+sum((i+1)*ord(c) for i,c in enumerate(tag))

def analyse(s1,climate,ITS,*,tree_provider=None,nperm=999,scenario_limit=52):
    by,original_sets,roster=source_panels(s1,climate)
    if tree_provider is None:
        from analysis.v4.native_phylogeny import restore_trees
        tree_provider=restore_trees
    original_trees=tree_provider()
    if len(original_trees)!=52 or not 1<=scenario_limit<=52:
        raise ValueError("Original SHA-pinned global 52 placements unavailable")
    raw=Path(ITS).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=ITS_SHA:
        raise ValueError("Real original GenBank nuclear ITS source tree changed")
    its=Phylo.read(io.StringIO(raw.decode()),"newick")
    true_tip_names=sorted(x.name for x in its.get_terminals())
    if len(true_tip_names)!=109 or not set(true_tip_names)<=set(by):
        raise ValueError("109 original world S1 ITS tree tip names changed")
    rows={r["original_source_name"]:r for r in roster}
    its_clean=[name for name in true_tip_names
        if rows[name]["WCVP_key_collides_source_names"]=="False"]
    if len(its_clean)!=105:raise ValueError("Unambiguous real ITS 105 source labels changed")
    world=original_sets["world_collision_free191"]
    world10=[x for x in world if int(by[x]["n_photo_complete_four_axis"])>=10]
    common97=sorted(set(world)&set(its_clean))
    if len(world10)!=143 or len(common97)!=97:
        raise ValueError("World photo >=10 or matched genuine ITS 97 cohort changed")
    panels={
      "world191_original_52_grafted_placement_scenarios":
        evaluate_panel(by,world,original_trees[:scenario_limit],nperm,stable_seed("world191")),
      "world191_internal_edges_only":
        evaluate_panel(by,world,original_trees[:scenario_limit],nperm,stable_seed("world191"),edge_scope="internal_edges"),
      "world143_photo_min10_placement_sensitivity":
        evaluate_panel(by,world10,original_trees[:scenario_limit],nperm,stable_seed("world143")),
      "real_genbank_nuclear_ITS105_multicopy_gene_tree":
        evaluate_panel(by,its_clean,[("GenBank_ITS_midpoint",its)],nperm,stable_seed("its105")),
      "matched97_global52_grafted_placements":
        evaluate_panel(by,common97,original_trees[:scenario_limit],nperm,stable_seed("shared97")),
      "matched97_real_GenBank_ITS_gene_tree":
        evaluate_panel(by,common97,[("GenBank_ITS_midpoint",its)],nperm,stable_seed("shared97"))
    }
    return {"schema":"AZAMI_WORLD_S1_BM_EAZAMI_STYLE_FOUR_AXIS_BRANCH_TRAIT_ENV_NICHE_20261024_V1",
      "status":"WORLD_SOURCE_CONTINUOUS_PHOTO_ANGLE_BROWNIAN_BRANCH_NULL_COMPUTED_NOT_TRUE_U_D",
      "original_world_S1_251_source_names":251,
      "source_original_total_native_head_images":170066,
      "source_new_four_climate_all_complete_head_images":169020,
      "matched_52_grafted_world_source_names_including_aliases":202,
      "primary_world_source_names_collision_free":191,
      "world_source_names_photo_at_least10":143,
      "real_ITS_gene_tree_source_names":109,
      "real_ITS_gene_tree_source_names_alias_clean":105,
      "same_original_world_placement_and_ITS_strict_taxa":len(common97),
      "exposures":list(FOUR),
      "positive_direction_for_focal_BIO18_fixed_from_original_S1_species_signal":True,
      "no_posthoc_favourable_signs_for_BIO05_RSDS_elevation":True,
      "why_NOT_original_EAzami_CTMC":"Original EAzami has independent genuine botanical binary tip U/D and count-fixed discrete CTMC maps; global source has continuous photo-frame head angle and no independently calibrated categorical posture for all 251, so only its Brownian signed-branch-change / tip permutation analogy is scientifically legitimate.",
      "52_placements_not_true_global_multilocus_species_tree":True,
      "ITS_is_multicopy_gene_tree_not_true_nuclear_multilocus_species_tree":True,
      "not_a_true_gravitational_nodding_or_U_to_D_origin_test":True,
      "no_historical_climate_selection_or_adaptation_inferred":True,
      "source_sha256":{"S1_251":ORIGINAL_SHA,"BIO18_BIO05_RSDS_ELEV_251":FOUR_SHA,"ITS109":ITS_SHA},
      "panels":panels}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("original251","newfour","its-tree","out"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--permutations",type=int,default=999)
    p.add_argument("--scenario-limit",type=int,default=52)
    args=p.parse_args()
    r=analyse(args.original251,args.newfour,args.its_tree,nperm=args.permutations,scenario_limit=args.scenario_limit)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    for key,x in r["panels"].items():
        print("WORLD_S1_BRANCH_FOUR="+json.dumps({"panel":key,"n":x["n_original_source_names"],"score_BIO18":x["scores_standardised_climatic_shift_per_signed_degree"]["chelsa_bio18"],"BIO18_one_sided_tip_permutation_fraction":x["BIO18_positive_direction_one_sided_fraction"],"BIO18_4axis_maxT_fraction":x["four_axis_maxT_fraction"]["chelsa_bio18"],"n_positive_52_placement_scenarios":x["n_scenarios_BIO18_direction_positive"]},sort_keys=True))
if __name__=="__main__":
    main()
