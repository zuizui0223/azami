"""T2: does the ORIGINAL S1 continuous photo head angle cluster on genuine trees?

Source-pinned measured original 251 taxa, 170066 photo observations;
published NA63 nuclear tree pruned to exactly non-collided NA42;
genuine newly inferred GenBank nuclear-rDNA ITS109 gene tree; posthoc110
Q20 108-candidate-reference-locus gene tree (original111 QC FAILED);
both gene trees pruned separately to the EXACT SAME 58 photographic taxa.

For each unchanged 0-180 S1 photo angle (mean + median) and log(BIO12),
compare average closest-genetic-tip absolute differences with 9999
permutations of tip values on the SAME molecular topology. Lower one-sided
randomization P = nearer relatives are MORE similar than arbitrary label
placement. This is a descriptive topological signal, NOT causal history,
adaptive gains, a gravity-validated head tilt or a species-level trait model.
Visualization roots ITS/Q20 at midpoint for display only; patristic nearest
neighbor distances are unrooted-invariant. Branch lengths between sources
are NEVER compared numerically.
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import io

import numpy as np
from Bio import Phylo
from scipy.stats import false_discovery_control

N_S1=251
PHOTO_S1=170066
S1_SHA="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
ITS_SHA="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"
Q20_SHA="541f21801d79a7dda872a020835e6838ea0fd002668324a157d24bca16b9fa08"
NA63_GIT_BLOB_SHA="6417a38793427bc7c8562dc0e2cab59b74d8c316"
NA63_AUTHOR_COMMIT="c925d4af86160064448f794a9715f70b917a2786"
TREE_COHS={"NA42_real_authored_nuclear":42,"ITS109_GenBank_nuclear_rDNA":109,
           "Q20_posthoc110_candidate_loci":110,
           "same58_ITS_gene_tree":58,"same58_Q20_candidate_tree":58}
MIN_PHOTOS=(1,10,30)
ANGLE_ENDPOINTS=("mean_image_angle_deg","median_image_angle_deg")
SEED=20261024

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_sources(s1, na63, q20_source, q20_receipt, q20_tree, its_tree):
    from analysis.v4.t2_110_multilocus_vs_109_its_photo_niche import verify_sources
    original, qt, it, mapping, shared, _ = verify_sources(
        s1,q20_source,q20_receipt,q20_tree,its_tree)
    if len(original)!=N_S1 or sum(int(r["original_photo_count"]) for r in original.values())!=PHOTO_S1:
        raise ValueError("S1 original photo census or taxon names drift")
    na_bytes=Path(na63).read_bytes()
    git_blob=hashlib.sha1(b"blob "+str(len(na_bytes)).encode()+b"\0"+na_bytes).hexdigest()
    if git_blob!=NA63_GIT_BLOB_SHA:
        raise ValueError("Original published NA nuclear tree Git blob drift")
    na=Phylo.read(io.StringIO(na_bytes.decode("utf-8-sig")),"newick")
    if len(na.get_terminals())!=63:
        raise ValueError("Published North American nuclear phylogeny not 63 tips")
    na_map={x.name:x.name for x in na.get_terminals()}
    na_valid={x for x in na_map.values() if x in original and
              original[x]["WCVP_key_collides_source_names"]=="False"}
    if len(na_valid)!=42:
        raise ValueError("Must exclude both genuine original WCVP name collisions to get NA42")
    it_map={x.name:x.name for x in it.get_terminals()}
    if set(it_map)!=set(it_map)&set(original) or len(it_map)!=109:
        raise ValueError("Original GenBank ITS109 lacks exact source name matching")
    if len(shared)!=58 or set(mapping.values())&set(it_map)!=shared:
        raise ValueError("Real 58 source-name intersection changed")
    return original,[
        ("NA42_real_authored_nuclear",na,na_map,na_valid,"original_published_63tip_nuclear_pruned_to42"),
        ("ITS109_GenBank_nuclear_rDNA",it,it_map,set(it_map),"original_ITS109_multicopy_linked_marker_gene_tree"),
        ("Q20_posthoc110_candidate_loci",qt,mapping,set(mapping.values()),
         "original111_FAILED_technical_gate_posthoc110_reference_copy_gene_tree"),
        ("same58_ITS_gene_tree",it,it_map,set(shared),"original_ITS109_pruned_to_exact_same58"),
        ("same58_Q20_candidate_tree",qt,mapping,set(shared),
         "original_Q20_posthoc110_pruned_to_exact_same58"),
    ]

def sub_tree(tree, mapping, include):
    clone=deepcopy(tree)
    if len(include)!=len(set(include)):
        raise ValueError("Duplicate original source taxa")
    keep={k for k,v in mapping.items() if v in include}
    if len(keep)!=len(include):
        raise ValueError("Incomplete molecular source-name mapping")
    for leaf in list(clone.get_terminals()):
        if leaf.name not in keep:
            clone.prune(leaf)
    names=[mapping[x.name] for x in clone.get_terminals()]
    if len(names)!=len(include) or len(names)!=len(set(names)):
        raise ValueError("One S1 photographic source taxon required per genuine gene tip")
    return clone,names

def nearest_phylo_neighbor(tree):
    tips=tree.get_terminals()
    if len(tips)<7:
        raise ValueError("Too few true molecular tips")
    D=np.array([[float(tree.distance(a,b)) for b in tips] for a in tips])
    if not np.isfinite(D).all() or np.any(D<-1e-9):
        raise ValueError("Invalid actual DNA branch distance")
    np.fill_diagonal(D,np.inf)
    return D.argmin(axis=1),D

def trait_clustering(values, nearest, label, *, reps=9999, seed=SEED):
    """Exact data, permutation null; use original scalar units throughout."""
    vals=np.asarray(values,float)
    if vals.ndim!=1 or not np.isfinite(vals).all() or np.std(vals)<=0:
        raise ValueError("No valid continuous original trait")
    n=len(vals)
    if len(nearest)!=n or len(set(range(n)) - set(nearest))==n:
        raise ValueError("Nearest relatives not defined")
    obs=float(np.mean(np.abs(vals-vals[nearest])))
    key=hashlib.sha256(label.encode()).hexdigest()[:8]
    rng=np.random.default_rng(seed+int(key,16))
    draws=np.empty(reps,float)
    for start in range(0,reps,250):
        count=min(250,reps-start)
        order=np.array([rng.permutation(n) for _ in range(count)])
        shuffled=vals[order]
        draws[start:start+count]=np.mean(
            np.abs(shuffled-shuffled[:,nearest]),axis=1)
    expect=float(np.mean(draws))
    p_low=float((1+np.count_nonzero(draws<=obs))/(reps+1))
    return {
        "n_taxa":n,
        "value_difference_units":("degrees_original_camera_photo_angle" if label.endswith(("mean_image_angle_deg","median_image_angle_deg")) else "absolute_log_annual_precipitation_difference"),
        "observed_closest_tip_mean_absolute_difference":obs,
        "unrelated_tip_label_permutation_mean":expect,
        "observed_to_null_ratio_below1_indicates_clustering":obs/expect,
        "permutation_null_95_interval":[float(x) for x in np.quantile(draws,[.025,.975])],
        "single_sided_clustering_p":p_low,
        "n_permutations":reps,
        "null":"source_tip_label_shuffle_on_exact_fixed_actual_sequence_tree",
        "topology_branch_error_not_propagated":True,
        "not_a_test_of_rain_causal_evolution_or_directional_nodding":True,
    }

def calculate_tree(source, cohort, min_photos, *, reps,seed):
    label,tree,mapping,allowed,kind=cohort
    sample={name for name in allowed if
        int(source[name]["original_photo_count"])>=min_photos}
    sub,names=sub_tree(tree,mapping,sample)
    nearest,D=nearest_phylo_neighbor(sub)
    analyses={}
    for endpoint in ANGLE_ENDPOINTS:
        values=[float(source[name][endpoint]) for name in names]
        analyses[endpoint]=trait_clustering(
            values,nearest,label+"-"+str(min_photos)+"-"+endpoint,
            reps=reps,seed=seed)
    analyses["log_BIO12"]=trait_clustering(
        np.log([float(source[name]["chelsa_bio12"]) for name in names]),
        nearest,label+"-"+str(min_photos)+"-log_BIO12",reps=reps,seed=seed)
    pairs=[]
    seen=set()
    for i,j in enumerate(nearest):
        key=tuple(sorted((names[i],names[int(j)])))
        if key in seen: continue
        seen.add(key)
        a,b=key
        pairs.append({
            "tip_A":a,"tip_B":b,
            "DNA_patristic_pair_distance_in_original_tree_units":float(D[i,int(j)]),
            "S1_image_angle_mean_A_deg":float(source[a]["mean_image_angle_deg"]),
            "S1_image_angle_mean_B_deg":float(source[b]["mean_image_angle_deg"]),
            "absolute_angle_mean_difference_deg":
                abs(float(source[a]["mean_image_angle_deg"])-float(source[b]["mean_image_angle_deg"])),
            "BIO12_A_mm":float(source[a]["chelsa_bio12"]),
            "BIO12_B_mm":float(source[b]["chelsa_bio12"]),
            "absolute_log_BIO12_difference":
                abs(np.log(float(source[a]["chelsa_bio12"])/float(source[b]["chelsa_bio12"]))),
        })
    # This summary intentionally chooses pairs by genetic patristic distance,
    # NOT by a desirable head-angle or precipitation result.
    pairs.sort(key=lambda z:(z["DNA_patristic_pair_distance_in_original_tree_units"],
                             z["tip_A"],z["tip_B"]))
    return {
        "n_real_tree_tips_matched_to_original_S1":len(names),
        "n_original_S1_photos_represented":
            sum(int(source[x]["original_photo_count"]) for x in names),
        "min_S1_original_photo_observations_per_taxon":min_photos,
        "tree_source_type":kind,
        "three_original_S1_traits":analyses,
        "first_12_closest_genetic_neighbor_pairs_by_DNA_branch_distance":pairs[:12],
        "not_biological_ancestor_tilt_reconstruction":True,
    },(sub,names,nearest,pairs)

def make_figure(out, source, tree, names, title, full_note):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    shown=deepcopy(tree)
    if "NA42" not in title:
        shown.root_at_midpoint()  # DISPLAY ONLY; pairwise distance unchanged
    tips=shown.get_terminals()
    y={tip:(i+1) for i,tip in enumerate(tips)}
    def assign(clade):
        if clade in y: return y[clade]
        vals=[assign(x) for x in clade.clades]
        y[clade]=float(np.mean(vals))
        return y[clade]
    assign(shown.root)
    dep=shown.depths()
    maxdep=max(float(dep[t]) for t in tips)
    if maxdep<=0: raise ValueError("No molecular branch variation")
    n=len(tips)
    fig,ax=plt.subplots(figsize=(12.0,max(7,n*.22)))
    def draw(parent):
        x=dep[parent]/maxdep
        if parent.clades:
            vals=[y[k] for k in parent.clades]
            ax.plot([x,x],[min(vals),max(vals)],lw=.8,color="0.40")
            for ch in parent.clades:
                child=dep[ch]/maxdep
                ax.plot([x,child],[y[ch],y[ch]],lw=.8,color="0.40")
                draw(ch)
    draw(shown.root)
    a_norm=Normalize(0,180)
    prec_all=np.array([float(x["chelsa_bio12"]) for x in source.values()])
    b_norm=LogNorm(vmin=max(1,float(prec_all.min())),
                   vmax=float(prec_all.max()))
    name_map={x.name:x.name for x in tips} if "Q20" not in title else {
        x.name:full_note[x.name] for x in tips}
    angle_cmap=matplotlib.colormaps["cividis"]
    prec_cmap=matplotlib.colormaps["viridis"]
    for tip in tips:
        name=name_map[tip.name]
        if name not in source:
            raise ValueError("Original DNA tip not matched to S1 photo source")
        yy=y[tip]
        ax.text(1.03,yy,name.replace("Cirsium_","C. ").replace("_"," "),
                fontsize=6.0,va="center")
        ax.add_patch(Rectangle((1.83,yy-.41),.14,.82,
            facecolor=angle_cmap(a_norm(float(source[name]["mean_image_angle_deg"]))),
            edgecolor="none"))
        ax.add_patch(Rectangle((2.02,yy-.41),.14,.82,
            facecolor=prec_cmap(b_norm(float(source[name]["chelsa_bio12"]))),
            edgecolor="none"))
        ax.text(2.23,yy,str(source[name]["original_photo_count"]),
                fontsize=5.7,va="center")
    ax.text(1.9,0,"ANGLE",fontsize=7,ha="center",weight="bold")
    ax.text(2.09,0,"BIO12",fontsize=7,ha="center",weight="bold")
    ax.text(2.27,0,"n",fontsize=7,ha="center",weight="bold")
    ax.set_xlim(-.015,2.41)
    ax.set_ylim(n+1,.3)
    ax.set_axis_off()
    ax.set_title(title,fontsize=12,loc="left",pad=24)
    ax.text(.01,1.005,"Original S1 photographic angle / annual precipitation | DNA branch length within panel only",
            fontsize=7,transform=ax.transAxes)
    fig.savefig(out,dpi=160,bbox_inches="tight")
    plt.close(fig)

def run(s1,na63,q20_source,receipt,q20_tree,its_tree,out,*,reps=9999,seed=SEED):
    source,cohorts=load_sources(s1,na63,q20_source,receipt,q20_tree,its_tree)
    output=Path(out);output.mkdir(parents=True,exist_ok=True)
    full={}
    figs={}
    for cohort in cohorts:
        name=cohort[0]
        reports={}
        first_display=None
        for gate in MIN_PHOTOS:
            report,display=calculate_tree(source,cohort,gate,reps=reps,seed=seed)
            reports["min_n"+str(gate)]=report
            if gate==1:first_display=display
        full[name]=reports
        if name in ("NA42_real_authored_nuclear","ITS109_GenBank_nuclear_rDNA",
                    "Q20_posthoc110_candidate_loci"):
            path=output/("T2_PHYLO_S1_ANGLE_BIO12_"+name+".png")
            t,names,_,_=first_display
            make_figure(path,source,t,names,
                name+" | "+str(len(names))+" actual S1-photo nuclear tips",
                cohort[2])
            figs[name]=path.name
    # Five overlapping molecular panels, not five independent new studies.
    names=list(TREE_COHS)
    for trait in ("mean_image_angle_deg","log_BIO12"):
        pv=[full[n]["min_n1"]["three_original_S1_traits"][trait]["single_sided_clustering_p"]
            for n in names]
        q=false_discovery_control(np.array(pv),method="bh")
        for n,qq in zip(names,q):
            full[n]["min_n1"]["three_original_S1_traits"][trait][
                "BH_q_5_nonindependent_tree_panels_exploratory"]=float(qq)
    result={
        "schema":"v4_T2_original_S1_continuous_true_sequence_phylogenetic_clustering_v1",
        "status":"REAL_GENE_AND_NUCLEAR_PHYLOGENETIC_NEIGHBOR_CLUSTER_PERMUTATION_COMPLETE",
        "S1_original_source_taxa":N_S1,
        "S1_original_photo_count":PHOTO_S1,
        "primary_response":"S1_original_uncategorized_continuous_head_image_angle_mean_degrees",
        "original_image_angle_not_gravity_verified":True,
        "original_botanical_U_D_not_analyzed":True,
        "NA63_original_nuclear_published_Git_blob":NA63_GIT_BLOB_SHA,
        "NA42_excluded_source_WCVP_concept_collision":True,
        "GenBank_ITS_gene_tree_not_accepted_species_tree":True,
        "Q20_110_from_failed_original111_primary_gate_NOT_species_tree":True,
        "161_source_name_union_NOT_a_161tip_tree":True,
        "identical_58_taxon_set_both_real_geneaologies":True,
        "min_taxon_original_image_eligibility_gates":list(MIN_PHOTOS),
        "permutation_seed":seed,
        "permutation_trials_per_test":reps,
        "five_nonindependent_molecular_tree_panels":full,
        "published_original_S1_data_sha256":S1_SHA,
        "authentic_ITS109_tree_sha256":ITS_SHA,
        "posthoc_Q20_tree_sha256":Q20_SHA,
        "tree_visuals":figs,
        "main_interpretation_ceiling":{
            "phylogenetic_clustering_is_tip_similarity_not_ancestral_adaptive_transition":True,
            "fixed_source_tree_topology_branch_no_bootstrap_uncertainty":True,
            "species_environmental_geography_not_independent_in_tip_shuffle":True,
            "photo_source_and_gene_source_sampling_are_nonrandom":True,
            "n_taxa_not_n_photographs_is_phylogenetic_replication":True,
        }
    }
    (output/"T2_S1_CONTINUOUS_REAL_TREE_PHYLOGENETIC_CLUSTERING.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    for name in names:
        r=full[name]["min_n1"]
        a=r["three_original_S1_traits"]["mean_image_angle_deg"]
        b=r["three_original_S1_traits"]["log_BIO12"]
        print(f"T2_PHOTOGRAPHIC_PHYLOGENY {name} n={r['n_real_tree_tips_matched_to_original_S1']} "
              f"angle_ratio={a['observed_to_null_ratio_below1_indicates_clustering']:.3f} "
              f"angle_p={a['single_sided_clustering_p']:.4f} "
              f"rain_ratio={b['observed_to_null_ratio_below1_indicates_clustering']:.3f} "
              f"rain_p={b['single_sided_clustering_p']:.4f}",flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("s1","NA63-original-tree","q20-source110","receipt","Q20-tree","ITS109-tree","out"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--permutations",type=int,default=9999)
    a=p.parse_args()
    if a.permutations<999:
        raise ValueError("At least 999 random permutations per genuine tree endpoint")
    run(a.s1,a.NA63_original_tree,a.q20_source110,a.receipt,a.Q20_tree,
        a.ITS109_tree,a.out,reps=a.permutations)
if __name__=="__main__":main()
