"""S1 continuous head PHOTO angle vs BIO12 on TWO authentic nuclear trees.

Primary genuine multilocus molecular trees (NOT chloroplast, ITS109 or posthoc
Q20 candidate110): Siniscalchi 2023 NA63 nuclear species tree and EAzami
Japan38 176-rootable-locus nuclear reconstruction. Match S1 ALL original251
source binomials with no accepted-key collisions or infraspecific collapse.
The outcome is never a botanical nodding/erect category. Do not combine
unrelated trees into one global tree or count evolutionary U->D events.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import csv
import hashlib
import io
import json
import re
import urllib.request
from pathlib import Path
import numpy as np
from Bio import Phylo

S1_SHA = "47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
NA_REPO = "carol-siniscalchi/cirsium"
NA_REF = "c925d4af86160064448f794a9715f70b917a2786"
NA_TREE = "traits/cirsium.phyluce95p.ultrametric.tre"
NA_GIT_BLOB = "6417a38793427bc7c8562dc0e2cab59b74d8c316"
J_REPO = "zuizui0223/EAzami"
J_REF = "c9870e91604ed1514344ddc8f765eae5237dbf27"
J_TREE = "data/evidence/japan38_comp1061_primary_tree_v1.nwk"
J_MAP = "data/evidence/japan38_comp1061_concept_map_v1.csv"
J_MANIFEST = "data/evidence/japan38_comp1061_primary_tree_acceptance_v1.json"
J_TREE_SHA = "ba592d3a9f52b5481314ab20e3202048505cdf72a3fc15134d5dca261b9d9b7e"
NA_EXPECTED = 42
J_EXPECTED = 12
SEED = 20261024
B = 9999
TRAITS = ("mean_image_angle_deg","median_image_angle_deg","log_BIO12")
PHOTO_GATES = (1,10,30)

def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def upstream(repo,sha,path):
    url=f"https://raw.githubusercontent.com/{repo}/{sha}/{path}"
    with urllib.request.urlopen(url,timeout=90) as response:
        return response.read()

def read_sources(original251, na_raw, j_raw, j_map_raw, j_manifest_raw):
    original=Path(original251).read_bytes()
    if hashlib.sha256(original).hexdigest()!=S1_SHA:
        raise ValueError("S1 real measured original 251 source SHA mismatch")
    rows=list(csv.DictReader(io.StringIO(original.decode("utf-8-sig"))))
    if len(rows)!=251 or len({r["accepted_key"] for r in rows})!=244 or sum(
        int(r["original_photo_count"]) for r in rows)!=170066:
        raise ValueError("S1 original census 251/244/170066 altered")
    if any(r["botanical_gravity_posture_inferred_from_photo_angle"]!="False" for r in rows):
        raise ValueError("S1 camera angle must not be assumed botanical gravity tilt")
    by={r["original_source_name"]:r for r in rows}
    if len(by)!=251:raise ValueError("S1 exact original source taxa duplicated")
    if git_blob(na_raw)!=NA_GIT_BLOB:raise ValueError("Published authentic NA nuclear tree Git blob drift")
    na=Phylo.read(io.StringIO(na_raw.decode()),"newick")
    if len(na.get_terminals())!=63:raise ValueError("Wrong NA published nuclear tree")
    if hashlib.sha256(j_raw).hexdigest()!=J_TREE_SHA:
        raise ValueError("Japan published nuclear tree SHA drift")
    jmanifest=json.loads(j_manifest_raw)
    if (jmanifest["tree_sha256"]!=J_TREE_SHA
        or jmanifest["tree_tip_count"]!=40
        or jmanifest["rootable_loci"]!=176
        or jmanifest["root_tip"]!="OUTGROUP_saff"):
        raise ValueError("Japan genuine multi-locus nuclear phylogeny provenance drift")
    jtree=Phylo.read(io.StringIO(j_raw.decode()),"newick")
    if len(jtree.get_terminals())!=40 or "OUTGROUP_saff" not in [x.name for x in jtree.get_terminals()]:
        raise ValueError("Japan source nuclear tree tip coverage drift")
    jmap=list(csv.DictReader(io.StringIO(j_map_raw.decode("utf-8-sig"))))
    if len(jmap)!=38 or len({r["paper_japan_member_id"] for r in jmap})!=38:
        raise ValueError("Japan nuclear species concept mapping not full 38")
    name_na={x.name:x.name for x in na.get_terminals() if
        x.name in by and by[x.name]["WCVP_key_collides_source_names"]=="False"}
    name_j={}
    for r in jmap:
        label=r["paper_taxon_concept"]
        tip=r["tip_ids"]
        base="_".join(label.split()[:2])
        if (r["trait_asr_primary_allowed"]!="true"
            or r["n_biological_samples"]!="1"
            or "|" in tip
            or re.search(r"\b(var\.|subsp\.|ssp\.|forma\b)",label)):
            continue
        if base in by and by[base]["WCVP_key_collides_source_names"]=="False":
            if base in name_j.values():
                raise ValueError("Many Japan DNA specimens collapsed to same S1 binomial")
            name_j[tip]=base
    if len(name_na)!=NA_EXPECTED or len(name_j)!=J_EXPECTED:
        raise ValueError(f"Strict full species-only native source match changed {len(name_na)}/{len(name_j)}")
    if set(name_na.values()) & set(name_j.values()):
        raise ValueError("Separate geographic trees must not have duplicate source species")
    return by, {"NA42_genuine_published_nuclear":(na,name_na),
                "Japan12_genuine_176_loci_nuclear":(jtree,name_j)},{
        "NA_tree_Git_blob":git_blob(na_raw),
        "Japan_tree_SHA256":hashlib.sha256(j_raw).hexdigest(),
        "S1_251_CSV_SHA256":hashlib.sha256(original).hexdigest(),
        "Japan_accepted_loci":jmanifest["rootable_loci"]}

def prune(tree,tipmap):
    target=deepcopy(tree)
    for t in list(target.get_terminals()):
        if t.name not in tipmap:
            target.prune(t)
    tips=target.get_terminals()
    if set(x.name for x in tips)!=set(tipmap):raise ValueError("Real tree did not retain exact source tips")
    if any(not np.isfinite(x.branch_length or 0) or (x.branch_length or 0)<0
           for x in target.find_clades() if x is not target.root):
        raise ValueError("Original nuclear phylogeny branch lengths invalid")
    return target

def neighbor_kernel(tree):
    tips=tree.get_terminals()
    n=len(tips); distances=np.full((n,n),np.inf)
    for i in range(n):
        for j in range(i):
            d=float(tree.distance(tips[i],tips[j]))
            distances[i,j]=distances[j,i]=d
    if np.any(~np.isfinite(np.min(distances,axis=1))) or np.any(np.min(distances,axis=1)<=0):
        raise ValueError("True phylogenetic patristic pair distances missing")
    dmin=distances.min(axis=1)
    neighbors=np.isclose(distances,dmin[:,None],rtol=1e-10,atol=1e-12)
    weights=neighbors/neighbors.sum(axis=1,keepdims=True)
    return tips,weights,distances

def observed_nn_stat(x,weights):
    v=np.asarray(x,float)
    if v.ndim!=1 or weights.shape!=(len(v),len(v)) or not np.isfinite(v).all():
        raise ValueError("Exact source taxon angles/annual rainfall need scalar finite array")
    return float(np.sum(weights*np.abs(v[:,None]-v[None,:]))/len(v))

def permutation_test(x,weights,nrep,seed):
    value=observed_nn_stat(x,weights)
    gen=np.random.default_rng(seed)
    null=np.empty(nrep)
    for i in range(nrep):
        null[i]=observed_nn_stat(gen.permutation(x),weights)
    return {"source_nearest_neighbor_mean_absolute_difference":value,
            "source_tip_label_null_mean":float(null.mean()),
            "ratio_observed_to_null":float(value/null.mean()),
            "null_95_interval":[float(v) for v in np.quantile(null,[.025,.975])],
            "one_sided_clustering_p":float((1+np.count_nonzero(null<=value))/(nrep+1)),
            "n_permutations":nrep,"null":"real_nuclear_tree_fixed_tip_label_randomization",
            "geographic_spatial_autocorrelation_not_removed":True,
            "tree_topology_uncertainty_not_propagated":True}

def trait_value(r,k):
    if k=="log_BIO12":
        z=float(r["chelsa_bio12"])
        if z<=0:raise ValueError("Original native annual precipitation invalid")
        return float(np.log(z))
    v=float(r[k])
    if not 0<=v<=180:raise ValueError("S1 original camera angle invalid")
    return v

def analyze_tree(tree,tipmap,by,*,seed,replicates=B):
    report={"n_original_nuclear_tree_tips":len(tree.get_terminals()),
        "source_S1_matched_taxa":len(tipmap),
        "source_S1_matched_original_photos":sum(int(by[k]["original_photo_count"]) for k in tipmap.values()),
        "camera_relative_degree_not_true_downward":True,
        "botanical_literature_U_D_used_for_taxon_selection":False,
        "nuclear_genealogy_only_no_ancestral_state_or_adaptive_gain":True,
        "gates":{}}
    source_tags=[]
    for gate in PHOTO_GATES:
        keep={x:y for x,y in tipmap.items() if int(by[y]["original_photo_count"])>=gate}
        if len(keep)<8:
            report["gates"][f"min_original_photos_{gate}"]={
                "n_source_species":len(keep),"status":"NOT_EVALUABLE_TOO_FEW_NUCLEAR_TIPS",
                "no_formal_permutation_test":True}
            continue
        t=prune(tree,keep)
        tips,w,dist=neighbor_kernel(t)
        fit={"n_source_species":len(tips),
             "n_original_S1_photos":sum(int(by[keep[x.name]]["original_photo_count"]) for x in tips),
             "n_nuclear_nearest_neighbor_tip_ties":int(np.sum(w.max(axis=1)<1)),
             "traits":{}}
        for offset,trait in enumerate(TRAITS):
            x=np.array([trait_value(by[keep[t.name]],trait) for t in tips])
            fit["traits"][trait]=permutation_test(x,w,replicates,seed+1000*gate+offset)
        report["gates"][f"min_original_photos_{gate}"]=fit
        if gate==1:
            for i,t in enumerate(tips):
                neighbor_names=[keep[tips[j].name] for j,z in enumerate(w[i]) if z>0]
                source_tags.append({"DNA_tip":t.name,"S1_original_source_name":keep[t.name],
                    "S1_original_photo_observations":int(by[keep[t.name]]["original_photo_count"]),
                    "head_photo_angle_mean_deg":float(by[keep[t.name]]["mean_image_angle_deg"]),
                    "head_photo_angle_median_deg":float(by[keep[t.name]]["median_image_angle_deg"]),
                    "native_BIO12_mm":float(by[keep[t.name]]["chelsa_bio12"]),
                    "nearest_nuclear_tips_S1_source_names":neighbor_names,
                    "nearest_patristic_distance_same_branch_units_only":float(dist[i].min())})
    return report,source_tags

def bh(vals):
    v=np.asarray(vals,float)
    order=np.argsort(v); z=v[order]*len(v)/np.arange(1,len(v)+1)
    z=np.minimum.accumulate(z[::-1])[::-1]
    result=np.empty_like(z);result[order]=np.minimum(z,1)
    return result

def plot_tree(tree,tipmap,by,output,title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.colors import Normalize
    t=prune(tree,tipmap)
    tips=t.get_terminals(); pos={tip:i for i,tip in enumerate(tips)}
    row={node:None for node in t.find_clades()}
    def ycoord(n):
        if n in pos:return pos[n]
        rows=[ycoord(z) for z in n.clades]
        return float(np.mean(rows))
    def draw_node(ax,node,parent_x):
        x=parent_x+max(0,float(node.branch_length or 0))
        ys=[ycoord(t) for t in node.clades]
        if ys:ax.plot([x,x],[min(ys),max(ys)],color="0.45",lw=0.8)
        if node is not t.root: ax.plot([parent_x,x],[ycoord(node),ycoord(node)],color="0.45",lw=.8)
        for c in node.clades:draw_node(ax,c,x)
    fig=plt.figure(figsize=(12,max(3,0.27*len(tips)+1.8)),layout="constrained")
    gr=fig.add_gridspec(1,4,width_ratios=[5.2,2.2,.6,.6],wspace=.05)
    ax=fig.add_subplot(gr[0,0]);an=fig.add_subplot(gr[0,1],sharey=ax)
    ha=fig.add_subplot(gr[0,2],sharey=ax);hr=fig.add_subplot(gr[0,3],sharey=ax)
    draw_node(ax,t.root,0)
    for i,tip in enumerate(tips):
        label=tipmap[tip.name].replace("Cirsium_","C. ").replace("_"," ")
        an.text(.01,i,label,ha="left",va="center",fontsize=8)
        r=by[tipmap[tip.name]]
        ha.add_patch(Rectangle((0,i-.4),1,.8,facecolor=plt.get_cmap("viridis")(float(r["mean_image_angle_deg"])/180),lw=0))
        rain=np.log(float(r["chelsa_bio12"]))
        rclip=np.clip((rain-np.log(200))/(np.log(5000)-np.log(200)),0,1)
        hr.add_patch(Rectangle((0,i-.4),1,.8,facecolor=plt.get_cmap("cividis")(rclip),lw=0))
    for axis in [ax,an,ha,hr]:
        axis.set_ylim(len(tips)-.2,-.8)
        axis.spines[:].set_visible(False)
        axis.set_yticks([])
    ax.set_xlabel("Published nuclear branch lengths (within panel only)",fontsize=8)
    ax.tick_params(axis="x",labelsize=7)
    an.set_xlim(0,1);an.set_xticks([])
    for axis,label in [(ha,"S1 head\nangle"),(hr,"BIO12")]:
        axis.set_xlim(0,1);axis.set_xticks([]);axis.set_title(label,fontsize=9)
    fig.suptitle(title+"\nTip colours are own S1 IMAGE angles, NOT botanical erect/nodding classes",fontsize=11)
    fig.savefig(output,dpi=170,bbox_inches="tight")
    plt.close(fig)

def run(original251,outdir,*,na_raw,j_raw,j_map_raw,j_manifest_raw,replicates=B):
    by,panels,source=read_sources(original251,na_raw,j_raw,j_map_raw,j_manifest_raw)
    results={}
    for i,(label,(tree,tipmap)) in enumerate(panels.items()):
        table,ledger=analyze_tree(tree,tipmap,by,seed=SEED+i*100000,replicates=replicates)
        results[label]=table
        if outdir is not None:
            plot_tree(tree,tipmap,by,Path(outdir)/f"T2_S1_CONTINUOUS_{label}_HEAD_BIO12.png",label.replace("_"," "))
            with (Path(outdir)/f"T2_S1_{label}_SOURCE_TIP_LEDGER.csv").open("w",newline="") as fd:
                writer=csv.DictWriter(fd,fieldnames=list(ledger[0]))
                writer.writeheader()
                for v in ledger:
                    v={**v,"nearest_nuclear_tips_S1_source_names":"|".join(v["nearest_nuclear_tips_S1_source_names"])}
                    writer.writerow(v)
    primary=[]
    for label,stats in results.items():
        primary.append((label,"mean_image_angle_deg",
                        stats["gates"]["min_original_photos_1"]["traits"]["mean_image_angle_deg"]))
        primary.append((label,"log_BIO12",
                        stats["gates"]["min_original_photos_1"]["traits"]["log_BIO12"]))
    corrected=bh([entry[2]["one_sided_clustering_p"] for entry in primary])
    for (lab,endpoint,d),q in zip(primary,corrected):
        d["BH_q_across_4_primary_nuclear_panel_trait_tests"]=float(q)
    return {"schema":"azami_V4_S1_CONTINUOUS_GENUINE_NA42_JAPAN12_MULTILOCUS_NUCLEAR_CLUSTERING_v1",
        "status":"TWO_AUTHENTIC_SEPARATE_PUBLISHED_NUCLEAR_PHYLOGENIES_S1_PHOTO_ANGLE_TESTED",
        "original_S1_source_taxon_count":251,
        "original_S1_eligible_photo_observations":170066,
        "genuinely_independent_nuclear_panels_separate_not_grafted":2,
        "authentic_multi_locus_nuclear_matched_source_species_total":54,
        "chloroplast_tree_not_used":True,
        "ITS_multicopy_rDNA_not_primary_nuclear_species_tree":True,
        "Q20_posthoc_failed_primary_QC_NOT_USED":True,
        "no_botanical_erect_nodding_states_used":True,
        "primary_trait_S1_continuous_camera_angle_not_verified_gravity":True,
        "climate_niche_not_ancestral_environment":True,
        "tree_source_provenance":source,
        "fixed_seed":SEED,"label_permutations_each_nuclear_test":replicates,
        "primary_two_geographies":{"North_America":"NA42_genuine_published_nuclear",
                                   "Japan":"Japan12_genuine_176_loci_nuclear"},
        "nuclear_panels":results}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original251",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    p.add_argument("--NA-tree",type=Path)
    p.add_argument("--Japan-tree",type=Path)
    p.add_argument("--Japan-map",type=Path)
    p.add_argument("--Japan-manifest",type=Path)
    p.add_argument("--permutations",type=int,default=B)
    a=p.parse_args()
    if a.permutations<999 or a.permutations>20000:raise ValueError("Insufficient or excessive diagnostic label permutations")
    local=(a.NA_tree,a.Japan_tree,a.Japan_map,a.Japan_manifest)
    if any(x is not None for x in local) and not all(x is not None for x in local):
        raise ValueError("Must provide all pinned nuclear source files, not a mixed-source tree")
    if all(x is not None for x in local):
        raws=[x.read_bytes() for x in local]
    else:
        raws=[upstream(NA_REPO,NA_REF,NA_TREE),upstream(J_REPO,J_REF,J_TREE),
              upstream(J_REPO,J_REF,J_MAP),upstream(J_REPO,J_REF,J_MANIFEST)]
    a.out_dir.mkdir(parents=True,exist_ok=True)
    result=run(a.original251,a.out_dir,
               na_raw=raws[0],j_raw=raws[1],j_map_raw=raws[2],j_manifest_raw=raws[3],
               replicates=a.permutations)
    (a.out_dir/"T2_S1_CONTINUOUS_GENUINE_TWO_NUCLEAR_PHYLO_CLUSTERING.json").write_text(
        json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    for label,v in result["nuclear_panels"].items():
        primary=v["gates"]["min_original_photos_1"]["traits"]
        print("GENUINE_NUCLEAR_S1_ANGLE_CLUSTER="+json.dumps({
          "tree":label,"n":v["source_S1_matched_taxa"],
          "S1_angle_nn_p":primary["mean_image_angle_deg"]["one_sided_clustering_p"],
          "logBIO12_nn_p":primary["log_BIO12"]["one_sided_clustering_p"],
          "photos":v["source_S1_matched_original_photos"]},sort_keys=True))

if __name__=="__main__":
    main()
