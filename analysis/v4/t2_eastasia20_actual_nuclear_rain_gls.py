"""East Asian genuine 20-tip nuclear tree × botanical U/D × GBIF-CHELSA rainfall.

NOT photo-derived: source head postures independently described in EAzami
botanical authorities, and climate comes from a separate GBIF occurrence
CHELSA environmental source atlas, NOT image locations/count selection.

Of 20 nuclear ingroup tips, 17 have original U/D botanical source coding,
3 ambiguous states remain UNKNOWN. The GBIF climatological data cover
14 unequivocally classified exact biological concept-name candidates.
True EAzami nuclear ML tree has substitution/site branches and a safflower
outgroup; it is not a dated tree or complete genus/global tree.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from collections import Counter
import urllib.request

import numpy as np
from Bio import Phylo

from analysis.v4.t2_pgls_placement_screen import git_blob_sha, parsimony
from analysis.v4.t2_na63_actual_tree_rain_gls import (
    brownian_shared_ancestry, fit_models, internal_marginals,
    equal_variance_independent_OLS,
)

EA_REPO="zuizui0223/EAzami"
EA_COMMIT="c9870e91604ed1514344ddc8f765eae5237dbf27"
EA_ROOT=f"https://raw.githubusercontent.com/{EA_REPO}/{EA_COMMIT}/data/evidence/"
TREE_FILE="full20_comp1061_primary_tree_v1.nwk"
TRAIT_FILE="orientation_comp1061_20tip_source_crosswalk_v1.csv"
TREE_BLOB="1e3e6ffaf128ccd439287070b8c3c81385f5df34"
TRAIT_BLOB="f655544e77f1659323c1e915567dacbb51b01cc3"
LOCAL_CLIMATE=Path(__file__).resolve().parent/"results"/"20261022_t2_public_niche_phylogeny"/"source_pinned_eazami_niche_centroids.csv"
CLIMATE_BLOB="2216152d185b450958ebf59c9de444dd3f81f0b3"
SOURCE_JP_TW_CLIMATE_N=14
ORIGINAL_NUCLEAR_TIPS=20
SOURCE_OBSERVED=17
SOURCE_UNKNOWN=3
SOURCE_DOWN=8
TREE_OUTGROUP="OUTGROUP_saff"


def fetch_fixed(path, blob):
    url=EA_ROOT+path
    with urllib.request.urlopen(url,timeout=90) as handle:
        payload=handle.read()
    if git_blob_sha(payload)!=blob:
        raise ValueError("EAzami original published nuclear tree/flower source SHA drift")
    return payload


def rooted_ingroup(tree,outgroup):
    """Source nuclear root with correct removal of outgroup-only long stem."""
    if outgroup not in {t.name for t in tree.get_terminals()}:
        raise ValueError("Published original nuclear outgroup was lost")
    tree.root_with_outgroup(outgroup)
    tree.prune(outgroup)
    # Biopython joins the long outgroup edge to the ingroup root edge.
    # A removed outgroup is NOT an ancestral ingroup branch and must not
    # inflate Brownian shared-ancestry covariance across ingroup tips.
    tree.root.branch_length=0.0
    return tree


def read_sources(*,fetch=fetch_fixed,climate=LOCAL_CLIMATE):
    tree_bytes=fetch(TREE_FILE,TREE_BLOB)
    trait_bytes=fetch(TRAIT_FILE,TRAIT_BLOB)
    tree=Phylo.read(io.StringIO(tree_bytes.decode("utf-8-sig")),"newick")
    if set(t.name for t in tree.get_terminals()) is None:
        raise ValueError("Unrecognized original nuclear tree")
    if len(tree.get_terminals())!=21 or TREE_OUTGROUP not in [t.name for t in tree.get_terminals()]:
        raise ValueError("Required original nuclear outgroup and 20 source tips missing")
    original=list(csv.DictReader(io.StringIO(trait_bytes.decode("utf-8-sig"))))
    if len(original)!=ORIGINAL_NUCLEAR_TIPS or len({r["tip_id"] for r in original})!=ORIGINAL_NUCLEAR_TIPS:
        raise ValueError("Original source botanical tip list changed")
    if not {r["tip_id"] for r in original} <= {t.name for t in tree.get_terminals()}:
        raise ValueError("Molecular tree to independent flora tip identity changed")
    if any(r["analysis_state"] not in ("U","D","?") or
           r["same_primary_voucher_state_verified"]!="false" for r in original):
        raise ValueError("Unknown botanical head phenotype incorrectly promoted to exact specimen")
    coded={r["tip_id"]:{"U":0,"D":1}[r["analysis_state"]] for r in original if r["analysis_state"]!="?"}
    count=Counter(r["analysis_state"] for r in original)
    if count!={"U":9,"D":SOURCE_DOWN,"?":SOURCE_UNKNOWN} or len(coded)!=SOURCE_OBSERVED:
        raise ValueError("True East Asia 20 published morphology changed")
    tree=rooted_ingroup(tree,TREE_OUTGROUP)
    if set(coded)-{z.name for z in tree.get_terminals()}:
        raise ValueError("Outgroup pruning lost original botanical source tip")
    climate_bytes=Path(climate).read_bytes()
    if git_blob_sha(climate_bytes)!=CLIMATE_BLOB:
        raise ValueError("Public original GBIF-CHELSA taxon centroids changed")
    rows=list(csv.DictReader(io.StringIO(climate_bytes.decode("utf-8-sig"))))
    if len(rows)!=SOURCE_JP_TW_CLIMATE_N:
        raise ValueError("No verified 14-taxon GBIF climatological source")
    name_from_row=lambda z:z.replace(" var. ","_var_").replace(" ","_")
    bytaxon={}
    for r in rows:
        tip=name_from_row(r["taxon"])
        if tip in bytaxon or tip not in coded or r["region_source"] not in ("japan","taiwan"):
            raise ValueError("Species not exact binary source morph / nuclear taxon or duplicated GBIF")
        if (int(r["n"])<1 or float(r["BIO12_mean"])<=0 or
            not np.isfinite(float(r["BIO1_mean"]))):
            raise ValueError("Source GBIF taxon climate invalid")
        bytaxon[tip]={"BIO12":float(r["BIO12_mean"]),
                      "BIO1":float(r["BIO1_mean"]),
                      "BIO15":float(r["BIO15_mean"]),
                      "region":r["region_source"],
                      "occurrence_n":int(r["n"]),
                      "taxon":r["taxon"]}
    if len(bytaxon)!=SOURCE_JP_TW_CLIMATE_N or Counter(coded[t] for t in bytaxon)!={0:8,1:6}:
        raise ValueError("Original East Asia nuclear-authority GBIF climate matches changed")
    return tree,coded,bytaxon,{
        "original_nuclear_tree_git_blob":TREE_BLOB,
        "original_morphology_git_blob":TRAIT_BLOB,
        "GBIF_CHELSA_taxon_source_git_blob":CLIMATE_BLOB,
        "actual_nuclear_tree_sha256":hashlib.sha256(tree_bytes).hexdigest(),
        "morphology_source_sha256":hashlib.sha256(trait_bytes).hexdigest(),
        "GBIF_CHELSA_source_sha256":hashlib.sha256(climate_bytes).hexdigest(),
    }


def fitted_tip_subtree(tree, tips):
    from copy import deepcopy
    t=deepcopy(tree)
    for clade in list(t.get_terminals()):
        if clade.name not in tips:
            t.prune(clade)
    if len(t.get_terminals())!=len(tips):
        raise ValueError("Authentic source nuclear subtree lost a real tip")
    return t


def analyze(tree,coded,climate,checks):
    accepted={tip.name for tip in tree.get_terminals()}
    if len(accepted)!=20 or len(coded)!=17 or len(climate)!=14:
        raise ValueError("Original genus nuclear tree/plant trait/GBIF data changed")
    if not set(coded)<=accepted or not set(climate)<=set(coded):
        raise ValueError("Botanical states or GBIF environmental niche mismatched to nuclear taxa")
    pruned=fitted_tip_subtree(tree,set(climate))
    tips=[x.name for x in pruned.get_terminals()]
    D=np.array([coded[t] for t in tips],float)
    rain=np.array([climate[t]["BIO12"] for t in tips],float)
    temp=np.array([climate[t]["BIO1"] for t in tips],float)
    region=np.array([climate[t]["region"]=="taiwan" for t in tips],float)
    if len(D)!=14 or D.sum()!=6 or len(set(region))!=2:
        raise ValueError("No real comparison of six nodding eight erect two areas")
    true_cov=brownian_shared_ancestry(pruned,tips)
    ztemp=(temp-temp.mean())/temp.std(ddof=0)
    X=np.column_stack((np.ones(len(D)),D))
    X_temp=np.column_stack((X,ztemp))
    X_full=np.column_stack((X,ztemp,region))
    models={
        "unadjusted_BIO12_mm":fit_models(rain,X,true_cov),
        "BIO1_adjusted_BIO12_mm":fit_models(rain,X_temp,true_cov),
        "BIO1_and_Taiwan_indicated_BIO12_mm":fit_models(rain,X_full,true_cov),
    }
    unrelated_IID={
        "unadjusted_BIO12_mm":equal_variance_independent_OLS(rain,X),
        "BIO1_adjusted_BIO12_mm":equal_variance_independent_OLS(rain,X_temp),
        "BIO1_and_Taiwan_indicated_BIO12_mm":equal_variance_independent_OLS(rain,X_full),
    }
    topology=fitted_tip_subtree(tree,set(coded))
    # Three unresolved phenotypes NOT included in the 17-tip parsimony;
    # never assign an inherited posture to them.
    mpo=parsimony(topology,coded)
    ancestor={
        r:internal_marginals(topology,coded,root_prior=r)
        for r in ("erect","stationary")}
    return {
        "schema":"azami_v4_T2_eastasia20_true_nuclear_botanical_GBIF_CHELSA_rain_GLs_20261023_v1",
        "status":"EAST_ASIA_ORIGINAL20_TREE_17_BOTANICAL_STATES_14_GBIF_NICHES_REAL_PHYLOGENY",
        "source":checks,
        "source_original_tree_total_ingroup":20,
        "known_independent_botanical_posture_tips":17,
        "unknown_or_developmental_traits_not_imputed":3,
        "source_GBIF_CHELSA_climate_matched_nuclear_tips":14,
        "GBIF_climate_coded_D":6,"GBIF_climate_coded_U":8,
        "GBIF_mean_annual_rain_down_mm":float(rain[D==1].mean()),
        "GBIF_mean_annual_rain_upright_mm":float(rain[D==0].mean()),
        "GBIF_raw_D_minus_U_rain_mm":float(rain[D==1].mean()-rain[D==0].mean()),
        "GBIF_n_occurrences_by_species":{t:climate[t]["occurrence_n"] for t in tips},
        "tree_is_same_original_20_independently_assembled_Comp1061":True,
        "root_outgroup_original_safflower_used_before_pruning":True,
        "tree_branch_length_unit":"substitutions/site; NOT calibrated absolute time",
        "GLS_on_same14_tip_covariate_panels":models,
        "exact_same_14tip_equal_variance_independent_OLS":unrelated_IID,
        "pagel_lambda_zero_is_weighted_diagonal_Brownian_not_IID_for_unscaled_EA_tree":True,
        "author_botanical_17_tip_minimal_total_changes_not_acquisition_count":mpo,
        "conditional_two_state_Mk_ancestral_state_marginals_17_tips":ancestor,
        "true_external_species_concept_to_same_voucher_phenotype_confirmed":False,
        "phenotype_scores_source_morphology_not_image":True,
        "climate_source_independent_GBIF_CHELSA_not_photo_occurrence":True,
        "reconstructed_ancestral_BIO12_at_U_D_origins":False,
        "whole_Cirsium_global_Moreyra_350_tree_obtained":False,
        "genuine_21tip_nuclear_root_alone_proves_adaptive_transition":False,
        "trait_rain_correlations_not_evidence_of_historical_rain_selection":True,
    }


def render_phylogram(tree,coded,climate,path):
    """Actual branch-proportional topology; separate columns for names and rain."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    tips=tree.get_terminals()
    coords=tree.depths()
    positions={t:float(i) for i,t in enumerate(tips)}
    for node in tree.find_clades(order="postorder"):
        if not node.is_terminal():
            positions[node]=sum(positions[c] for c in node.clades)/len(node.clades)
    longest=max(float(coords[t]) for t in tips)
    if not longest>0:
        raise ValueError("Original EAzami nuclear branch lengths not positive")
    scaled_x=lambda node: 5.0*float(coords[node])/longest
    label_col=5.30
    rain_start=10.65
    rain_width=4.2
    maxrain=max(rec["BIO12"] for rec in climate.values())
    fig,ax=plt.subplots(figsize=(17,12),dpi=180)
    for node in tree.find_clades(order="preorder"):
        x=scaled_x(node)
        if not node.is_terminal():
            ys=[positions[c] for c in node.clades]
            ax.plot([x,x],[min(ys),max(ys)],c="#87959b",lw=.85,zorder=1)
        for child in node.clades:
            ax.plot([x,scaled_x(child)],[positions[child],positions[child]],
                    c="#87959b",lw=.85,zorder=1)
    for t in tips:
        name=t.name.replace("Cirsium_","C. ").replace("_"," ")
        state=coded.get(t.name,None)
        color="#b8583e" if state==1 else "#577688" if state==0 else "#999999"
        y=positions[t]
        ax.scatter([scaled_x(t)],[y],s=28,c=color,zorder=5)
        ax.text(label_col,y,name,fontsize=9,
                fontweight="bold" if state==1 else "normal",
                color=color,va="center",clip_on=False)
        if t.name in climate:
            mm=climate[t.name]["BIO12"]
            ax.barh(y,rain_width*mm/maxrain,left=rain_start,height=.65,
                    color=color,alpha=.72,zorder=3)
            ax.text(rain_start+rain_width+.23,y,f"{mm:,.0f}",
                    fontsize=8,va="center")
        else:
            ax.text(rain_start+.15,y,"no GBIF niche",fontsize=8,
                    va="center",color="#999999")
    ax.axvline(rain_start-.18,c="#dddddd",lw=0.8)
    ax.set_xlim(-.3,rain_start+rain_width+1.6)
    ax.set_ylim(len(tips)+.35,-1.7)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)
    ax.text(0,-1.2,"Genuine East Asian 20-tip nuclear phylogeny",
            fontsize=16,fontweight="bold")
    ax.text(0,-.38,"Author botanical capitulum posture • original Comp1061 nuclear branches",
            fontsize=9)
    ax.text(rain_start,-.35,"GBIF-CHELSA BIO12 (mm/year)",
            fontsize=9,fontweight="bold")
    ax.legend(handles=[
        Line2D([0],[0],marker="o",color="#b8583e",linestyle="none",label="Nodding (D): 8"),
        Line2D([0],[0],marker="o",color="#577688",linestyle="none",label="Erect (U): 9"),
        Line2D([0],[0],marker="o",color="#999999",linestyle="none",label="Unknown: 3"),
    ],loc="lower left",frameon=False,fontsize=9,
      bbox_to_anchor=(0.01,.03))
    fig.text(.10,.02,"BIO12: modern taxon range niche, not ancestral rain. No photographs or EXIF classification.",
             fontsize=9)
    fig.savefig(path,bbox_inches="tight",dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir",required=True,type=Path)
    parser.add_argument("--no-figure",action="store_true")
    a=parser.parse_args()
    tree,coded,climate,manifest=read_sources()
    report=analyze(tree,coded,climate,manifest)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"T2_EASTASIA20_REAL_NUCLEAR_BOTANICAL_GBIF_RAIN_COMPARISON.json").write_text(
        json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    with (a.out_dir/"T2_EA20_ORIGINAL_NUCLEAR_BOTANICAL_TIP_GBIF_BIO12.csv").open("w",newline="",encoding="utf-8") as stream:
        w=csv.DictWriter(stream,fieldnames=[
            "source_original_nuclear_tip","original_independent_botanical_state",
            "GBIF_CHELSA_source_tips","BIO12_mean_mm","BIO1_mean_source","region","occurrence_count",
            "same_botanical_voucher_confirmed"])
        w.writeheader()
        for t in tree.get_terminals():
            orig=coded.get(t.name,None)
            geo=climate.get(t.name)
            w.writerow({
                "source_original_nuclear_tip":t.name,
                "original_independent_botanical_state":("D" if orig==1 else "U" if orig==0 else "?"),
                "GBIF_CHELSA_source_tips":bool(geo),
                "BIO12_mean_mm":geo["BIO12"] if geo else "",
                "BIO1_mean_source":geo["BIO1"] if geo else "",
                "region":geo["region"] if geo else "",
                "occurrence_count":geo["occurrence_n"] if geo else "",
                "same_botanical_voucher_confirmed":False})
    if not a.no_figure:
        render_phylogram(tree,coded,climate,
            a.out_dir/"T2_EASTASIA20_REAL_NUCLEAR_TREE_BOTANICAL_POSTURE_GBIF_RAIN.png")
    print("T2_EA20_ACTUAL_TREE_NOT_PHOTOS="+json.dumps({
        "n":20,"binary":17,"GBIF":14,
        "rain_D_minus_U":report["GBIF_raw_D_minus_U_rain_mm"],
        "gls_BM":report["GLS_on_same14_tip_covariate_panels"]["BIO1_and_Taiwan_indicated_BIO12_mm"]["Brownian_complete_phylogeny"],
        "gls_lambda":report["GLS_on_same14_tip_covariate_panels"]["BIO1_and_Taiwan_indicated_BIO12_mm"]["estimated_lambda"],
        "parsimony":report["author_botanical_17_tip_minimal_total_changes_not_acquisition_count"]},
        sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
