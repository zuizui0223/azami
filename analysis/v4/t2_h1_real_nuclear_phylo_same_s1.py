"""Matched-original-S1 rainfall versus botanical U/D under two INDEPENDENT nuclear trees.

North America: exact 44 S1 species on author Siniscalchi 2023 NA63 nuclear tree.
Japan: 12 exact S1 species on independently reconstructed EAzami Japan38
176-locus nuclear tree (names/vouchers not yet fully reconciled).
Each tree is fitted separately under a symmetric Mk neutral null conditional
on its observed number of nodding taxa; H1 is the one-sided positive contrast.
No global 66-tip species tree, clock, historical rain, causal selection, or
validated independent gain count is asserted.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import csv
import hashlib
import io
from itertools import combinations, product
import json
from pathlib import Path
import urllib.request

from Bio import Phylo
import numpy as np
import pandas as pd

from analysis.v4.t2_phylogenetic_climate_null import conditional_mk_null

NA_REPO = "https://raw.githubusercontent.com/carol-siniscalchi/cirsium/c925d4af86160064448f794a9715f70b917a2786/traits/"
EA_REPO = "https://raw.githubusercontent.com/zuizui0223/EAzami/c9870e91604ed1514344ddc8f765eae5237dbf27/data/evidence/"
NA_TREE = "cirsium.phyluce95p.ultrametric.tre"
EA_TREE = "japan38_comp1061_primary_tree_v1.nwk"
EA_MAP = "japan38_comp1061_concept_map_v1.csv"
NA_GIT_BLOB = "6417a38793427bc7c8562dc0e2cab59b74d8c316"
EA_TREE_SHA256 = "ba592d3a9f52b5481314ab20e3202048505cdf72a3fc15134d5dca261b9d9b7e"
EA_MAP_GIT_BLOB = "cc093a8e0269b5dfe8086f938c3f29c9e02a9866"
SEED = 20261023


def git_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def download(url, digest, *, git=True):
    with urllib.request.urlopen(url, timeout=90) as response:
        payload=response.read()
    observed=git_blob(payload) if git else hashlib.sha256(payload).hexdigest()
    if observed!=digest:
        raise ValueError("External published sequence tree or concept-map SHA drift")
    return payload


def source_panels(csv_path):
    source=pd.read_csv(csv_path)
    if len(source)!=251 or len(set(source["accepted_key"]))!=244 or sum(source["original_photo_count"])!=170066:
        raise ValueError("Original 251 source names, 244 WCVP keys or 170066 photographs changed")
    na=source[source["independent_botanical_authority_lane"]=="NA63_author_scored"].copy()
    na["phylogenetic_tip"]=na["original_source_name"].str.replace(" ","_",regex=False)
    if len(na)!=44 or sum(na.botanical_U_D_state=="D")!=3:
        raise ValueError("NA44 authentic H1 original photo/author state sample changed")
    original_na=download(NA_REPO+NA_TREE,NA_GIT_BLOB)
    tree_na=Phylo.read(io.StringIO(original_na.decode("utf-8-sig")),"newick")
    if len(tree_na.get_terminals())!=63:
        raise ValueError("Published NA nuclear tree has changed taxon count")

    map_raw=download(EA_REPO+EA_MAP,EA_MAP_GIT_BLOB)
    cmap=list(csv.DictReader(io.StringIO(map_raw.decode("utf-8-sig"))))
    if len(cmap)!=38:
        raise ValueError("Published Japan38 concept map changed")
    exact={}
    for r in cmap:
        parts=r["paper_taxon_concept"].split()
        if (r["trait_asr_primary_allowed"]!="true"
                or r["n_biological_samples"]!="1"
                or len(r["tip_ids"].split("|"))!=1
                or any(q in parts for q in ("var.","subsp.","ssp.","f."))
                or len(parts)<2 or parts[0]!="Cirsium"):
            continue
        key=" ".join(parts[:2])
        if key in exact:
            raise ValueError("Duplicate Japan species concept")
        exact[key]=r["tip_ids"]
    east=source[source["independent_botanical_authority_lane"]=="East_Asia_botanical_authority"].copy()
    jp=east[east.original_source_name.isin(exact)].copy()
    jp["phylogenetic_tip"]=jp["original_source_name"].map(exact)
    if len(jp)!=12 or sum(jp.botanical_U_D_state=="D")!=5 or jp.phylogenetic_tip.nunique()!=12:
        raise ValueError("JP12 exact original S1/Molecular core not supported")
    original_jp=download(EA_REPO+EA_TREE,EA_TREE_SHA256,git=False)
    tree_jp=Phylo.read(io.StringIO(original_jp.decode("utf-8-sig")),"newick")
    if len(tree_jp.get_terminals())!=40 or "OUTGROUP_saff" not in {z.name for z in tree_jp.get_terminals()}:
        raise ValueError("Original Japan38 nuclear tree or outgroup missing")
    tree_jp.root_with_outgroup("OUTGROUP_saff")
    return [
        make_panel("NA44",na,tree_na),
        make_panel("Japan12",jp,tree_jp),
    ]


def make_panel(name,df,full_tree):
    tree=deepcopy(full_tree)
    accepted=set(df.phylogenetic_tip)
    if not accepted<={z.name for z in tree.get_terminals()}:
        raise ValueError("Source molecular tree lacks exact original S1 names")
    for tip in list(tree.get_terminals()):
        if tip.name not in accepted:
            tree.prune(tip)
    tips=[t.name for t in tree.get_terminals()]
    if len(tips)!=len(df):
        raise ValueError("Distinct original tree tips lost")
    data=df.set_index("phylogenetic_tip").loc[tips]
    y=(data.botanical_U_D_state=="D").astype(int).to_numpy()
    rain=data.chelsa_bio12.to_numpy(dtype=float)
    temperature=data.chelsa_bio01.to_numpy(dtype=float)
    rain_z=(rain-rain.mean())/rain.std(ddof=0)
    temp_z=(temperature-temperature.mean())/temperature.std(ddof=0)
    rem=rain_z-temp_z*(rain_z@temp_z)/(temp_z@temp_z)
    rem/=rem.std(ddof=0)
    X=np.column_stack((rain_z,rem))
    n=len(y);k=int(y.sum())
    # Allow explicitly labeled small-n alternative panels with >=3 in each state;
    # the registered NA44/Japan12 principal comparison remains unchanged.
    if k<3 or n-k<3:
        raise ValueError("Source state composition too thin")
    return {"name":name,"tree":tree,"tips":tips,"data":data,
            "y":y,"X":X,"n":n,"k":k,"weight":2*k*(n-k)/n,
            "observed":(X[y==1].mean(axis=0)-X[y==0].mean(axis=0)),
            "rain":rain,"temp":temperature}


def neutral_label_null(p):
    n,k=p["n"],p["k"]
    picks=np.asarray(list(combinations(range(n),k)),dtype=np.int16)
    XX=p["X"]
    # Direct exact k-label assignments, not draws from a non-phylogenetic model.
    return (n*XX[picks].mean(axis=1)-n*XX.mean(axis=0))/(n-k)


def phylo_null(p,mode,seed,draws):
    tree=p["tree"]
    observed={tip:int(state) for tip,state in zip(p["tips"],p["y"])}
    null,meta=conditional_mk_null(
        tree,observed,p["X"],root_mode=mode,seed=seed,
        target_draws=draws,max_proposals=1000000,min_accepted=draws)
    if len(null)!=draws:
        raise ValueError("Not enough exact D-prevalence matched real phylogenetic null draws")
    return null,meta


def summarize(data_path,draws=30000):
    pn,pj=source_panels(data_path)
    panels={}
    simulation={}
    naive={}
    for i,p in enumerate((pn,pj)):
        name=p["name"]
        obs=p["observed"]
        nc=neutral_label_null(p)
        naive[name]=nc
        record={
            "nuclear_tips_matched":p["n"],"D":p["k"],"U":p["n"]-p["k"],
            "source_name_to_tree_tip":dict(zip(p["data"].original_source_name,p["tips"])),
            "mean_BIO12_down_mm":float(p["rain"][p["y"]==1].mean()),
            "mean_BIO12_erect_mm":float(p["rain"][p["y"]==0].mean()),
            "rain_mm_down_minus_erect":float(p["rain"][p["y"]==1].mean()-p["rain"][p["y"]==0].mean()),
            "observed_standardized_BIO12_then_residual_BIO1":obs.tolist(),
            "naive_exact_permutations":len(nc),
            "naive_H1_one_sided_p":np.mean(nc>=obs-1e-12,axis=0).tolist(),
            "phylogeny_symmetric_Mk_conditional":{},
        }
        for j,root in enumerate(("erect","stationary")):
            null,details=phylo_null(p,root,SEED+i*200+j*500,draws)
            simulation[(name,root)]=null
            record["phylogeny_symmetric_Mk_conditional"][root]={
                "n_conditioned_states":len(null),
                "fitted_Mk_rate":details["fitted_Mk_rate_tree_units"],
                "null_mean":null.mean(axis=0).tolist(),
                "H1_one_sided_p":((np.sum(null>=obs-1e-12,axis=0)+1)/(draws+1)).tolist(),
            }
        panels[name]=record

    weight_sum=pn["weight"]+pj["weight"]
    observed=(pn["weight"]*pn["observed"]+pj["weight"]*pj["observed"])/weight_sum
    n1,n2=naive["NA44"],naive["Japan12"]
    exact_n=len(n1)*len(n2)
    exact_exceed=np.zeros(2,dtype=np.int64)
    for i in range(0,len(n1),512):
        block=(pn["weight"]*n1[i:i+512,None,:]+pj["weight"]*n2[None,:,:])/weight_sum
        exact_exceed+=(block>=observed-1e-12).sum(axis=(0,1))
    calibrated={}
    for root1,root2 in product(("erect","stationary"),repeat=2):
        n=(pn["weight"]*simulation[("NA44",root1)]+
           pj["weight"]*simulation[("Japan12",root2)])/weight_sum
        calibrated[root1+"__"+root2]={
            "n_joint_conditional_null":draws,
            "H1_one_sided_p":((np.sum(n>=observed-1e-12,axis=0)+1)/(draws+1)).tolist(),
            "null_mean":n.mean(axis=0).tolist()}
    return {
        "schema":"azami_T2_H1_REAL_NUCLEAR_NA44_Japan12_climate_permutation_v1",
        "source_exact_S1_photo_taxa":251,"n_exact_nuclear_S1_joint":56,
        "number_of_all_66_author_scored_unmatched_to_these_two_trees":10,
        "environment_variables":["standardized_BIO12","standardized_BIO12_linear_residual_on_BIO1"],
        "published_source_git_blob_NA":NA_GIT_BLOB,
        "EAzami_tree_sha256":EA_TREE_SHA256,
        "EAzami_map_git_blob":EA_MAP_GIT_BLOB,
        "photo_input_sha256":hashlib.sha256(Path(data_path).read_bytes()).hexdigest(),
        "regions":panels,
        "separately_phylogeny_calibrated_combination":{
            "weights":[pn["weight"],pj["weight"]],
            "estimated_effect":observed.tolist(),
            "nonphylo_exact_joint_permutations":exact_n,
            "nonphylo_one_sided_p":(exact_exceed/exact_n).tolist(),
            "phylo_Mk_H1_root_sensitivity":calibrated,
        },
        "limitations":{
            "not_full66_global_phylogeny":True,
            "Japanese_exact_names_not_same_biological_voucher_verified":True,
            "unmatched_10_taxa_not_set_to_erect":True,
            "NA_time_tree_and_Japan_substitution_tree_fit_Mk_separately":True,
            "phylogenetic_null_not_rain_dependent_evolution":True,
            "no_geography_adjustment":True,
            "root_and_rates_fixed_per_scenario":True,
            "no_independent_origin_number_or_ancestral_precipitation_inferred":True,
            "posthoc_directional_H1_sensitivity_not_confirmation":True,
        },
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original251",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--draws",type=int,default=30000)
    a=parser.parse_args()
    if not 2000<=a.draws<=100000:
        raise ValueError("Conditional Mk reproducibility draws out of range")
    result=summarize(a.original251,draws=a.draws)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_SAME_56_S1_PHYLOGENY_RAIN_H1="+json.dumps({
        "n":result["n_exact_nuclear_S1_joint"],
        "observed":result["separately_phylogeny_calibrated_combination"]["estimated_effect"],
        "naive":result["separately_phylogeny_calibrated_combination"]["nonphylo_one_sided_p"],
        "phylo":{k:v["H1_one_sided_p"] for k,v in
          result["separately_phylogeny_calibrated_combination"]["phylo_Mk_H1_root_sensitivity"].items()},
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
