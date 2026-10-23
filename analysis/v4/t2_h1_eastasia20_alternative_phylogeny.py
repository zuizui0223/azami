"""Alternative East Asia 20-tip genuine nuclear topology screening on exact 8 S1 photo taxa.

Six AU non-rejected topologies belong to ONE source-data study; not independent
convergent evolution events. The EA20 eight original S1 taxa overlap Japan38
two-taxon names; do not append eight to Japan12 in the same combined model.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

from Bio import Phylo
import pandas as pd

from analysis.v4.t2_h1_real_nuclear_phylo_same_s1 import (
    make_panel, neutral_label_null, phylo_null,
)

SOURCE_EAZAMI_COMMIT = "c9870e91604ed1514344ddc8f765eae5237dbf27"
SOURCE_ARTIFACT_SHA256 = "04bad63dafae6bf3f91e0b44be8cd5977eb30f14ca2294cd783a3b0f660c9b74"
ORIGINAL_SOURCE_CROSSWALK_GIT_BLOB = "f655544e77f1659323c1e915567dacbb51b01cc3"
N_TOPOLOGIES=6
N_NUCLEAR_TIPS=20
N_MATCHED_S1=8
D_COUNT=5
U_COUNT=3
SOURCE_NAME_REQUIRE={
    "Cirsium tatakaense","Cirsium morii","Cirsium yezoense",
    "Cirsium suffultum","Cirsium kamtschaticum","Cirsium kawakamii",
    "Cirsium brevicaule","Cirsium irumtiense",
}
EXACT_TREE_SENSITIVITY_FOCUS=True


def git_blob(content):
    return hashlib.sha1(b"blob "+str(len(content)).encode()+b"\0"+content).hexdigest()


def assess(s1,author_source,trees,draws=30000):
    original=pd.read_csv(s1)
    if len(original)!=251 or sum(original.original_photo_count)!=170066:
        raise ValueError("Original S1 photographs/named taxon count drift")
    body=Path(author_source).read_bytes()
    if git_blob(body)!=ORIGINAL_SOURCE_CROSSWALK_GIT_BLOB:
        raise ValueError("Pinned EAzami independent morphology source changed")
    crosswalk=list(csv.DictReader(io.StringIO(body.decode("utf-8-sig"))))
    if len(crosswalk)!=20:
        raise ValueError("Original 20 EAzami phylogenetic tips missing")
    codes={row["accepted_taxon"]:row for row in crosswalk}
    matched=set(original.original_source_name)&set(codes)
    if matched!=SOURCE_NAME_REQUIRE:
        raise ValueError("Original eight molecular-photo source names drift")
    frame=original[original.original_source_name.isin(SOURCE_NAME_REQUIRE)].copy()
    for row in frame.itertuples():
        assert codes[row.original_source_name]["analysis_state"]==row.botanical_U_D_state
    frame["phylogenetic_tip"]=frame.original_source_name.str.replace(" ","_",regex=False)
    if len(frame)!=N_MATCHED_S1 or sum(frame.botanical_U_D_state=="D")!=D_COUNT:
        raise ValueError("Eight S1 photo taxa (D5 U3) were not retained")
    dir_path=Path(trees)
    paths=sorted(dir_path.glob("*.nwk"))
    if len(paths)!=N_TOPOLOGIES:
        raise ValueError("Six separately archived AU nonrejected source topology hypotheses required")
    results=[]
    for j,path in enumerate(paths):
        original_phy=Phylo.read(path,"newick")
        ids={t.name for t in original_phy.get_terminals()}
        if len(ids)!=21 or "OUTGROUP_saff" not in ids:
            raise ValueError("Original independently reconstructed twenty-taxon tree invalid")
        original_phy.root_with_outgroup("OUTGROUP_saff")
        panel=make_panel("EastAsia20_S1_8",frame,original_phy)
        nonphylo=neutral_label_null(panel)
        obs=panel["observed"]
        root_cases={}
        for k,root in enumerate(("erect","stationary")):
            null,meta=phylo_null(panel,root,20261023+501*j+81*k,draws)
            root_cases[root]={
                "fitted_Mk_rate_per_source_substitution_length":meta["fitted_Mk_rate_tree_units"],
                "n_conditional_prevalence_matched_null":len(null),
                "one_sided_positive_H1_p":((null>=obs-1e-12).sum(axis=0)+1).astype(float).tolist(),
                "p_denominator":draws+1,
            }
            root_cases[root]["one_sided_positive_H1_p"]=[
                (x)/(draws+1) for x in root_cases[root]["one_sided_positive_H1_p"]]
        results.append({
            "topology":path.name,
            "original_NA_like_source_tree_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "n_matched_S1":len(frame),"n_D":D_COUNT,"n_U":U_COUNT,
            "mean_BIO12_D":float(panel["rain"][panel["y"]==1].mean()),
            "mean_BIO12_U":float(panel["rain"][panel["y"]==0].mean()),
            "mean_BIO12_D_minus_U":float(panel["rain"][panel["y"]==1].mean()-panel["rain"][panel["y"]==0].mean()),
            "observed_standardized_BIO12_and_temperature_residual":obs.tolist(),
            "naive_exact_assignments":len(nonphylo),
            "naive_H1_one_sided_positive_p":(nonphylo>=obs-1e-12).mean(axis=0).tolist(),
            "phylo_Mk":root_cases,
        })
    return {
        "schema":"azami_T2_H1_EA20_six_AU_tree_8original_S1_taxa_v1",
        "source_EAzami_commit":SOURCE_EAZAMI_COMMIT,
        "source_recovered_artifact_digest":SOURCE_ARTIFACT_SHA256,
        "source_morphological_crosswalk_git_blob":ORIGINAL_SOURCE_CROSSWALK_GIT_BLOB,
        "actual_source_N_AU_topologies":N_TOPOLOGIES,
        "actual_source_independent_nuclear_tips":N_NUCLEAR_TIPS,
        "exact_matched_original_photo_taxa":N_MATCHED_S1,
        "matched_taxon_names":sorted(SOURCE_NAME_REQUIRE),
        "outcomes":results,
        "boundaries":{
            "alternative_8_overlaps_original_Japan12":True,
            "do_not_add_as_eight_independent_taxa_to_japan12":True,
            "six_topology_scenarios_not_replicate_adaptive_origins":True,
            "not_full_66_taxa_phylogeny":True,
            "no_causal_rain_protection_or_ancestral_state_reconstruction":True,
            "posthoc_sample_coverage_sensitivity":True,
        }
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original251",type=Path,required=True)
    p.add_argument("--author-crosswalk",type=Path,required=True)
    p.add_argument("--au-tree-directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--draws",type=int,default=30000)
    args=p.parse_args()
    report=assess(args.original251,args.author_crosswalk,args.au_tree_directory,args.draws)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    a=report["outcomes"]
    print("T2_H1_EA20_TRUE_NUCLEAR_ALTERNATIVE_ORIGINAL_S1="+json.dumps({
        "n":N_MATCHED_S1,"D":D_COUNT,
        "rain_D_minus_U":a[0]["mean_BIO12_D_minus_U"],
        "root_sensitivity_ranges":{
            mode:[
                float(min(x["phylo_Mk"][mode]["one_sided_positive_H1_p"][0] for x in a)),
                float(max(x["phylo_Mk"][mode]["one_sided_positive_H1_p"][0] for x in a))]
            for mode in ("erect","stationary")
        },
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
