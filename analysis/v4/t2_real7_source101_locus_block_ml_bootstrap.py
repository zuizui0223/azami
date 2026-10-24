"""Real seven-taxon whole-locus 300-replicate *ML* bootstrap, separate from site UFBoot.

Resamples, with replacement, 101 entirely empirically Q20 observed candidate
nuclear loci (>=100 same-coordinate sites all seven specimens) in each draw.
Re-estimates GTR+G4 IQ-TREE 2 maximum-likelihood topology, NOT UFBoot or
Newick simulation. Compare 300 locus-level ML replicates with the original
IQ-TREE technical 7-source topology and independent NA63 4-tip quartet.
A gene-locus bootstrap is different from both IQ-TREE site UFBoot and
earlier 300 whole-locus *NJ* bootstrap and is not a species-tree posterior.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import random
import shutil
import subprocess
import time

from Bio import Phylo

from analysis.v4 import t2_real7_author101_iqtree_ml as ml
from analysis.v4 import t2_seven_original110_Q20_nj_diagnostic as nj
from analysis.v4 import t2_real7_NA63_quartet_crosscheck as original

N_DRAWS=300
N_SHARDS=2
SOURCE_SEED=20261023
CANDIDATE_LOCI=101
SOURCE_RUNS=7
N_INTERNAL_SPLITS=4
EXPECTED_TRUE_NA_QUARTET=["Cirsium_flodmanii","Cirsium_tracyi"]


def boot_names_and_loci(source125,chosen,targets,receipts):
    original_rows=ml.read_csv(source125)
    chosen_rows=ml.read_csv(chosen)
    refs=ml.source.targets(ml.read_csv(targets))
    names,records,loci,primary,counts=ml.build_matrices(
         original_rows,chosen_rows,refs,receipts)
    if len(primary)!=CANDIDATE_LOCI or counts["n_primary_observed_only_sites"]!=28726:
        raise ValueError("Actual seven-run 101 source observed loci changed")
    return names,primary,counts


def bootstrap_alignment(names,loci,seed):
    if len(names)!=SOURCE_RUNS or len(loci)!=CANDIDATE_LOCI:
        raise ValueError("Expected seven original run sources and 101 source loci")
    rng=random.Random(seed)
    chosen=[loci[rng.randrange(len(loci))] for _ in range(len(loci))]
    # Each locus is one independent bootstrap block, including invariant bases.
    concatenated=ml.alignment_by_samples(names,chosen)
    if len(set(len(x) for x in concatenated.values()))!=1:
        raise ValueError("Variable observed site count across real source samples")
    return concatenated,chosen


def split_labels(tree):
    splits=nj.splits_unrooted(tree)
    if len(splits)!=N_INTERNAL_SPLITS:
        raise ValueError("Unresolved seven-source unrooted ML tree")
    return sorted("|".join(x) for x in splits)


def infer_replicate(i,names,loci,workdir,iqtree="iqtree2",run=subprocess.run):
    seq,locus_rows=bootstrap_alignment(names,loci,SOURCE_SEED+i)
    folder=Path(workdir)/f"rep_{i:04d}"
    folder.mkdir(parents=True,exist_ok=True)
    try:
        fa=folder/"shared_observed_Q20.fasta"
        fa.write_text(ml.fasta_payload(seq),encoding="ascii")
        prefix=folder/"ml"
        cmd=[iqtree,"-s",str(fa),"-m",ml.MODEL,"-T","1",
             "--seed",str(SOURCE_SEED+i),"-pre",str(prefix)]
        proc=run(cmd,check=True,capture_output=True,timeout=120)
        treefile=folder/"ml.treefile"
        if not treefile.is_file():
            raise ValueError("IQ-TREE did not produce a real seven-SRA ML tree")
        phy=Phylo.read(str(treefile),"newick")
        if {z.name for z in phy.get_terminals()}!=set(names):
            raise ValueError("One original source sample vanished in bootstrapped ML tree")
        return {
           "bootstrap_locus_draw_id":i,
           "source_ML_loci_with_replacement":len(locus_rows),
           "n_original_unique_candidate_loci_in_draw":
                len({z["locus"] for z in locus_rows}),
           "n_source_complete_case_nucleotide_sites":len(next(iter(seq.values()))),
           "unrooted_ML_split_labels":split_labels(phy),
           "original_7_SRA_source_SNPs_only":True,
           "full_locus_sequences_reconstructed":False,
        }
    finally:
        shutil.rmtree(folder,ignore_errors=True)


def run_shard(names,loci,shard_id,shard_count,workdir,runner=infer_replicate,draws=N_DRAWS):
    if not(0<=shard_id<shard_count<=N_SHARDS) or draws<100:
        raise ValueError("Invalid source real ML locus-bootstrap shard settings")
    results=[]
    for i in range(shard_id,draws,shard_count):
        trial=runner(i,names,loci,workdir)
        if trial["bootstrap_locus_draw_id"]!=i or \
           len(trial["unrooted_ML_split_labels"])!=N_INTERNAL_SPLITS:
            raise ValueError("Original source read ML bootstrap replicate missing")
        results.append(trial)
        if len(results)%25==0:
            print("REAL_Q20_SOURCE_ML_LOCUS_BOOTSTRAP_PROGRESS="+str(len(results)),flush=True)
    return {
      "schema":"azami_real7_source_101_locus_ML_block_bootstrap_shard_v1",
      "shard_id":shard_id,"shard_count":shard_count,
      "total_predeclared_locus_bootstrap_replicates":draws,
      "n_original_real_SRA_samples":len(names),
      "n_original_101_candidate_loci":len(loci),
      "n_results":len(results),"source_only_no_simulated_bases":True,
      "not_accepted_species_tree":True,
      "replicates":results,
    }


def aggregate(shards,original_ml_tree,independent_na_tree):
    if len(shards)!=N_SHARDS or {x["shard_id"] for x in shards}!={0,1}:
        raise ValueError("Missing paired source ML block bootstrap shard")
    results={}
    for r in shards:
        if (r["n_original_real_SRA_samples"]!=SOURCE_RUNS or
            r["n_original_101_candidate_loci"]!=CANDIDATE_LOCI or
            not r["source_only_no_simulated_bases"]):
            raise ValueError("Source real seven-runs Q20 bootstrap contract broken")
        for x in r["replicates"]:
            i=x["bootstrap_locus_draw_id"]
            if i in results or i%N_SHARDS!=r["shard_id"]:
                raise ValueError("Duplicate or wrong original locus-block bootstrap index")
            if not x["original_7_SRA_source_SNPs_only"] or \
               x["full_locus_sequences_reconstructed"]:
                raise ValueError("Bootstrap falsely claims simulated or ortholog gene data")
            results[i]=x
    if set(results)!=set(range(N_DRAWS)):
        raise ValueError("At least one of the 300 original-locus ML repetitions lost")
    primary_splits=set(split_labels(original_ml_tree))
    count=Counter(s for r in results.values() for s in r["unrooted_ML_split_labels"])
    result_rows=[]
    for label in sorted(primary_splits):
        result_rows.append({
            "source_ML_primary_unrooted_split":label,
            "n_300_real_Q20_whole_locus_ML_replicates_with_split":count[label],
            "ML_whole_locus_bootstrap_fraction":count[label]/N_DRAWS,
            "unit_of_resampling":"original_101_candidate_nuclear_loci_not_independent_sites",
            "not_species_tree_posterior":True,
        })
    quartets=[]
    for i,r in results.items():
        # Bootstrap replicate metadata lists only split labels. For a quartet,
        # determine whether F+T vs S+V split can be induced from the 7-tip
        # bipartitions; avoids a second IQ-TREE inference.
        fq={"Cirsium_flodmanii","Cirsium_tracyi"}
        sq={"Cirsium_scopulorum","Cirsium_vulgare"}
        names=set(z.name for z in original_ml_tree.get_terminals())
        matching=False
        for s in r["unrooted_ML_split_labels"]:
            side=set(s.split("|"))
            complement=names-side
            if (fq<=side and sq<=complement) or (fq<=complement and sq<=side):
                matching=True
                break
        quartets.append(matching)
    published=original.quartet_split(independent_na_tree)
    if published!=EXPECTED_TRUE_NA_QUARTET:
        raise ValueError("Author independent NA63 quartet no longer matches original proof")
    return {
       "schema":"azami_v4_source7_101_candidate_Q20_locus_block_ML_bootstrap_300_v1",
       "status":"REAL_PUBLIC_SRA_ML_LOCUS_BLOCK_SUPPORT_TECHNICAL_ONLY",
       "n_original_public_SRA_samples":SOURCE_RUNS,
       "n_shared_candidate_nuclear_loci":CANDIDATE_LOCI,
       "n_total_independently_reestimated_ML_replicates":N_DRAWS,
       "fixed_seed":SOURCE_SEED,
       "ML_inferred_original_Q20_observed_bases_only":True,
       "quartet_scored_exact_four_against_independent_author_NA63":True,
       "n_block_ML_replicates_consistent_with_published_NA63_quartet":sum(quartets),
       "ML_locus_block_published_NA63_quartet_fraction":sum(quartets)/N_DRAWS,
       "original_101_ML_splits_n":len(primary_splits),
       "whole_locus_ML_branch_frequency":result_rows,
       "uncertainty":{
          "not_site_UFBOOT":True,
          "not_independent_data_against_same_7_SRA_read_source":True,
          "orthology_and_homeology_not_verified":True,
          "species_phylogeny_rooted_or_dated":False,
          "head_posture_transition_or_rainfall_adaptation_estimated":False,
          "global_EAzami294_species_phylogeny_completed":False,
       }
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode",choices=("shard","aggregate"),required=True)
    p.add_argument("--sra-source",type=Path)
    p.add_argument("--chosen",type=Path)
    p.add_argument("--targets",type=Path)
    p.add_argument("--receipts",type=Path)
    p.add_argument("--shard",type=int)
    p.add_argument("--scratch",type=Path)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--shard-dir",type=Path)
    p.add_argument("--original-ml",type=Path)
    args=p.parse_args()
    if args.mode=="shard":
        if not all([args.sra_source,args.chosen,args.targets,args.receipts,args.scratch]) or args.shard is None:
            p.error("Require original public SRA run and locus receipts for ML block bootstrap")
        names,loci,_=boot_names_and_loci(args.sra_source,args.chosen,args.targets,args.receipts)
        output=run_shard(names,loci,args.shard,N_SHARDS,args.scratch)
    else:
        if not args.shard_dir or not args.original_ml:
            p.error("Require original seven-SRA ML tree and paired actual bootstrap shards")
        paths=sorted(args.shard_dir.rglob("blockML_*.json"))
        if len(paths)!=N_SHARDS:
            raise ValueError("Exactly two original source ML block bootstrap shards required")
        shards=[json.loads(path.read_text()) for path in paths]
        observed=Phylo.read(str(args.original_ml),"newick")
        author_tree=original.load_na63()
        output=aggregate(shards,observed,author_tree)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    if args.mode=="shard":
        print("REAL_SRA_Q20_101_LOCUS_ML_BLOCK_SHARD="+json.dumps({
           "shard":args.shard,"n":output["n_results"]}),flush=True)
    else:
        print("REAL_SRA_Q20_101_LOCUS_ML_BLOCK_300_COMPLETE="+json.dumps({
          "frequency":output["whole_locus_ML_branch_frequency"],
          "NA63_quartet_fraction":output["ML_locus_block_published_NA63_quartet_fraction"]
        },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
