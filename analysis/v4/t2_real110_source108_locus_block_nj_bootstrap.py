"""Whole-locus resampling diagnostic on real110 S1-linked nuclear read prefixes.

Uses SAME immutable 111-source original 89/111 masked Q20 nucleotide positions
(108 reference-copy candidates, 29310 bases, 110 observed SRA BioSamples
after one named sparse specimen is excluded). Per-locus pairwise counts are
aggregated; 200 NJ trees are inferred after sampling 108 whole candidate
loci WITH replacement. Compare their unrooted splits to original same-source
IQ-TREE GTR+G4 ML and NJ.

This is NJ *locus-block* resampling, not IQ-TREE UFBoot, ML locus bootstrap,
species-tree posterior, orthology, climate-response or independent nodding
gain evidence. Several candidate loci may be paralogous/homeologous.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import statistics

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix,DistanceTreeConstructor

from analysis.v4.t2_real110_ml_vs_nj_same_q20_topology import split_set
from analysis.v4.t2_real110_same_Q20_nj_robustness import read_fasta
from analysis.v4.t2_real110_posthoc_one_sparse_source_sensitivity import (
    ORIGINAL_SHARED_SITES,RECOVERED_LOCI,AFTER_SAMPLES,PREDECLARED_111_QUORUM
)

SEED=20261023
REPLICATES=200
MIN_PAIR_Q20_SHARED=500


def csv_rows(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError("Missing public source nuclear gene partitions/distances")
    return rows


def partitions_check(rows):
    if len(rows)!=RECOVERED_LOCI or len({r["original_locus"] for r in rows})!=RECOVERED_LOCI:
        raise ValueError("Missing original 108 reference-copy Q20 gene units")
    ordered=sorted(rows,key=lambda r:int(r["start_1based"]))
    prev=0
    intervals=[]
    for x in ordered:
        a=int(x["start_1based"])-1
        b=int(x["end_1based"])
        if a!=prev or b<=a or x["not_single_copy_Cirsium_orthology"].lower()!="true":
            raise ValueError("Original whole-locus boundaries or source orthology caveat changed")
        intervals.append((a,b,x["original_locus"]))
        prev=b
    if prev!=ORIGINAL_SHARED_SITES:
        raise ValueError("Original 29310 masked Q20 nuclear positions changed")
    return intervals


def source_pair_counts(dna,intervals,original_distances):
    names=sorted(dna)
    if len(names)!=AFTER_SAMPLES:
        raise ValueError("Original source 110 real nuclear BioSamples missing")
    arr=np.stack([np.frombuffer(dna[x].encode("ascii"),dtype="S1") for x in names])
    if arr.shape!=(AFTER_SAMPLES,ORIGINAL_SHARED_SITES):
        raise ValueError("Q20 110x29310 observed source base alignment altered")
    pairs=[(j,i) for i in range(len(names)) for j in range(i)]
    n_pairs=len(pairs)
    source_map={(r["tip_A"],r["tip_B"]):r for r in original_distances}
    if len(source_map)!=n_pairs:
        raise ValueError("Missing original 5995 same-source observed Q20 pair distances")
    matches=np.zeros((len(intervals),n_pairs),dtype=np.uint16)
    compared=np.zeros_like(matches)
    for k,(begin,end,_name) in enumerate(intervals):
        locus=arr[:,begin:end]
        for p,(j,i) in enumerate(pairs):
            one=locus[j];two=locus[i]
            valid=(one!=b"N")&(two!=b"N")
            n=int(np.count_nonzero(valid))
            compared[k,p]=n
            matches[k,p]=int(np.count_nonzero(one[valid]!=two[valid]))
    paired=compared.astype(np.uint32).sum(axis=0)
    variants=matches.astype(np.uint32).sum(axis=0)
    if int(np.min(paired))<MIN_PAIR_Q20_SHARED:
        raise ValueError("Original complete site-pair source data below declared Q20 overlap gate")
    max_error=0.
    for p,(j,i) in enumerate(pairs):
        key=(names[j],names[i])
        if key not in source_map:
            raise ValueError("Wrong original source BioSample pair key")
        r=source_map[key]
        if (int(r["n_shared_Q20_positions"])!=int(paired[p]) or
            int(r["n_differing_Q20_positions"])!=int(variants[p])):
            raise ValueError("Per-locus reconstruction of original Q20 distances differs from published NJ receipt")
        observed=float(r["pairwise_deletion_p_distance"])
        max_error=max(max_error,abs(float(variants[p]/paired[p])-observed))
    if max_error>1e-12:
        raise ValueError("Real reference-copy pair distances changed at Q20 source positions")
    return names,compared,matches,pairs,max_error


def distance_matrix(names,pairs,pdist):
    rows=[[] for _ in names]
    for (j,i),dist in zip(pairs,pdist):
        rows[i].append(float(dist))
    for i,row in enumerate(rows):
        if len(row)!=i or any(not 0<=x<=1 for x in row):
            raise ValueError("Original public nuc Q20 distances invalid")
        row.append(0.0)
    return DistanceMatrix(names,rows)


def bootstrap_analysis(names,counts,differences,pairs,ml_tree,nj_tree,n=REPLICATES,seed=SEED):
    if n<100:
        raise ValueError("No exploratory tiny resample on locus-block phylo uncertainty")
    origin=set(names)
    ml_names,ml_splits=split_set(Phylo.read(str(ml_tree),"newick"))
    nj_names,nj_splits=split_set(Phylo.read(str(nj_tree),"newick"))
    if ml_names!=origin or nj_names!=origin or len(ml_splits)!=len(names)-3 or len(nj_splits)!=len(names)-3:
        raise ValueError("Real110 ML and NJ trees do not share the exact same unrooted biological source tips")
    rng=np.random.default_rng(seed)
    support=Counter()
    boot_pair_min=[]
    overlaps=[]
    for r in range(n):
        sampled=rng.integers(0,len(counts),size=len(counts))
        observed=counts[sampled].astype(np.uint32).sum(axis=0)
        changes=differences[sampled].astype(np.uint32).sum(axis=0)
        if np.min(observed)<MIN_PAIR_Q20_SHARED:
            raise ValueError("One locus resample made original Q20 distance unidentifiable")
        pdist=changes/observed
        tree=DistanceTreeConstructor().nj(distance_matrix(names,pairs,pdist))
        the_names,splits=split_set(tree)
        if the_names!=origin or len(splits)!=len(names)-3:
            raise ValueError("Whole-locus resample lost original tree source tip/edge scope")
        support.update(splits)
        boot_pair_min.append(int(np.min(observed)))
        overlaps.append(len(splits&ml_splits))
    by_split=[]
    for item in sorted(ml_splits):
        by_split.append({
           "ML_unrooted_bipartition_small_canonical_side":"|".join(item),
           "source_taxa_on_one_side":len(item),
           "NJ_original_same_mask_split":item in nj_splits,
           "whole_locus_resampled_NJ_count":support[item],
           "whole_locus_resampled_NJ_frequency":support[item]/n,
           "is_ML_IQTREE_UFBoot":False,
           "Q20_taxon_concept_and_homeolog_phylogeny_validated":False,
        })
    frequencies=[x["whole_locus_resampled_NJ_frequency"] for x in by_split]
    return by_split,{
       "schema":"azami_T2_original111_posthoc110_Q20_108locus_200_NJ_block_bootstrap_v1",
       "status":"EMPIRICAL_PUBLIC_Q20_WHOLE_REFERENCE_LOCUS_NJ_RESAMPLING_NOT_GENUINE_ORTHOLOG_SPECIES_TREE",
       "original111_failed_taxon_gate_preserved":True,
       "posthoc110_source_tips":len(names),
       "original_original111_89tip_site_quorum_unchanged":PREDECLARED_111_QUORUM,
       "original_same_copy_Q20_loci":len(counts),
       "original_masked_Q20_sites":ORIGINAL_SHARED_SITES,
       "resampling_unit":"one original author candidate reference-copy locus",
       "n_whole_locus_resampled_NJ_trees":n,
       "seed":seed,
       "source_IQTREE_ML_splits":len(ml_splits),
       "source_NJ_fulldata_splits":len(nj_splits),
       "ML_splits_also_present_in_original_NJ":len(ml_splits&nj_splits),
       "ML_splits_ge75pct_whole_locus_NJ_resamples":sum(x>=.75 for x in frequencies),
       "ML_splits_ge95pct_whole_locus_NJ_resamples":sum(x>=.95 for x in frequencies),
       "median_ML_split_locus_resample_NJ_frequency":statistics.median(frequencies),
       "median_per_resampled_NJ_ML_split_overlap":statistics.median(overlaps),
       "minimum_pair_source_Q20_bases_in_bootstrap":min(boot_pair_min),
       "site_only_bootstrap_or_IQTREE_ML_bootstrap_claimed":False,
       "public_108_loci_independently_verified_single_copy":False,
       "EAzami294_whole_genus_nuclear_species_tree_completed":False,
       "U_to_D_events_or_rainfall_adaptive_radiation_tested":False,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("fasta","partitions","original-nj-pairs","ml-tree","nj-tree","out-dir"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    dna=read_fasta(a.fasta)
    intervals=partitions_check(csv_rows(a.partitions))
    names,coverage,differences,pairs,err=source_pair_counts(
        dna,intervals,csv_rows(a.original_nj_pairs))
    scores,summary=bootstrap_analysis(names,coverage,differences,pairs,a.ml_tree,a.nj_tree)
    summary["source_real_Q20_FASTA_SHA256"]=hashlib.sha256(a.fasta.read_bytes()).hexdigest()
    summary["all_5995_original_NJ_pair_distances_reproduced_with_absolute_error_le1e_12"]=True
    summary["maximum_original_distance_recompute_error"]=err
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"REAL110_ML_BRANCH_SUPPORT_BY_WHOLE_Q20_CANDIDATE_LOCUS_NJ.csv").open(
            "w",newline="",encoding="utf-8") as fp:
        writer=csv.DictWriter(fp,fieldnames=list(scores[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(scores)
    (a.out_dir/"REAL110_WHOLE_REFERENCE_LOCUS_NJ_200_SOURCE_Q20_BOOTSTRAP.json").write_text(
        json.dumps(summary,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    print("REAL110_SOURCE_SAME_Q20_WHOLE_LOCUS_200_NJ_SUPPORT="+json.dumps(summary,
          sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
