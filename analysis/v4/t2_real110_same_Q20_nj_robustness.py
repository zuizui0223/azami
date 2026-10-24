"""Independent algorithmic sensitivity: observed-only 110-tip Q20 nuclear NJ.

Use the *identical, source-hash checked* 108-locus/29310-base original111
mask, with one post hoc low-coverage sample excluded. Pairwise distances are
calculated only where BOTH public BioSamples have observed Q20 nucleotides.
Missing bases never become the ancestral/reference allele. Neighbour-joining
is an unrooted technical crosscheck for IQ-TREE, NOT a dated or orthology-
validated global Cirsium species tree, evolutionary radiation model or BM/OU.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

N_TIPS=110
SOURCE_TIPS_BEFORE=111
LOCUS_COUNT=108
SITE_COUNT=29310
MIN_SOURCE_Q20_SITES=1000
MIN_PAIRED_CALLABLE=500
EXCLUDED="WCVP206740970"


def read_fasta(path):
    out={}
    current=None
    for raw in Path(path).read_text(encoding="ascii").splitlines():
        row=raw.strip()
        if not row: continue
        if row.startswith(">"):
            label=row[1:].split()[0]
            if label in out or not label.startswith("WCVP"):
                raise ValueError("Duplicate or non-source original WCVP nuclear tip")
            out[label]=""
            current=label
        else:
            if current is None or set(row)-set("ACGTN"):
                raise ValueError("Invalid Q20 observed-only public nuclear sequences")
            out[current]+=row
    if len(out)!=N_TIPS or len({len(x) for x in out.values()})!=1 or (
        len(next(iter(out.values())))!=SITE_COUNT):
        raise ValueError("Original110 source posthoc fixed-mask FASTA shape drift")
    if EXCLUDED in out:
        raise ValueError("Original low-Q20 specimen illegally retained")
    if min(len(x)-x.count("N") for x in out.values())<MIN_SOURCE_Q20_SITES:
        raise ValueError("Unexpected second low-Q20 original tip")
    return out


def pairwise_observed_pdist(seqs,threshold=MIN_PAIRED_CALLABLE):
    ids=sorted(seqs)
    arrays=[np.frombuffer(seqs[tip].encode("ascii"),dtype="S1") for tip in ids]
    matrix=[];scores=[]
    for i,one in enumerate(arrays):
        row=[]
        for j,two in enumerate(arrays[:i]):
            usable=(one!=b"N")&(two!=b"N")
            shared=int(np.count_nonzero(usable))
            if shared<threshold:
                raise ValueError("Too few joint observed source Q20 bases for this original pair")
            differences=int(np.count_nonzero(one[usable]!=two[usable]))
            dist=differences/shared
            if not 0<=dist<=1:
                raise ValueError("Pairwise original DNA distance invalid")
            scores.append({"tip_A":ids[j],"tip_B":ids[i],
                           "n_shared_Q20_positions":shared,
                           "n_differing_Q20_positions":differences,
                           "pairwise_deletion_p_distance":dist})
            row.append(dist)
        row.append(0.)
        matrix.append(row)
    if len(scores)!=(len(ids)*(len(ids)-1)//2):
        raise ValueError("Not all original Q20 pairwise sample distances computed")
    return ids,matrix,scores


def check_manifest(path):
    r=json.loads(Path(path).read_text(encoding="utf-8"))
    if (r.get("posthoc_sensitivity_taxa")!=N_TIPS or
        r.get("original_all111_Q20_successful_real_SRA_BioSamples")!=SOURCE_TIPS_BEFORE or
        r.get("original111_technical_gate_passed") is not False or
        r.get("original_111_source_108_retained_loci")!=LOCUS_COUNT or
        r.get("observed_Q20_reference_positions_reused_unchanged")!=SITE_COUNT or
        r.get("original_111_predeclared_89_of111_quorum_unmodified")!=89 or
        r.get("original_111_reference_position_selection_recomputed_on110") is not False or
        r.get("original_S1_WCVP_master_accepted_keys")!=244 or
        r.get("original_S1_source_name_operational_units")!=251 or
        r.get("source111_single_species_Q20_outlier",{}).get("tip")!=EXCLUDED or
        r.get("copy_number_paralogy_polyploidy_review_completed") is not False):
        raise ValueError("Original111 predeclared failed gate or posthoc110 source mask altered")
    return r


def audit_nj(seq,manifest):
    ids,matrix,dist=pairwise_observed_pdist(seq)
    tree=DistanceTreeConstructor().nj(DistanceMatrix(ids,matrix))
    tips={x.name for x in tree.get_terminals()}
    if tips!=set(ids):
        raise ValueError("Neighbour joining lost actual 110 original public BioSamples")
    branches=[float(x.branch_length) for x in tree.find_clades()
              if x.branch_length is not None]
    if len(branches)<N_TIPS-1:
        raise ValueError("No meaningful original real source NJ branch lengths")
    pairs=[x["n_shared_Q20_positions"] for x in dist]
    return tree,dist,{
       "schema":"azami_T2_original111_posthoc110_same_Q20_distance_NJ_sensitivity_v1",
       "status":"REAL110_ORIGINAL111_MASK_INDEPENDENT_ALGORITHM_TECHNICAL_NJ_NOT_SPECIES_PHYLOGENY",
       "source_public_BioSamples":SOURCE_TIPS_BEFORE,
       "source_posthoc_Q20_BioSamples":N_TIPS,
       "original_loci":LOCUS_COUNT,
       "original_Q20_reference_positions":SITE_COUNT,
       "Q20_site_mask_refitted_on_110":False,
       "original_failed111_quality_gate_relabelled_as_success":False,
       "n_original_public_BioSample_pairs":len(dist),
       "joint_observed_Q20_positions_min":min(pairs),
       "joint_observed_Q20_positions_median":statistics.median(pairs),
       "joint_observed_Q20_positions_max":max(pairs),
       "technical_NJ_positive_branches":sum(x>0 for x in branches),
       "technical_NJ_negative_branches":sum(x<0 for x in branches),
       "all_244_original_accepted_keys_and_251_named_units_preserved":True,
       "nuclear_copy_number_single_orthology_verified":False,
       "global_EAzami294_true_tree_ready":False,
       "nodding_ancestral_gains_or_rain_adaptation_inferred":False,
       "quartets_are_not_independent_replicates":True,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta",type=Path,required=True)
    p.add_argument("--original-receipt",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    receipts=check_manifest(a.original_receipt)
    dna=read_fasta(a.fasta)
    tree,dist,summary=audit_nj(dna,receipts)
    summary["source_fasta_sha256"]=hashlib.sha256(a.fasta.read_bytes()).hexdigest()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    Phylo.write(tree,str(a.out_dir/"REAL110_SOURCE_NUCLEAR_ORIGINALMASK_UNROOTED_NJ.nwk"),"newick")
    with (a.out_dir/"REAL110_PUBLIC_SOURCE_Q20_PAIRWISE_DISTANCES.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(dist[0]),extrasaction="raise")
        w.writeheader();w.writerows(dist)
    (a.out_dir/"REAL110_NJ_SOURCE_ROBUSTNESS_SUMMARY.json").write_text(
        json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n",
        encoding="utf-8")
    print("REAL_SOURCE110_DIFFERENT_METHOD_NJ_DISTANCES="+json.dumps(summary,
          sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
