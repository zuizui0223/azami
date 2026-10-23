"""Prepare a genuine observed-only 111-source Cirsium nucleotide ML input.

Identical original Compositae1061 reference copies were fixed a priori for
118 author QC loci based on >=80% of 111 public source taxa having 3x/25%
Q20 read-depth breadth. Per-site base calls require >=3 Q20 reads and 90%
majority from first120k source reads. Uncalled bases remain N.

For a second, independently set completeness check, only same-copy
coordinates called in >=80% of individually certified source samples
contribute. Locus chosen if >=100 such sites; taxon chosen for a technical
unrooted tree if >=25% of retained base positions are observed. The single
source accession SRR30887230 with conflicting published species labels is
retained in the ALL111 evidence registry, but never automatically promoted
to the strictly named phylogenetic tree. No genotypes are verified as single
copy orthologs; this is not an accepted species tree / BM OU input.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path

FULL=111
LOCI=118
NUCLEOTIDE_QUORUM=.80
MIN_SITES_PER_LOCUS=100
MIN_SAMPLE_FRACTION=.25
MIN_TAXA_FOR_TECHNICAL_ML=60
MIN_LOCI_FOR_TECHNICAL_ML=30
MIN_BASES_FOR_TECHNICAL_ML=5000
MIN_VARIABLE_SITES_FOR_TECHNICAL_ML=50
CONFLICTING_RUN={"SRR30887230"}


def source_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError("Missing source-pinned original 111 tip table")
    return rows


def source_targets(samples,targets):
    if len(samples)!=FULL or len({z["run"] for z in samples})!=FULL:
        raise ValueError("Original Q20 111 source sample inventory changed")
    if len(targets)!=LOCI or len({z["locus"] for z in targets})!=LOCI:
        raise ValueError("Predeclared 118 author target references missing")
    wanted={z["locus"]:z["original_compositae_reference"] for z in targets}
    if any(not name.endswith("-"+gene) for gene,name in wanted.items()):
        raise ValueError("Reference copy and original nuclear locus disagree")
    return {z["run"]:z for z in samples},wanted


def read_receipts(root,samples,targets):
    index,refs=source_targets(samples,targets)
    paths=sorted(Path(root).rglob("shard_*.json"))
    if len(paths)!=10:
        raise ValueError("All 10 independently run empirical 120k nucleotide shards required")
    scopes=[]
    allout={}
    for file in paths:
        stage=json.loads(file.read_text(encoding="utf-8"))
        if stage.get("n_shards")!=10 or stage.get("no_original_FASTQ_or_BAM_preserved") is not True:
            raise ValueError("Bad original public read metadata boundary")
        if stage["failed"] or stage["passed"]!=stage["requested"]:
            raise ValueError("Source nucleotide site calls incomplete at one or more runs")
        scopes.append(stage["shard"])
        for r in stage["outcomes"]:
            run=r["run"]
            if run not in index or run in allout:
                raise ValueError("Unrecognized or duplicate original public run in base-call receipts")
            src=index[run]
            for k in ("biosample","accepted_key","phylogenetic_tip","taxon_name"):
                if r[k]!=src[k]:
                    raise ValueError("Original 111 SRA physical source/voucher mismatch")
            if (r["status"]!="real_source111_120k_Q20_identical_copy_base_calls" or
                r["source_public_ENA_paired_read_count"]!=120000 or
                r["n_predeclared_80percent_reference_loci"]!=LOCI or
                set(r["observed_original_reference_q20_site_calls"])!=set(refs)):
                raise ValueError("Real 118 source reference-copy DNA receipts incomplete")
            if bool(r["source_run_label_conflict_unresolved"])!=(run in CONFLICTING_RUN):
                raise ValueError("Scientific name conflict source exception changed")
            if not r["uncalled_bases_not_filled_by_reference"] or r["original_single_copy_orthology_verified"]:
                raise ValueError("Claimed true orthology or imputed reference at missing DNA site")
            for locus,rec in r["observed_original_reference_q20_site_calls"].items():
                if rec["reference"]!=refs[locus]:
                    raise ValueError("Original same-copy target changed across taxa")
                seen=set()
                for pos,base,depth,fraction in rec["q20_majority_bases_1based"]:
                    if (pos in seen or not 1<=pos<=rec["reference_length_nt"]
                        or base not in "ACGT" or depth<3 or not .9<=fraction<=1):
                        raise ValueError("Invalid genuine DNA call, duplicate coordinate or base quality")
                    seen.add(pos)
            allout[run]=r
    if set(scopes)!=set(range(10)) or set(allout)!=set(index):
        raise ValueError("Original all111 120k public run batches not complete")
    return allout,refs


def site_and_taxa_alignment(data,refs):
    ordered=sorted(data.values(),key=lambda z:z["run"])
    strict=[x for x in ordered if not x["source_run_label_conflict_unresolved"]]
    n=len(strict)
    need=math.ceil(n*NUCLEOTIDE_QUORUM)
    if n<MIN_TAXA_FOR_TECHNICAL_ML:
        raise ValueError("Not enough provenance-unambiguous S1 source taxa")
    sequences={r["phylogenetic_tip"]:[] for r in strict}
    locus_result=[]
    all_sites=0
    informative=0
    variants=0
    # Fixed-locus logic: do not choose loci using observed topology or phenotype.
    for locus,reference in sorted(refs.items()):
        observed={}
        coords=Counter()
        for r in strict:
            src=r["observed_original_reference_q20_site_calls"][locus]
            if src["reference"]!=reference:
                raise ValueError("Source sample sequence was mapped to wrong original reference")
            calls={int(pos):base for pos,base,_dp,_fraction in src["q20_majority_bases_1based"]}
            if len(calls)!=len(src["q20_majority_bases_1based"]):
                raise ValueError("Duplicate observed source base calls")
            observed[r["phylogenetic_tip"]]=calls
            coords.update(calls)
        positions=sorted(pos for pos,k in coords.items() if k>=need)
        valid=len(positions)>=MIN_SITES_PER_LOCUS
        var=info=0
        if valid:
            for pos in positions:
                alleles=Counter(observed[s["phylogenetic_tip"]].get(pos,"N") for s in strict)
                alleles.pop("N",None)
                if len(alleles)>1:
                    var+=1
                if sum(k>=2 for k in alleles.values())>=2:
                    info+=1
            first=all_sites+1
            all_sites+=len(positions)
            for r in strict:
                key=r["phylogenetic_tip"]
                sequences[key].append("".join(observed[key].get(pos,"N") for pos in positions))
            informative+=info
            variants+=var
        else:
            first=0
        locus_result.append({
            "locus":locus,"original_fixed_reference_copy":reference,
            "required_source_taxa_with_real_base_at_one_position":need,
            "n_sites_called_by_80pct_source_taxa":len(positions),
            "n_sites_retained_in_phylo_alignment":len(positions) if valid else 0,
            "non_imputed_variable_sites":var,
            "non_imputed_parsimony_informative_sites":info,
            "alignment_start":first,
            "alignment_end":all_sites if valid else 0,
            "site_selection_only_coverage_not_phenotype_or_tree":True,
            "true_single_copy_orthology_not_confirmed":True,
        })
    concat={k:"".join(v) for k,v in sequences.items()}
    if any(len(seq)!=all_sites for seq in concat.values()):
        raise ValueError("Source genuine site and locus alignment length mismatch")
    sampled=[]
    technical={}
    for r in strict:
        seq=concat[r["phylogenetic_tip"]]
        count=len(seq)-seq.count("N")
        fraction=count/all_sites if all_sites else 0
        item={
            "phylogenetic_tip":r["phylogenetic_tip"],
            "accepted_key":r["accepted_key"],"taxon_name":r["taxon_name"],
            "run":r["run"],"biosample":r["biosample"],
            "available_Q20_sites":count,
            "fraction_observed_same_copy_DNA_sites":fraction,
            "technical_input_min25pct_covered":bool(fraction>=MIN_SAMPLE_FRACTION),
            "provenance_taxonomic_voucher_concept_unverified":True,
        }
        sampled.append(item)
        if fraction>=MIN_SAMPLE_FRACTION:
            technical[r["phylogenetic_tip"]]=seq
    for r in ordered:
        if r["source_run_label_conflict_unresolved"]:
            sampled.append({
                "phylogenetic_tip":r["phylogenetic_tip"],
                "accepted_key":r["accepted_key"],"taxon_name":r["taxon_name"],
                "run":r["run"],"biosample":r["biosample"],
                "available_Q20_sites":None,
                "fraction_observed_same_copy_DNA_sites":None,
                "technical_input_min25pct_covered":False,
                "excluded_from_primary_due_original_taxon_label_conflict":True,
                "provenance_taxonomic_voucher_concept_unverified":True,
            })
    status={
        "schema":"azami_v4_source111_fixed118_realQ20_site_alignment_20261023_v1",
        "source_full_S1_244_accepted_keys_kept_in_original_registry":True,
        "source_original_S1_251_named_units_kept_in_original_registry":True,
        "original_real120k_SRA_source_runs":len(ordered),
        "original_source_label_conflict_run_count":len(ordered)-n,
        "original_unambiguous_taxon_sample_count":n,
        "original_author241_candidates":241,
        "source_fixed_same_original_reference_loci":LOCI,
        "same_copy_source_site_completeness_threshold":NUCLEOTIDE_QUORUM,
        "source_number_for_individual_site":need,
        "candidate_loci_with_100plus_matching_real_sites":
            sum(x["n_sites_retained_in_phylo_alignment"]>0 for x in locus_result),
        "primary_technical_alignment_real_nucleotide_sites":all_sites,
        "primary_technical_alignment_variable_positions":variants,
        "primary_technical_alignment_informative_positions":informative,
        "technical_source_taxa_with25pct_sites":len(technical),
        "initial_all111_alignment_created_from_real_DNA_not_reference_imputation":True,
        "ML_readiness":(len(technical)>=MIN_TAXA_FOR_TECHNICAL_ML and
            sum(x["n_sites_retained_in_phylo_alignment"]>0 for x in locus_result)
            >=MIN_LOCI_FOR_TECHNICAL_ML and
            all_sites>=MIN_BASES_FOR_TECHNICAL_ML and
            variants>=MIN_VARIABLE_SITES_FOR_TECHNICAL_ML),
        "minimum_predeclared_ML_taxa":MIN_TAXA_FOR_TECHNICAL_ML,
        "minimum_predeclared_ML_loci":MIN_LOCI_FOR_TECHNICAL_ML,
        "minimum_predeclared_ML_bases":MIN_BASES_FOR_TECHNICAL_ML,
        "minimum_predeclared_ML_variable":MIN_VARIABLE_SITES_FOR_TECHNICAL_ML,
        "unrooted_technical_tree_only_no_validated_species_tree":True,
        "biological_accepted_taxon_concepts_or_orthologs_validated":False,
        "T2_nodding_origins_or_rain_adaptation_supported_by_this_preflight":False,
        "missing_S1_species_not_reclassified_as_no_gene_or_erect":True,
        "nucleotides_are_from_partial_120k_read_prefix_not_full_SRA":True,
    }
    return sampled,locus_result,concat,technical,status


def write_csv(path,records):
    with Path(path).open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=sorted(set().union(*(z.keys() for z in records))),
                         extrasaction="raise")
        w.writeheader();w.writerows(records)


def fasta(path,seq):
    with Path(path).open("w",encoding="utf-8") as f:
        for tip,bases in sorted(seq.items()):
            f.write(">"+tip+"\n")
            for i in range(0,len(bases),80):
                f.write(bases[i:i+80]+"\n")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples",type=Path,required=True)
    p.add_argument("--targets",type=Path,required=True)
    p.add_argument("--receipts",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    samples=source_csv(a.samples)
    targets=source_csv(a.targets)
    data,refs=read_receipts(a.receipts,samples,targets)
    audited,loci,full,technical,summary=site_and_taxa_alignment(data,refs)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"Q20_111_SRA_SOURCE_SAMPLE_AND_MISSINGNESS.csv",audited)
    write_csv(a.out_dir/"Q20_118_FIXED_LOCUS_REAL_SITE_OCCUPANCY.csv",loci)
    if full:
        fasta(a.out_dir/"REAL111_Q20_ALL_PROVENANCE_CLEAN_ALIGNED_NUCLEOTIDES.fasta",full)
    if technical:
        fasta(a.out_dir/"REAL_SOURCE_Q20_TECHNICAL_ML_PROVISIONAL_NOT_ORTHOLOGS.fasta",
              technical)
    (a.out_dir/"REAL111_NUCLEOTIDE_SITE_PHYLOGENY_INPUT_GATE.json").write_text(
        json.dumps(summary,sort_keys=True,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
        encoding="utf-8")
    print("T2_REAL111_ORIGINAL118_Q20_SOURCE_SITES="+json.dumps(summary,
          sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
