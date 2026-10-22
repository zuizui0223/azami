"""Actual seven Cirsium 120k-run Q20 majority-base evidence at identical positions.

Input is exactly 110 original author241 targets with seven real 120k Q20
specimens sharing the *same* original Asteraceae reference copy.
Run samtools/BWA on real public ENA partial FASTQ, then pysam pileup at
mapQ>=20, baseQ>=20, n>=3, >=90% single-base agreement. All other positions
are ABSENT/uncallable, never silently imputed as the Asteraceae reference.

Resulting 7-sample site calls are screening evidence, NOT validated genotypes:
unresolved polyploidy, heterozygosity, reference-copy bias, paralogous mapping,
contamination, allelic dropout, and taxon identity remain review gates.
No rooted phylogeny, phenotype transition or adaptation is inferred.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path
import shutil
import subprocess

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as original
from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as readops
from analysis.v4 import t2_moreyra_all125_author241_site_depth as depth

N_ORIGINAL_RUNS=125
N_DEEP_SAMPLES=7
N_SHARED_TARGETS=110
READ_PAIRS=120000
Q20=20
MIN_SUPPORT=3
MIN_MAJOR_FRACTION=.9
SEVEN={
 "SRR30887180","SRR30887204","SRR30887201",
 "SRR30887197","SRR30887312","SRR30887284","SRR25265665",
}

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError("Missing same-run real public original sequence source")
    return rows

def targets(rows):
    if len(rows)!=241 or len({r["locus"] for r in rows})!=241:
        raise ValueError("Source original 241 no-warning locus names are missing")
    out={}
    for r in rows:
        flag=str(r["deep7_all_Q20_dp3_ge25pct_same_ref"]).lower()
        if flag not in ("true","false"):
            raise ValueError("Frozen all7 one-copy source high-depth flag unspecified")
        if flag!="true":continue
        if str(r["best_reference_copy_identical_all7"]).lower()!="true":
            raise ValueError("Reference copy count contradicts seven original source agreement")
        names=r["original_ref_names"].split("|")
        if len(names)!=1 or not names[0].startswith(("saff-","lett-","sunf-")):
            raise ValueError("Original shared Asteraceae target reference label invalid")
        if not names[0].endswith("-"+r["locus"]):
            raise ValueError("Original locus and single reference do not agree")
        out[r["locus"]]=names[0]
    if len(out)!=N_SHARED_TARGETS or len(set(out.values()))!=N_SHARED_TARGETS:
        raise ValueError("Exactly 110 original same-copy 120k shared targets required")
    return out

def source_sample(original125,chosen,accession):
    if len(original125)!=N_ORIGINAL_RUNS or len({r["run"] for r in original125})!=N_ORIGINAL_RUNS:
        raise ValueError("Source original 125 SRA run/BioSamples not complete")
    if len(chosen)!=N_DEEP_SAMPLES or {x["run"] for x in chosen}!=SEVEN:
        raise ValueError("Seven source-pinned high-depth runs changed")
    src={r["run"]:r for r in original125}
    candidates={r["run"]:r for r in chosen}
    if accession not in SEVEN:
        raise ValueError("Run outside seven source-authenticated nuclear samples")
    row=src.get(accession);s=candidates[accession]
    if row is None or any(row[k]!=s[k] for k in ("accepted_key","biosample")):
        raise ValueError("Source BioSample/accepted taxon changed after source 120k scan")
    name=lambda x:"_".join(str(x).replace(" ","_").split("_")).strip("_")
    if name(s["source_taxon"])!=name(row["S1_source_taxon"]):
        raise ValueError("Original source taxon and selected 7-run name mismatch")
    if int(s["same_run_120k_dp3_ge50pct"])<100:
        raise ValueError("Public seven source coverage class is no longer >=100 loci")
    return row

def majority_base_call(raw_bases):
    bases=[str(x).upper() for x in raw_bases if str(x).upper() in {"A","C","G","T"}]
    if len(bases)<MIN_SUPPORT:
        return None
    common,count=Counter(bases).most_common(1)[0]
    fraction=count/len(bases)
    if fraction<MIN_MAJOR_FRACTION:
        return None
    return common,len(bases),round(fraction,6)

def call_sites_from_pysam(bam,contigs):
    """Never infer reference base at a site without empirical Q20 observations."""
    import pysam
    outcome={}
    with pysam.AlignmentFile(str(bam),"rb") as input_bam:
        known=set(input_bam.references)
        for locus,contig in contigs.items():
            if contig not in known:
                raise ValueError("One of original 110 authentic reference contigs missing")
            length=input_bam.get_reference_length(contig)
            if length<50:
                raise ValueError("Unexpected short original nuclear target")
            rows=[]
            for col in input_bam.pileup(
                contig,0,length,truncate=True,stepper="samtools",
                min_base_quality=Q20,min_mapping_quality=Q20,
                ignore_overlaps=True,ignore_orphans=False,
                flag_filter=0xF04,max_depth=100000):
                bases=[]
                for pile in col.pileups:
                    if pile.is_del or pile.is_refskip or pile.query_position is None:
                        continue
                    read=pile.alignment
                    q=pile.query_position
                    if read.mapping_quality<Q20 or read.is_secondary or \
                       read.is_supplementary or read.is_duplicate or read.is_qcfail:
                        continue
                    if read.query_qualities is None or read.query_qualities[q]<Q20:
                        continue
                    bases.append(read.query_sequence[q])
                called=majority_base_call(bases)
                if called:
                    base,n,consistency=called
                    rows.append([col.reference_pos+1,base,n,consistency])
            if len({rec[0] for rec in rows})!=len(rows):
                raise ValueError("Duplicate source reference-coordinate major-base call")
            outcome[locus]={
                "reference":contig,
                "reference_length_nt":length,
                "n_conservatively_callable_sites":len(rows),
                "q20_majority_bases_1based":rows,
                "not_a_validated_allelic_genotype":True,
            }
    if set(outcome)!=set(contigs):
        raise ValueError("Original source accepted 110 reference targets lost")
    return outcome

def fetch_align_and_call(src,ref,contigs,original241,scratch):
    accession=src["run"]
    folder=Path(scratch)/accession
    folder.mkdir(parents=True,exist_ok=True)
    try:
        mates,_meta=readops.acquire_larger_prefix(src)
        if len(mates)!=2 or any(len(m)!=READ_PAIRS*4 for m in mates):
            raise ValueError("Source public run did not yield exact first120k paired reads")
        # Full original author241 depth check is a baseline identity check;
        # no genotype or variant sites are silently filled by the reference.
        if len(original241)!=241 or not set(contigs).issubset(set(original241)):
            raise ValueError("Original author241 panel must contain 110 shared reference targets")
        stats=readops.map_read_pairs(src,ref,original241,mates,folder)
        if stats.get("input_read_pairs")!=READ_PAIRS:
            raise ValueError("Incorrect original source nested read count")
        calls=call_sites_from_pysam(folder/"aligned.sorted.bam",contigs)
        return {
          "status":"real_source120k_Q20_majority_nucleotide_position_pilot",
          "run":accession,"biosample":src["biosample"],
          "accepted_key":src["accepted_key"],
          "S1_source_taxon":src["S1_source_taxon"],
          "n_original_shared_reference_targets":N_SHARED_TARGETS,
          "n_source_pairs":READ_PAIRS,
          "n_Q20_callable_positions_total":sum(x["n_conservatively_callable_sites"] for x in calls.values()),
          "targets":calls,
          "variant_positions_not_adjudicated":True,
          "heterozygous_sites_masked_not_imputed":True,
          "no_ortholog_alignment_or_nuclear_tree_built":True,
          "original_FASTQ_BAM_deleted_before_artifact":True,
        }
    finally:
        shutil.rmtree(folder,ignore_errors=True)

def validate_single(result,source,contigs):
    if (result["run"]!=source["run"] or
        result["biosample"]!=source["biosample"] or
        result["accepted_key"]!=source["accepted_key"] or
        result["n_source_pairs"]!=READ_PAIRS or
        result["n_original_shared_reference_targets"]!=N_SHARED_TARGETS or
        set(result["targets"])!=set(contigs)):
        raise ValueError("Incomplete original real SRA site-level inference receipt")
    for locus,rec in result["targets"].items():
        if rec["reference"]!=contigs[locus] or not rec["not_a_validated_allelic_genotype"]:
            raise ValueError("Source target locus and unvalidated-allele status drift")
        seen=set()
        for pos,base,depth_num,fraction in rec["q20_majority_bases_1based"]:
            if (pos in seen or not(1<=pos<=rec["reference_length_nt"]) or
                base not in "ACGT" or depth_num<MIN_SUPPORT or
                not(MIN_MAJOR_FRACTION<=fraction<=1)):
                raise ValueError("Malformed original Q20 conservative base-call evidence")
            seen.add(pos)
    return True

def aggregate(source_rows,contigs):
    if len(source_rows)!=N_DEEP_SAMPLES or {x["run"] for x in source_rows}!=SEVEN:
        raise ValueError("All seven independent original SRA runs required")
    by_locus={locus:[] for locus in contigs}
    for row in source_rows:
        if row["n_original_shared_reference_targets"]!=N_SHARED_TARGETS or set(row["targets"])!=set(contigs):
            raise ValueError("Different reference coordinate universes in source samples")
        for locus in contigs:
            r=row["targets"][locus]
            if r["reference"]!=contigs[locus]:
                raise ValueError("Asteraceae source reference copy differs between taxa")
            by_locus[locus].append({int(pos):base for pos,base,_depth,_fraction in r["q20_majority_bases_1based"]})
    locus_summaries=[]
    for locus,seven in by_locus.items():
        intersection=set(seven[0])
        for call in seven[1:]:
            intersection.intersection_update(call)
        shared=sorted(intersection)
        variable=0
        informative=0
        for pos in shared:
            bases=[r[pos] for r in seven]
            counts=Counter(bases)
            if len(counts)>1:variable+=1
            if sum(count>=2 for count in counts.values())>=2:informative+=1
        locus_summaries.append({
            "locus":locus,"original_reference":contigs[locus],
            "n_shared_Q20_majority_sites_all7":len(shared),
            "n_shared_variable_sites_all7":variable,
            "n_shared_parsimony_informative_sites_all7":informative,
            "n_species_Q20_majority_calls":N_DEEP_SAMPLES,
            "eligible_for_future_orthology_review":len(shared)>=100,
            "not_accepted_phylogenetic_alignment":True,
        })
    report={
       "schema":"azami_T2_seven_source120k_original110_Q20_shared_reference_sites_v1",
       "status":"EMPIRICAL_REAL_SEVEN_SOURCE_SAMPLES_NUCLEOTIDE_SITE_INTERSECTION_NOT_SPECIES_TREE",
       "source_S1_WCVP_master_taxon_keys":244,
       "source_S1_name_unit_master":251,
       "source_all113_nuclear_name_candidate_taxa_preserved":True,
       "sample_count":N_DEEP_SAMPLES,
       "same_run_120k_source_sequences":True,
       "original_low_warning_241_locus_panel_not_moreyra_final350":True,
       "source_original241_shared_copy_target_candidates":N_SHARED_TARGETS,
       "n_targets_all7_shared_Q20_majority_atleast100_sites":sum(z["eligible_for_future_orthology_review"] for z in locus_summaries),
       "n_original_loci_with_one_or_more_variable_Q20_position":sum(z["n_shared_variable_sites_all7"]>0 for z in locus_summaries),
       "n_original_loci_with_parsimony_informative_site":sum(z["n_shared_parsimony_informative_sites_all7"]>0 for z in locus_summaries),
       "total_shared_Q20_majority_positions_all7":sum(z["n_shared_Q20_majority_sites_all7"] for z in locus_summaries),
       "total_shared_variable_positions_all7":sum(z["n_shared_variable_sites_all7"] for z in locus_summaries),
       "total_shared_parsimony_informative_positions_all7":sum(z["n_shared_parsimony_informative_sites_all7"] for z in locus_summaries),
       "minimum_call_depth":MIN_SUPPORT,
       "minimum_base_quality":Q20,
       "minimum_mapping_quality":Q20,
       "required_major_base_frequency":MIN_MAJOR_FRACTION,
       "species_tree_inferred":False,
       "single_copy_orthology_or_homeology_validated":False,
       "unvalidated_reference_copy_and_voucher_warnings":True,
       "raw_reads_or_alignment_BAM_preserved":False,
    }
    return locus_summaries,report

def csv_write(path,rows):
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode",choices=["one","combine"],required=True)
    ap.add_argument("--run-id")
    ap.add_argument("--source",type=Path)
    ap.add_argument("--chosen",type=Path)
    ap.add_argument("--author-loci",type=Path)
    ap.add_argument("--targets",type=Path)
    ap.add_argument("--reference",type=Path)
    ap.add_argument("--scratch",type=Path)
    ap.add_argument("--reports",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    if a.mode=="one":
        if not all((a.run_id,a.source,a.chosen,a.targets,a.author_loci,a.reference,a.scratch)):
            raise ValueError("Real one-run mode requires frozen source and original ref inputs")
        allrows=original.load_source(a.source)
        chosen=read_csv(a.chosen)
        candidates=targets(read_csv(a.targets))
        sample=source_sample(allrows,chosen,a.run_id)
        original241=depth.loci_from_pinned_author_list(a.author_loci)
        result=fetch_align_and_call(sample,a.reference,candidates,original241,a.scratch)
        validate_single(result,sample,candidates)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
        print("T2_ORIGINAL110_REAL120K_Q20_SOURCE_SITE_CALLS="+json.dumps({
            "run":a.run_id,"n_callable_sites":result["n_Q20_callable_positions_total"]
        }),flush=True)
    else:
        if not a.reports or not a.targets:raise ValueError("Seven source receipts and original target locus list required")
        contigs=targets(read_csv(a.targets))
        paths=sorted(a.reports.rglob("SRR*.json"))
        if len(paths)!=N_DEEP_SAMPLES:
            raise ValueError("Must have seven independently acquired original run call receipts")
        src=[json.loads(p.read_text(encoding="utf-8")) for p in paths]
        for record in src:
            validate_single(record,record,contigs)
        loci,report=aggregate(src,contigs)
        a.output.mkdir(parents=True,exist_ok=True)
        csv_write(a.output/"ORIGINAL110_REAL120K_SEVEN_SPECIMENS_SHARED_Q20_SITES.csv",loci)
        (a.output/"ORIGINAL110_REAL120K_SEVEN_TAXON_Q20_SITE_SUMMARY.json").write_text(
            json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
        print("T2_SEVEN_NUCLEAR_REAL120K_SHARED_NUCLEOTIDE_SITES="+json.dumps(report,sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
