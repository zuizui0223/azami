"""Post-result 110-taxon technical sensitivity, NOT an a priori 111-tip result.

The source-authenticated original111 Q20 120k-read comparison succeeded for
all 111 biological samples after exactly three ENA transport retries, but the
a priori >=1,000 observed positions per taxon gate failed because one specimen,
Cirsium tioganum WCVP206740970, had only 119/29310 masked Q20 positions.
This module DOES NOT lower the original gate or modify original 108 valid
receipts. It retains the 111-taxon same-reference, >=89/111 callable positions
and 108/118 original targets exactly; only removes one explicitly failed tip
for a separately flagged post-hoc technical ML source-sensitivity screen.

Q20 majority position calls are from real public source reads, not validated
single-copy orthologs or a dated, accepted genus-wide species tree. Exclude no
phenotype on the basis of evolutionary direction; both author-scored nodding
anchors remain. Do not infer rainfall-driven trait gains or adaptive radiation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from analysis.v4 import t2_real111_original118_Q20_nucleotide_calls as source
from analysis.v4 import t2_real111_original118_shared_Q20_aggregate as aggregate

UNMODIFIED_Q20_SOURCE_RUN=39889925464
MASTER_WCVP=244
MASTER_NAME_UNITS=251
BEFORE_SAMPLES=111
AFTER_SAMPLES=110
PREDECLARED_111_QUORUM=89
ORIGINAL_LOCI=118
RECOVERED_LOCI=108
ORIGINAL_SHARED_SITES=29310
ORIGINAL_PARSIMONY_SITES=1438
MIN_TIP_Q20_CALLABLE=1000
MIN_INFORMATIVE_Q20_SITES=50
EXCLUDED_TIP="WCVP206740970"
EXCLUDED_NAME="Cirsium tioganum"
EXCLUDED_ORIGINAL_Q20_BASES=119
REQUIRED_NODDING={"Cirsium conspicuum","Cirsium scopulorum"}
NA63_EXACT_FOUR={"Cirsium flodmanii","Cirsium scopulorum",
                 "Cirsium tracyi","Cirsium vulgare"}


def source_prehistory(shards,samples,targets):
    reports=aggregate.load_shards(shards)
    rows=source.read_csv(samples)
    loci=source.read_csv(targets)
    data,refs=aggregate.collect(rows,loci,reports)
    if len(data)!=BEFORE_SAMPLES or len(refs)!=ORIGINAL_LOCI:
        raise ValueError("Original full111 source Q20 sample/locus inventory changed")
    if len({r["source_original_run"] for r in data})!=BEFORE_SAMPLES:
        raise ValueError("Original independent source SRA runs duplicated")
    if {r["tip"] for r in data}!={r["phylogenetic_tip"] for r in rows}:
        raise ValueError("Original predeclared111 nuclear Q20 tip IDs changed")
    genes,aligned,breadth,report=aggregate.align_usable_sites(data,refs)
    if (len(genes)!=ORIGINAL_LOCI or
        report["n_loci_with_atleast100_site_positions_in_80pct_taxa"]!=RECOVERED_LOCI or
        report["n_full111_alignable_positions_called_by_80pct_taxa"]!=ORIGINAL_SHARED_SITES or
        report["n_parsimony_informative_positions_in_retained_loci"]!=ORIGINAL_PARSIMONY_SITES or
        report["technical_ML_alignment_source_gate_passed"]):
        raise ValueError("Predeclared real111 gate did not fail in its original source-defined way")
    return rows,data,genes,aligned,breadth,report


def sensitivity(rows,data,genes,aligned,breadth,report):
    source_tip={r["phylogenetic_tip"]:r for r in rows}
    source_details={r["tip"]:r for r in data}
    wrong=[(tip,x["n_Q20_bases_in_80pct_site_alignment"])
           for tip,x in sorted(breadth.items())
           if x["n_Q20_bases_in_80pct_site_alignment"]<MIN_TIP_Q20_CALLABLE]
    if wrong!=[(EXCLUDED_TIP,EXCLUDED_ORIGINAL_Q20_BASES)]:
        raise ValueError("Exactly one genuine original111 Q20 low-coverage specimen expected")
    if source_details[EXCLUDED_TIP]["taxon_name"]!=EXCLUDED_NAME:
        raise ValueError("Cannot mislabel low coverage reference-divergent original specimen")
    kept={tip:s for tip,s in aligned.items() if tip!=EXCLUDED_TIP}
    if len(kept)!=AFTER_SAMPLES or len(set(kept))!=AFTER_SAMPLES:
        raise ValueError("Post hoc sensitivity must keep exactly 110 distinct source BioSamples")
    if any(len(seq)!=ORIGINAL_SHARED_SITES or set(seq)-set("ACGTN") for seq in kept.values()):
        raise ValueError("Posthoc technical original reference Q20 nucleotide mask changed")
    if min(len(s)-s.count("N") for s in kept.values())<MIN_TIP_Q20_CALLABLE:
        raise ValueError("One of posthoc 110 original sources remains under original>=1000 Q20 cutoff")
    remaining={source_details[tip]["taxon_name"] for tip in kept}
    if not REQUIRED_NODDING<=remaining or not NA63_EXACT_FOUR<=remaining:
        raise ValueError("Independent original author nodding and NA63 quartet lost")
    # Do NOT refit locus selection to 110, which would be post hoc.
    chosen=[z for z in genes if z["retain_min100_Q20_positions"]]
    if len(chosen)!=RECOVERED_LOCI or sum(
        x["n_sites_same_reference_called_at_least80pct_taxa"] for x in chosen
    )!=ORIGINAL_SHARED_SITES:
        raise ValueError("Original111 108/29310 source locus/position mask changed")
    info=var=0
    counts=Counter()
    for j in range(ORIGINAL_SHARED_SITES):
        bases=Counter(seq[j] for seq in kept.values() if seq[j]!="N")
        k=sum(bases.values())
        if k<PREDECLARED_111_QUORUM-1:
            raise ValueError("An original source111 89-tip position lost more than one sample")
        counts[k]+=1
        info+=sum(v>=2 for v in bases.values())>=2
        var+=len(bases)>1
    if info<MIN_INFORMATIVE_Q20_SITES:
        raise ValueError("Posthoc 110 real original position alignment lacks independent SNP signal")
    sample_rows=[]
    for tip,seq in sorted(kept.items()):
        sourceinfo=source_details[tip]
        sample_rows.append({
            "source_tip":tip,
            "taxon_name":sourceinfo["taxon_name"],
            "original_run":sourceinfo["source_original_run"],
            "source_BioSample":sourceinfo["source_BioSample"],
            "n_actual_Q20_called":len(seq)-seq.count("N"),
            "fraction_observed":round((len(seq)-seq.count("N"))/len(seq),9),
            "source_read_count":sourceinfo["n_real_read_pairs"],
            "independently_verified_one_copy_orthology":False,
            "original_sample_selection_after_phenotype":False,
            "data_type":"same_111_predeclared_Q20_ref_coordinate_subset",
        })
    if len({row["source_BioSample"] for row in sample_rows})!=AFTER_SAMPLES:
        raise ValueError("Technical candidate duplicate source BioSample")
    boundaries=[];start=1
    for row in sorted(chosen,key=lambda x:x["locus"]):
        sites=row["n_sites_same_reference_called_at_least80pct_taxa"]
        boundaries.append({
            "original_locus":row["locus"],
            "same_Asteraceae_reference_copy":row["original_reference_copy"],
            "start_1based":start,"end_1based":start+sites-1,
            "original_quorum":PREDECLARED_111_QUORUM,
            "not_single_copy_Cirsium_orthology":True,
        })
        start+=sites
    summary={
        "schema":"azami_T2_original111_to_posthoc110_Q20_technical_sensitivity_v1",
        "status":"ACTUAL_PUBLIC_NUCLEAR_REFERENCE_COPY_Q20_POSTHOC_ONE_SPECIMEN_TECHNICAL_SENSITIVITY",
        "original_source_workflow":UNMODIFIED_Q20_SOURCE_RUN,
        "original_all111_Q20_successful_real_SRA_BioSamples":BEFORE_SAMPLES,
        "original111_technical_gate_passed":False,
        "source111_single_species_Q20_outlier":{"tip":EXCLUDED_TIP,
           "taxon":EXCLUDED_NAME,
           "n_real_Q20_positions":EXCLUDED_ORIGINAL_Q20_BASES,
           "all_original111_Q20_shared_positions":ORIGINAL_SHARED_SITES,
           "read_gene_absence_conclusion":False,
           "root_cause_is_unresolved_read_reference_compatibility":True},
        "posthoc_sensitivity_taxa":AFTER_SAMPLES,
        "original_S1_WCVP_master_accepted_keys":MASTER_WCVP,
        "original_S1_source_name_operational_units":MASTER_NAME_UNITS,
        "original_111_source_118_initial_loci":ORIGINAL_LOCI,
        "original_111_source_108_retained_loci":RECOVERED_LOCI,
        "original_111_predeclared_89_of111_quorum_unmodified":PREDECLARED_111_QUORUM,
        "original_111_reference_position_selection_recomputed_on110":False,
        "observed_Q20_reference_positions_reused_unchanged":ORIGINAL_SHARED_SITES,
        "original111_parsimony_informative_sites":ORIGINAL_PARSIMONY_SITES,
        "posthoc110_actual_parsimony_informative_sites":info,
        "posthoc110_variable_reference_positions":var,
        "posthoc110_min_Q20_calls_per_taxon":min(x["n_actual_Q20_called"] for x in sample_rows),
        "Q20_site_occupancy_after_one_specimen_omission_distribution":dict(sorted(counts.items())),
        "independently_author_coded_two_nodding_names_retained":sorted(REQUIRED_NODDING),
        "published_NA63_four_exact_name_overlap":sorted(NA63_EXACT_FOUR),
        "one_failed_original_taxon_still_tracked_as_missing_not_zero_genes":True,
        "posthoc_quality_based_tip_exclusion_not_a_preregistered111_tip_model":True,
        "source_original_BioSamples_not_pooled":True,
        "raw_FASTQ_BAM_in_artifact":False,
        "copy_number_paralogy_polyploidy_review_completed":False,
        "published_Moreyra_350_locus_tree_recreated":False,
        "global_EAzami294_true_species_tree_accepted":False,
        "true_Mk_nodding_transition_or_rainfall_adaptation_supported":False,
        "technical_GTRG4_ML_source_sensitivity_permitted":True,
    }
    return kept,sample_rows,boundaries,summary


def write_csv(path,rows):
    with Path(path).open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-shards",type=Path,required=True)
    p.add_argument("--repaired-shards",type=Path,required=True)
    p.add_argument("--source-samples",type=Path,required=True)
    p.add_argument("--source-targets",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    # Repaired shards are archived as derived distinct evidence, not overwrites.
    old=aggregate.load_shards(a.original_shards)
    new=aggregate.load_shards(a.repaired_shards)
    old_success={x["source_original_run"]:json.dumps(x,sort_keys=True)
                 for r in old for x in r["outcomes"] if x["status"]!="failed"}
    new_success={x["source_original_run"]:json.dumps(x,sort_keys=True)
                 for r in new for x in r["outcomes"] if x["status"]!="failed"}
    if len(old_success)!=108 or len(new_success)!=111 or any(
        new_success.get(run)!=s for run,s in old_success.items()
    ):
        raise ValueError("All 108 previously valid source Q20 genotype records must remain identical")
    samples,data,genes,seqs,breadth,report=source_prehistory(
        a.repaired_shards,a.source_samples,a.source_targets)
    filtered,rows,bounds,summary=sensitivity(samples,data,genes,seqs,breadth,report)
    summary["actual_repaired_shard_digest_sha256"]={
        x.name:hashlib.sha256(x.read_bytes()).hexdigest()
        for x in sorted(a.repaired_shards.rglob("shard_*.json"))}
    a.outdir.mkdir(parents=True,exist_ok=True)
    write_csv(a.outdir/"REAL110_SOURCE_TAXON_BIOSAMPLE_TECHNICAL_SENSITIVITY.csv",rows)
    write_csv(a.outdir/"REAL110_ORIGINAL111_LOCUS_108_Q20_PARTITIONS.csv",bounds)
    with (a.outdir/"REAL110_POSTHOC_ORIGINAL111_Q20_SITE_ALIGNMENT.fasta").open(
            "w",encoding="ascii") as handle:
        for tip,seq in sorted(filtered.items()):
            handle.write(f">{tip}\n")
            for j in range(0,len(seq),80):
                handle.write(seq[j:j+80]+"\n")
    summary["exact_observed_Q20_subset_FASTA_SHA256"]=hashlib.sha256(
        (a.outdir/"REAL110_POSTHOC_ORIGINAL111_Q20_SITE_ALIGNMENT.fasta").read_bytes()).hexdigest()
    (a.outdir/"REAL110_POSTHOC_ONLY_LOW_Q20_TIP_REMOVED_RECEIPT.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    print("T2_ORIGINAL111_TO_POSTHOC110_GENUINE_Q20_SENSITIVITY="+json.dumps({
        "original_source111":111,"posthoc_taxa":summary["posthoc_sensitivity_taxa"],
        "excluded":summary["source111_single_species_Q20_outlier"],
        "original_Q20_sites":summary["observed_Q20_reference_positions_reused_unchanged"],
        "posthoc_informative":summary["posthoc110_actual_parsimony_informative_sites"],
        "same_original_loci":summary["original_111_source_108_retained_loci"]
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
