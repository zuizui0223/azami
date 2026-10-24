"""Actual shared-Q20 nucleotide coordinate census across 111 real public source samples.

Accept ONLY 10 complete original120k first-prefix SRA read shard receipts and
118 fixed EAzami original Compositae1061 reference copies. A site is counted
as shareable when >=89/111 independently matched BioSamples have Q20 majority
A/C/G/T at the *same original reference coordinate*. Unobserved or ambiguous
positions are N, never copied from reference.

This detects technical nuclear data connectivity; the 118 regions are not
proved one-copy orthologs; tree, BM/OU, trait-history and rain causality require
independent admission and may NOT be inferred by this summary.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path

N_RUNS=111
N_LOCI=118
N_SHARDS=10
Q20_QUORUM=89
MIN_SITES_PER_LOCUS=100
MIN_LOCUS_FOR_PRELIMINARY_ML=30
MIN_TOTAL_SITES_FOR_PRELIMINARY_ML=10000
MIN_PARSIMONY_INFORMATIVE_SITES=50
MIN_CALLABLE_SITES_PER_TAXON=1000
BOTANICAL_NODDING={"Cirsium conspicuum","Cirsium scopulorum"}


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        records=list(csv.DictReader(f))
    if not records:
        raise ValueError("Missing actual S1 source molecular evidence")
    return records


def load_shards(folder):
    paths=sorted(Path(folder).rglob("shard_*.json"))
    if len(paths)!=N_SHARDS:
        raise ValueError("Ten complete original Moreyra 111 source read shards required")
    reports=[json.loads(p.read_text(encoding="utf-8")) for p in paths]
    if {r["shard"] for r in reports}!=set(range(N_SHARDS)):
        raise ValueError("Duplicate or missing empirical nuclear read shard ID")
    if any(r["n_shards"]!=N_SHARDS or r["raw_reads_or_BAM_publication"] or
           not r["single_original_reference_not_confirmed_ortholog"]
           for r in reports):
        raise ValueError("One original SRA shard violated real 120k/orthology claim boundary")
    return reports


def collect(sample_rows,target_rows,reports):
    if (len(sample_rows)!=N_RUNS or len(target_rows)!=N_LOCI or
        len({r["run"] for r in sample_rows})!=N_RUNS or
        len({r["locus"] for r in target_rows})!=N_LOCI):
        raise ValueError("Source 111 biological BioSamples or original118 loci missing")
    samples={r["run"]:r for r in sample_rows}
    refs={r["locus"]:r["original_compositae_reference"] for r in target_rows}
    results={}
    for shard in reports:
        if shard["requested"]!=len(shard["outcomes"]):
            raise ValueError("One run disappeared within its original shard")
        if shard["failure"]:
            raise ValueError("Failed source-run Q20 DNA base calls cannot be relabeled as gene absence")
        for entry in shard["outcomes"]:
            run=entry["source_original_run"]
            if run not in samples or run in results:
                raise ValueError("Unknown or duplicate authentic real run in Q20 base-calling")
            s=samples[run]
            if (entry["status"]!="real_Q20_majority_reference_site_calls" or
                entry["source_BioSample"]!=s["biosample"] or
                entry["source_S1_accepted_key"]!=s["accepted_key"] or
                entry["taxon_name"]!=s["taxon_name"] or
                entry["tip"]!=s["phylogenetic_tip"] or
                entry["n_real_read_pairs"]!=120000 or
                entry["n_original_common_reference_loci"]!=N_LOCI or
                set(entry["loci"])!=set(refs)):
                raise ValueError("Original 111 source Q20 sample/locus identity changed")
            for locus,val in entry["loci"].items():
                if val["reference"]!=refs[locus] or not val["not_a_validated_allelic_genotype"]:
                    raise ValueError("Q20 base calls mapped against a different source reference copy")
            results[run]=entry
    if len(results)!=N_RUNS or len({x["tip"] for x in results.values()})!=N_RUNS:
        raise ValueError("Not all original111 independent biological nuclear source tips retrieved")
    names={x["taxon_name"] for x in results.values()}
    if not BOTANICAL_NODDING<=names:
        raise ValueError("Original two independently botanically nodding source taxa missing")
    return sorted(results.values(),key=lambda x:x["tip"]),refs


def align_usable_sites(data,reference_names,quorum=Q20_QUORUM):
    if not 1<=quorum<=len(data):
        raise ValueError("Invalid a priori biological taxon site occupancy gate")
    taxa=[x["tip"] for x in data]
    if len(taxa)!=len(set(taxa)):
        raise ValueError("One public WCVP original tip entered twice")
    taxon_seq={x["tip"]:[] for x in data}
    gene_rows=[]
    diagnostic_no_seq=[]
    for locus in sorted(reference_names):
        ref=reference_names[locus]
        bypos=defaultdict(dict)
        target_length=None
        for s in data:
            c=s["loci"][locus]
            if c["reference"]!=ref:
                raise ValueError("Locus reference-copy identity changed")
            length=c["reference_length_nt"]
            if target_length is None:target_length=length
            if length!=target_length:
                raise ValueError("Different source reference lengths for identical contig")
            seen=set()
            for p,base,dp,major in c["q20_majority_bases_1based"]:
                p=int(p);dp=int(dp);major=float(major)
                if (p in seen or not 1<=p<=length or base not in {"A","C","G","T"} or
                    dp<3 or not .9<=major<=1):
                    raise ValueError("Invalid original high-confidence Q20 nucleotide coordinate")
                seen.add(p)
                bypos[p][s["tip"]]=base
        selected=sorted(p for p,orig in bypos.items() if len(orig)>=quorum)
        variable=informative=0
        for p in selected:
            observed=Counter(bypos[p].values())
            variable+=len(observed)>1
            informative+=sum(v>=2 for v in observed.values())>=2
        keep=len(selected)>=MIN_SITES_PER_LOCUS
        gene_rows.append({
            "locus":locus,"original_reference_copy":ref,
            "original_reference_length":target_length,
            "n_sites_callable_at_least1_taxon":len(bypos),
            "n_sites_same_reference_called_at_least80pct_taxa":len(selected),
            "n_variable_sites_in_80pct_called":variable,
            "n_parsimony_informative_sites_in_80pct_called":informative,
            "retain_min100_Q20_positions":keep,
            "reference_copy_single_ortholog_verified":False,
        })
        if keep:
            for tip in taxa:
                taxon_seq[tip].append("".join(bypos[p].get(tip,"N") for p in selected))
        else:
            diagnostic_no_seq.append(locus)
    ordered={tip:"".join(taxon_seq[tip]) for tip in taxa}
    lengths={len(x) for x in ordered.values()}
    if len(lengths)!=1:
        raise ValueError("Different effective reference-anchored alignment coordinate lengths")
    effective=next(iter(lengths))
    n_loci=sum(x["retain_min100_Q20_positions"] for x in gene_rows)
    informative=sum(x["n_parsimony_informative_sites_in_80pct_called"]
                    for x in gene_rows if x["retain_min100_Q20_positions"])
    taxon_breadth={
        tip:{"n_Q20_bases_in_80pct_site_alignment":len(seq)-seq.count("N"),
             "total_80pct_site_alignment_positions":len(seq),
             "Q20_site_call_fraction":(len(seq)-seq.count("N"))/len(seq) if seq else 0}
        for tip,seq in ordered.items()
    }
    preliminary=(n_loci>=MIN_LOCUS_FOR_PRELIMINARY_ML and
                 effective>=MIN_TOTAL_SITES_FOR_PRELIMINARY_ML and
                 informative>=MIN_PARSIMONY_INFORMATIVE_SITES and
                 min(x["n_Q20_bases_in_80pct_site_alignment"]
                     for x in taxon_breadth.values())>=MIN_CALLABLE_SITES_PER_TAXON)
    return gene_rows,ordered,taxon_breadth,{
        "all_original_111_S1_biological_samples":len(data),
        "n_original_source_candidate_loci":len(reference_names),
        "predeclared_site_occupancy_taxa":quorum,
        "n_loci_with_atleast100_site_positions_in_80pct_taxa":n_loci,
        "n_full111_alignable_positions_called_by_80pct_taxa":effective,
        "n_parsimony_informative_positions_in_retained_loci":informative,
        "technical_ML_alignment_source_gate_passed":preliminary,
        "strict_single_copy_orthology_gate_passed":False,
        "orientation_evolution_rainfall_test_permitted":False,
        "source_bioSample_duplicates_collapsed":True,
        "n_low_coverage_loci_preserved_in_audit":len(diagnostic_no_seq),
        "minimum_real_data_ML_site_per_locus":MIN_SITES_PER_LOCUS,
        "minimum_source_usable_taxa_at_locus":quorum,
        "minimum_total_Q20_sites_for_preliminary_ML":MIN_TOTAL_SITES_FOR_PRELIMINARY_ML,
        "minimum_parsimony_informative_sites":MIN_PARSIMONY_INFORMATIVE_SITES,
        "minimum_Q20_observed_sites_per_taxon":MIN_CALLABLE_SITES_PER_TAXON,
        "missing_values_are_N_not_original_reference_base":True,
        "not_accepted_genus_or_global_species_tree":True,
    }


def export_csv(path,records):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),extrasaction="raise")
        w.writeheader();w.writerows(records)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples",required=True,type=Path)
    ap.add_argument("--targets",required=True,type=Path)
    ap.add_argument("--shards",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    args=ap.parse_args()
    samples=read_csv(args.samples)
    targets=read_csv(args.targets)
    reports=load_shards(args.shards)
    data,refs=collect(samples,targets,reports)
    genes,sequences,breadths,summary=align_usable_sites(data,refs)
    args.outdir.mkdir(parents=True,exist_ok=True)
    export_csv(args.outdir/"REAL111_ORIGINAL118_NUCLEAR_SAME_COPY_80PCT_SITE_CENSUS.csv",genes)
    export_csv(args.outdir/"REAL111_Q20_80PCT_ALIGNMENT_TAXON_MISSINGNESS.csv",[
       {"tip":tip,"taxon":next(x["taxon_name"] for x in data if x["tip"]==tip),
        **stat} for tip,stat in sorted(breadths.items())])
    if summary["technical_ML_alignment_source_gate_passed"]:
        with (args.outdir/"REAL111_80PCT_SOURCE_118_ORIGINAL_REF_Q20_ALIGNMENT.fasta").open(
                "w",encoding="ascii") as f:
            for tip,seq in sorted(sequences.items()):
                f.write(f">{tip}\n")
                for i in range(0,len(seq),80):
                    f.write(seq[i:i+80]+"\n")
        summary["real_source_fasta_technical_candidate_written"]=True
    else:
        summary["real_source_fasta_technical_candidate_written"]=False
    summary["schema"]="azami_T2_original111_real120k_118loci_80pct_same_base_position_v1"
    summary["status"]="EMPIRICAL_REAL_PUBLIC_111_SOURCE_NUCLEAR_SITE_MATRIX_Q20_NOT_VERIFIED_ORTHOLOGS"
    summary["original_S1_WCVP_keys_retained"]=244
    summary["original_S1_source_name_operational_units_retained"]=251
    summary["source_Moreyra_original125_public_runs"]=125
    summary["source_short_120k_runs_missing_not_genetic_absence"]=2
    summary["two_original_author_independent_nodding_taxa_in_source"]=sorted(BOTANICAL_NODDING)
    summary["FASTQ_or_BAM_archived"]=False
    (args.outdir/"REAL111_SHARED_Q20_NUCLEOTIDE_SITE_SUMMARY.json").write_text(
       json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
       encoding="utf-8")
    print("T2_REAL111_ACTUAL_SHARED_Q20_SITE_CENSUS="+json.dumps(summary,
           ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
