"""Empirical 125-run -> 113-taxon Moreyra 241-locus 120k shared-target census.

After per-run high-quality depth measurement, use ONE original BioSample/run
per source taxon (max number of >=3x/>=50%-length source-241 loci, tie
>=25%, then lexicographic run). This is an optimistic QC selection for
phylogenetic DATA READINESS, not representative taxon evolution inference.
The 2 shorter original FASTQ run archives stay NULL / explicitly unsampled.

Compute taxon occupancy and SAME best-original-reference-copy occupancy for
25% and 50% site breadth. Original Asteraceae copy labels are compatibility
IDs and are NOT independently verified single-copy Cirsium orthologs.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path
import statistics

RUNS=125
TAXA=113
LOCI=241
NSHARDS=10
SHORT_RUNS={"SRR30887185","SRR30887248"}
PRIMARY_COVERAGE=.80
NODDING_SOURCE_NAMES={"Cirsium_conspicuum","Cirsium_scopulorum"}
ORIGINAL_SEVEN={"Cirsium_vulgare","Cirsium_tracyi","Cirsium_pulcherrimum",
                "Cirsium_italicum","Cirsium_flodmanii","Cirsium_latifolium",
                "Cirsium_scopulorum"}


def read_shards(root):
    filenames=sorted(Path(root).rglob("shard_*.json"))
    if len(filenames)!=NSHARDS:
        raise ValueError("Ten 120k full-run source file receipts required")
    reports=[json.loads(f.read_text(encoding="utf-8")) for f in filenames]
    if {r["shard"] for r in reports}!=set(range(NSHARDS)):
        raise ValueError("One original source 120k run shard missing or duplicated")
    if any(r["total_shards"]!=NSHARDS or r["raw_FASTQ_or_BAM_in_artifact"]
            for r in reports):
        raise ValueError("Source nuclear read artifact not allowed")
    return reports


def reconcile(source,reports):
    if len(source)!=RUNS or len({r["run"] for r in source})!=RUNS or \
       len({r["accepted_key"] for r in source})!=TAXA:
        raise ValueError("Complete original all125/113 S1 source census missing")
    byrun={r["run"]:r for r in source}
    seen={}
    for report in reports:
        if len(report["outcomes"])!=report["requested"]:
            raise ValueError("Source original run lost from shard")
        for r in report["outcomes"]:
            run=r["run"]
            if run not in byrun or run in seen:
                raise ValueError("Unexpected or repeated original public source run")
            src=byrun[run]
            if (r["accepted_key"],r["biosample"])!=(src["accepted_key"],src["biosample"]):
                raise ValueError("Original run/BioSample/WCVP accepted name mismatch")
            if r["status"] not in ("real_ENA_120k_same_run_author241_locus_Q20",
                                   "original_ENA_prefix_below_120k_pairs"):
                raise ValueError("Unresolved actual SRA run technical failure")
            if r["status"].startswith("real_"):
                if len(r["per_locus_241"])!=LOCI or r["original_120k_pairs"]!=120000:
                    raise ValueError("One full real sampled original241 locus matrix incomplete")
                if len({z["locus"] for z in r["per_locus_241"]})!=LOCI:
                    raise ValueError("Duplicate locus in original source run")
                for x in r["per_locus_241"]:
                    if not x["reference_choice_maximizes_coverage_not_orthology"]:
                        raise ValueError("Reference-copy quality scope incorrectly promoted")
                    if not (0<=x["dp5_breadth"]<=x["dp3_breadth"]<=x["dp1_breadth"]<=1):
                        raise ValueError("Biologically impossible Q20 site-breadth fractions")
            seen[run]=r
    if set(seen)!=set(byrun):
        raise ValueError("Missing source run from all 125 original FASTQ archives")
    missing={r["run"] for r in seen.values() if not r["status"].startswith("real_")}
    if missing!=SHORT_RUNS:
        raise ValueError("Previously confirmed two short original SRA runs drifted")
    observed=[r for r in seen.values() if r["status"].startswith("real_")]
    loci=set(observed[0]["per_locus_241"][i]["locus"] for i in range(LOCI))
    if any({z["locus"] for z in r["per_locus_241"]}!=loci for r in observed):
        raise ValueError("Different author241 locus universes across original SRA runs")
    return seen


def normalize(name):
    return " ".join(name.replace("_"," ").split())


def summarize(source,runs):
    if len(source)!=RUNS or len(runs)!=RUNS:
        raise ValueError("Cannot discard one original S1 nuclear accession")
    tax=defaultdict(list)
    for sr in source:
        tax[sr["accepted_key"]].append(sr)
    if len(tax)!=TAXA:
        raise ValueError("Original 113 WCVP source concepts missing")
    rows=[]
    current=[]
    for key, group in sorted(tax.items()):
        names={normalize(x["S1_source_taxon"]) for x in group}
        if len(names)!=1:
            raise ValueError("One accepted key carries multiple source names")
        srcname=next(iter(names))
        measured=[runs[x["run"]] for x in group if runs[x["run"]]["status"].startswith("real_")]
        if measured:
            # Repeated specimens are NOT joined at the sequence level.
            # Technical best-run selection is prior to across-species locus filter.
            chosen=sorted(measured,key=lambda x:(
                -x["n_loci_dp3_ge50pct"],-x["n_loci_dp3_ge25pct"],x["run"]))[0]
            observed=True
            current.append((key,srcname,chosen))
        else:
            chosen=None
            observed=False
        rows.append({
           "accepted_key":key,"original_source_taxon":srcname,
           "original_public_SRA_runs":len(group),
           "available_real120k_runs":len(measured),
           "selected_single_real_run":chosen["run"] if chosen else "",
           "selected_single_biosample":chosen["biosample"] if chosen else "",
           "best_single_run_Q20_3x_25pct_loci":chosen["n_loci_dp3_ge25pct"] if chosen else None,
           "best_single_run_Q20_3x_50pct_loci":chosen["n_loci_dp3_ge50pct"] if chosen else None,
           "real120k_depth_measured":observed,
           "other_original_runs_not_pooled":True,
           "source_name_not_voucher_verified":True,
        })
    if len(rows)!=TAXA or len(current)!=TAXA-len(SHORT_RUNS):
        raise ValueError("Expected 111 actual 120k public source name candidates")
    per_locus=defaultdict(list)
    for _key,_name,chosen in current:
        for z in chosen["per_locus_241"]:
            per_locus[z["locus"]].append(z)
    if len(per_locus)!=LOCI or any(len(v)!=len(current) for v in per_locus.values()):
        raise ValueError("Every one of 241 source loci must be measured in 111 chosen specimens")
    n=len(current)
    loci=[]
    for locus,records in sorted(per_locus.items()):
        qualifying25=[r for r in records if r["dp3_breadth"]>=.25]
        qualifying50=[r for r in records if r["dp3_breadth"]>=.50]
        best25=Counter(z["best_original_reference"] for z in qualifying25)
        best50=Counter(z["best_original_reference"] for z in qualifying50)
        count25=max(best25.values(),default=0)
        count50=max(best50.values(),default=0)
        fields={
          "locus":locus,
          "available_source_taxon_keys":n,
          "n_taxa_Q20_dp3_ge25pct":len(qualifying25),
          "n_taxa_Q20_dp3_ge50pct":len(qualifying50),
          "n_taxa_Q20_dp3_ge25pct_same_original_reference_copy":count25,
          "n_taxa_Q20_dp3_ge50pct_same_original_reference_copy":count50,
          "modal_Q20_dp3_ge25pct_original_ref":(
              sorted(best25,key=lambda x:(-best25[x],x))[0] if best25 else ""),
          "median_single_ref_Q20_dp3_breadth":
              statistics.median(z["dp3_breadth"] for z in records),
          "n_candidate_target_taxa_with_any_Q20_depth3":
              sum(z["dp3_callable_sites"]>0 for z in records),
          "not_orthologous_sequence_alignment":True,
        }
        loci.append(fields)
    scores={
      ">=50pct_taxa_Q20_dp3_ge25pct":sum(x["n_taxa_Q20_dp3_ge25pct"]>=math.ceil(n*.5) for x in loci),
      ">=65pct_taxa_Q20_dp3_ge25pct":sum(x["n_taxa_Q20_dp3_ge25pct"]>=math.ceil(n*.65) for x in loci),
      ">=80pct_taxa_Q20_dp3_ge25pct":sum(x["n_taxa_Q20_dp3_ge25pct"]>=math.ceil(n*.80) for x in loci),
      ">=90pct_taxa_Q20_dp3_ge25pct":sum(x["n_taxa_Q20_dp3_ge25pct"]>=math.ceil(n*.90) for x in loci),
      ">=80pct_taxa_Q20_dp3_ge25pct_same_ref_copy":sum(
          x["n_taxa_Q20_dp3_ge25pct_same_original_reference_copy"]>=math.ceil(n*.80) for x in loci),
      ">=80pct_taxa_Q20_dp3_ge50pct":sum(x["n_taxa_Q20_dp3_ge50pct"]>=math.ceil(n*.80) for x in loci),
      ">=80pct_taxa_Q20_dp3_ge50pct_same_ref_copy":sum(
          x["n_taxa_Q20_dp3_ge50pct_same_original_reference_copy"]>=math.ceil(n*.80) for x in loci),
    }
    name_map={normalize(name):val for _key,name,val in current}
    # Neither nodding state can be invented from photo vertical angle. The two
    # exact taxon names were author scored independently in Siniscalchi 2023.
    labels={normalize(x) for x in NODDING_SOURCE_NAMES}
    exact_both=all(name in name_map for name in labels)
    down_shared={k:0 for k in (
        "n_original241_two_author_nodding_names_both_Q20_dp3_ge25pct",
        "n_original241_two_author_nodding_names_both_Q20_dp3_ge25pct_same_ref")}
    if exact_both:
        nodding=[name_map[x] for x in sorted(labels)]
        for locus in per_locus:
            rec=[next(z for z in r["per_locus_241"] if z["locus"]==locus)
                 for r in nodding]
            if all(x["dp3_breadth"]>=.25 for x in rec):
                down_shared["n_original241_two_author_nodding_names_both_Q20_dp3_ge25pct"]+=1
                if len({x["best_original_reference"] for x in rec})==1:
                    down_shared["n_original241_two_author_nodding_names_both_Q20_dp3_ge25pct_same_ref"]+=1
    seven={normalize(x) for x in ORIGINAL_SEVEN}
    seven_present=seven<={normalize(n) for _k,n,_r in current}
    original7_shared25=None
    if seven_present:
        pool=[r for _k,name,r in current if normalize(name) in seven]
        original7_shared25=sum(all(
            next(z for z in r["per_locus_241"] if z["locus"]==locus)["dp3_breadth"]>=.25 for r in pool)
            for locus in per_locus)
    out={
      "schema":"azami_v4_real120k_S1_source125_113candidate_111measured_241loci_v1",
      "status":"EMPIRICAL_ALL_ELIGIBLE_SPECIES_SHARED_LOCUS_OCCUPANCY_NOT_SPECIES_TREE",
      "S1_original_WCVP_accepted_keys":244,
      "S1_original_source_name_units":251,
      "S1_original_eligible_image_observations":170066,
      "source_original_Moreyra_public_SRA_run_count":RUNS,
      "source_public_taxon_name_candidates":TAXA,
      "source_independently_measured_120k_taxa":n,
      "source_short_run_unavailable_not_zero_genes":len(SHORT_RUNS),
      "source_241_author_filtered_no_warning_loci":LOCI,
      "read_depth_quality":"samtools Q20 base and mapping; first120k authentic same-source ENA paired prefixes",
      "taxon_sample_choice":"max within-taxon number 3x/50% original241 loci; tie 25%, then exact run",
      "n_candidates_with_best_run_atleast100_3x50pct":sum(
           r["best_single_run_Q20_3x_50pct_loci"] is not None and
           r["best_single_run_Q20_3x_50pct_loci"]>=100 for r in rows),
      "cross_taxon_shared_locus_occupancy":scores,
      "two_original_NA_author_nodding_names_both_present":exact_both,
      **down_shared,
      "previous_seven_taxa_present_with_best_single_run_selection":seven_present,
      "n_loci_same_seven_name_subset_Q20_dp3_ge25pct":original7_shared25,
      "not_exact_published_350_loci":True,
      "same_original_reference_is_not_confirmed_orthology":True,
      "taxonomic_sample_name_is_not_voucher_identity":True,
      "no_294_tip_nuclear_tree_constructed":True,
      "no_orientation_evolution_or_BIO12_adaptation_inferred":True,
      "no_FASTQ_BAM_or_site_call_uploaded":True,
    }
    return rows,loci,out


def write_csv(path,rows):
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--receipts",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    source=[]
    with a.source.open(newline="",encoding="utf-8-sig") as f:source=list(csv.DictReader(f))
    reports=read_shards(a.receipts)
    outcomes=reconcile(source,reports)
    taxa,loci,summary=summarize(source,outcomes)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"T2_ALL113_SOURCE_TAXA_REAL120K_BEST_SINGLE_RUN.csv",taxa)
    write_csv(a.out_dir/"T2_AUTHOR241_SHARED_Q20_120K_ACROSS111_TAXA.csv",loci)
    (a.out_dir/"T2_REAL120K_AUTHOR241_ALL111_CROSS_TAXON_SUMMARY.json").write_text(
       json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_REAL120K_ALL113_SOURCE111_MEASURED_SHARED241="+json.dumps(summary,
            sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
