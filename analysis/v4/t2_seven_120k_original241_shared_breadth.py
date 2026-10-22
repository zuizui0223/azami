"""Re-measure original source 241-locus Q20 breadth in 7 same-run 120k ENA samples.

This independently recovers the *per-locus* 120k statistics omitted from the
full-125 numerical depth census, with the seven specimens selected by the
already archived 120k summary (>=100 loci, >=3x, >=50% reference length).

Candidate targets are for phylogenomic assembly QC only: the best one of
2–3 Asteraceae reference copies is not an orthology certificate; coverage of
one target does not prove that the same homologous sites/alleles can be aligned
across taxa. Do not identify any trait transition or infer rain adaptation.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path
import shutil

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as all125
from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as large
from analysis.v4 import t2_moreyra_all125_author241_site_depth as depth

N_SOURCE=125
N_SELECTED=7
N_LOCI=241
READ_PAIRS=120000
TAXA={
 "Cirsium_vulgare","Cirsium_tracyi","Cirsium_pulcherrimum",
 "Cirsium_italicum","Cirsium_flodmanii","Cirsium_latifolium",
 "Cirsium_scopulorum"
}
MINIMUM_DEPTH3_BREADTH=.25

def load_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as fp:
        rows=list(csv.DictReader(fp))
    if not rows:raise ValueError("Missing frozen public source run")
    return rows

def normalize(x):
    return "_".join(str(x).replace(" ","_").split("_")).strip("_")

def validate_sources(source,chosen,old120k):
    if len(source)!=N_SOURCE or len({z["run"] for z in source})!=N_SOURCE:
        raise ValueError("Original S1 125 fully identified Moreyra SRA runs missing")
    if len(chosen)!=N_SELECTED or {normalize(z["source_taxon"]) for z in chosen}!=TAXA:
        raise ValueError("Seven source-authenticated 120k high-quality focal taxa missing")
    allruns={r["run"]:r for r in source}
    pair={r["run"]:r for r in old120k}
    if len(pair)!=N_SOURCE or set(pair)!=set(allruns):
        raise ValueError("Original all125 same-run first120k Q20 receipts incomplete")
    if len({r["run"] for r in chosen})!=N_SELECTED or \
       len({r["accepted_key"] for r in chosen})!=N_SELECTED:
        raise ValueError("Seven must be distinct source runs and accepted keys")
    approved=[]
    for t in chosen:
        run=t["run"]
        match=allruns.get(run)
        earlier=pair.get(run)
        if match is None or earlier is None:
            raise ValueError("7-run pilot selected accession outside frozen 125 sources")
        for key in ("biosample","accepted_key"):
            if match[key]!=t[key] or match[key]!=earlier[key]:
                raise ValueError("Real SRA run/BioSample/WCVP accepted taxon source mismatch")
        if normalize(t["source_taxon"])!=normalize(match["S1_source_taxon"]):
            raise ValueError("Public run originally assigned a different source taxon")
        if earlier["source_30k_to_120k_valid"]!="True":
            raise ValueError("One of seven original 120k source FASTQ archives was unavailable")
        num=int(earlier["first120k_n_loci_dp3_ge50pct"])
        if int(t["same_run_120k_dp3_ge50pct"])!=num or num<100:
            raise ValueError("Seven-run source high-coverage selection no longer reproducible")
        approved.append((match,earlier))
    return approved

def analyse_one(source,original_241,ref,scratch,acquire=large.acquire_larger_prefix,
                mapper=large.map_read_pairs):
    folder=Path(scratch)/source["run"]
    folder.mkdir(parents=True,exist_ok=True)
    try:
        mates,_meta=acquire(source)
        if len(mates)!=2 or any(len(x)!=4*READ_PAIRS for x in mates):
            raise ValueError("Original ENA full source did not provide exact 120k paired reads")
        result=mapper(source,ref,original_241,mates,folder,include_per_locus=True)
        if result.get("input_read_pairs")!=READ_PAIRS:
            raise ValueError("Inaccurate original 120k read spot count")
        loci=result.get("per_locus_241")
        if not isinstance(loci,list) or len(loci)!=N_LOCI or \
            {x["locus"] for x in loci}!=set(original_241):
            raise ValueError("Missing source-authenticated author241 full per-locus site breadth")
        for locus in loci:
            if (locus["reference_choice_maximizes_coverage_not_orthology"] is not True
                    or not 0<=locus["dp5_breadth"]<=locus["dp3_breadth"]<=locus["dp1_breadth"]<=1):
                raise ValueError("Invalid measured one-reference quality breadth")
        return result
    finally:
        shutil.rmtree(folder,ignore_errors=True)

def combine(approved,original_241,source_results):
    if len(approved)!=N_SELECTED or len(source_results)!=N_SELECTED:
        raise ValueError("Seven original SRA runs must actually pass sequence depth scan")
    byrun={}
    for row in source_results:
        if row["run"] in byrun:
            raise ValueError("Public BioSample run duplicated in 7 source")
        byrun[row["run"]]=row
    per_locus=defaultdict(list)
    taxa=[]
    for source,earlier in approved:
        result=byrun.get(source["run"])
        if result is None or result.get("biosample")!=source["biosample"]:
            raise ValueError("Missing original 120k source run / specimen evidence")
        if result["read_pairs"]!=READ_PAIRS or len(result["loci"])!=N_LOCI:
            raise ValueError("Did not extract 241 original low-warning nuclear target loci")
        for key in ("n_loci_dp3_ge25pct","n_loci_dp3_ge50pct"):
            if result[key]!=int(earlier["first120k_"+key]):
                raise ValueError("Real same-run 120k depth differs from independent full125 screen")
        taxa.append({
            "source_name":source["S1_source_taxon"],
            "run":source["run"],"biosample":source["biosample"],
            "accepted_key":source["accepted_key"],"paired_reads":READ_PAIRS,
            "n_loci_dp3_ge25pct":result["n_loci_dp3_ge25pct"],
            "n_loci_dp3_ge50pct":result["n_loci_dp3_ge50pct"],
            "original_350_locus_tree_reproduced":False,
        })
        for entry in result["loci"]:
            per_locus[entry["locus"]].append(entry)
    if set(per_locus)!=set(original_241) or any(len(v)!=N_SELECTED for v in per_locus.values()):
        raise ValueError("Cross-taxon source target locus universes do not coincide")
    shared=[]
    for locus in original_241:
        records=per_locus[locus]
        dp25=sum(x["dp3_breadth"]>=.25 for x in records)
        dp50=sum(x["dp3_breadth"]>=.5 for x in records)
        copy_count=Counter(x["best_original_reference"] for x in records)
        one_ref=len(copy_count)==1
        shared.append({
           "locus":locus,
           "n_source_taxa":len(records),
           "n_Q20_dp3_ge25pct":dp25,
           "n_Q20_dp3_ge50pct":dp50,
           "best_reference_copy_identical_all7":one_ref,
           "modal_original_ref_n":max(copy_count.values()),
           "modal_original_ref":max(copy_count,key=lambda x:(copy_count[x],x)),
           "original_ref_names":"|".join(sorted(copy_count)),
           "deep7_all_Q20_dp3_ge25pct":dp25==N_SELECTED,
           "deep7_all_Q20_dp3_ge50pct":dp50==N_SELECTED,
           "deep7_all_Q20_dp3_ge25pct_same_ref":dp25==N_SELECTED and one_ref,
           "deep7_min_Q20_dp3_breadth":min(x["dp3_breadth"] for x in records),
           "nuclear_orthology_verified":False,
           "shared_nucleotide_site_positions_not_yet_ascertained":True,
        })
    summary={
       "schema":"azami_v4_deep7_real_120k_original_author241_Q20_joint_breadth_v1",
       "status":"ACTUAL_SOURCE_IDENTIFIED_MULTITAXON_120K_TARGET_SHARED_DEPTH_ONLY",
       "original_S1_accepted_keys_master":244,
       "original_S1_source_name_units_master":251,
       "original_source_125_run_registry_preserved":True,
       "source_run_count":N_SELECTED,
       "unique_WCVP_taxon_keys":N_SELECTED,
       "quality":"same read-acquisition, bwa mem, samtools depth baseQ20/mapQ20, 3x best single original reference copy",
       "original_241_no_warning_target_names":N_LOCI,
       "n_original_target_loci_seven_taxa_ge25pct_Q20":sum(x["deep7_all_Q20_dp3_ge25pct"] for x in shared),
       "n_original_target_loci_seven_taxa_ge50pct_Q20":sum(x["deep7_all_Q20_dp3_ge50pct"] for x in shared),
       "n_original_target_loci_seven_taxa_ge25pct_Q20_same_single_ref":sum(x["deep7_all_Q20_dp3_ge25pct_same_ref"] for x in shared),
       "n_original_target_loci_at_least6of7_ge50pct_Q20":sum(x["n_Q20_dp3_ge50pct"]>=6 for x in shared),
       "actual_allelic_genotypes_called":False,
       "high_coverage_does_not_prove_orthology":True,
       "same_reference_is_not_single_copy_paralog_audit":True,
       "nuclear_species_tree_constructed":False,
       "real_read_scratch_erased_and_not_uploaded":True,
       "source_seven_taxa":[x["source_name"] for x in taxa],
       "CAUTION":"Coverage records share reference locus identifiers; not complete homologous species sequences, and source taxon identity still needs specimen/voucher acceptance audit.",
    }
    return taxa,shared,summary

def run(source_path,chosen_path,pair_path,locus_path,ref,scratch,output,*,runner=analyse_one):
    source=all125.load_source(source_path)
    selected=load_csv(chosen_path)
    pairs=load_csv(pair_path)
    approved=validate_sources(source,selected,pairs)
    loci=depth.loci_from_pinned_author_list(locus_path)
    if len(loci)!=N_LOCI or not Path(ref).is_file():
        raise ValueError("Original locus names or Compositae DNA reference missing")
    collected=[]
    for s,_p in approved:
        print("T2_REAL_SEVEN120K_SOURCE_START="+s["run"],flush=True)
        result=runner(s,loci,ref,scratch)
        subset={
            "run":s["run"],"biosample":s["biosample"],
            "read_pairs":result["input_read_pairs"],
            "n_loci_dp3_ge25pct":result["n_loci_dp3_ge25pct"],
            "n_loci_dp3_ge50pct":result["n_loci_dp3_ge50pct"],
            "loci":result["per_locus_241"],
        }
        collected.append(subset)
        print("T2_REAL_SEVEN120K_SOURCE_DONE="+s["run"],flush=True)
    taxa,shared,report=combine(approved,loci,collected)
    Path(output).mkdir(parents=True,exist_ok=True)
    for name,rows in (("ORIGINAL241_SEVEN_SOURCE_120K_Q20_COMMON_TARGETS.csv",shared),
                      ("SEVEN_AUTHENTIC_PUBLIC_SRA_120K_SAMPLES.csv",taxa)):
        with (Path(output)/name).open("w",newline="",encoding="utf-8") as h:
            wr=csv.DictWriter(h,fieldnames=list(rows[0]),extrasaction="raise")
            wr.writeheader();wr.writerows(rows)
    (Path(output)/"SEVEN_REAL_120K_AUTHOR241_SHARED_LOCUS_SUMMARY.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
        encoding="utf-8")
    print("T2_REAL_SEVEN120K_AUTHOR241_SHARED_LOCUS="+json.dumps(report,ensure_ascii=False,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("source","chosen","pair120k","loci","reference","scratch","out-dir"):
        p.add_argument("--"+arg,required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.chosen,a.pair120k,a.loci,a.reference,a.scratch,a.out_dir)
if __name__=="__main__":
    main()
