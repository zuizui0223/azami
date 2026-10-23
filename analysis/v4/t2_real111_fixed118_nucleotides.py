"""Call real identical-reference Q20 nucleotide positions across 111 Cirsium source SRA runs.

Each accession is strictly one EAzami294-vouched Moreyra public BioSample and
the already-evaluated same first-120000 ENA paired reads. Original reference
copies and 118 predeclared 80%-occupancy loci come from separate source-verified
depth studies. No absent calls are filled from the reference; heterozygous and
poorly-supported calls are absent, not alternative phenotype states.

This is NOT a confirmed single-copy ortholog, accepted Cirsium species tree,
dated phylogeny, or rainfall-driven repeated head-nodding conclusion.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import time

from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as reads
from analysis.v4 import t2_moreyra_all125_ena_read_mapping as source
from analysis.v4.t2_seven_real120k_shared_sites import call_sites_from_pysam

SOURCE_RUNS=125
SAMPLE_COUNT=111
TARGET_COUNT=118
READ_PAIRS=120000
SHARDS=10
AMBIGUOUS_ORIGINAL_SOURCE_RUNS={"SRR30887230"}
SOURCE_TAXA_FIELDS={"accepted_key","taxon_name","phylogenetic_tip","run","biosample"}
TARGET_FIELDS={"locus","original_compositae_reference","empirically_supported_source_taxa"}
RUN_ORDER=("run","biosample","accepted_key")
KEEP_Q20_MAJORITY=True


def load_csv(path,required):
    with Path(path).open(newline="",encoding="utf-8-sig") as handle:
        rows=list(csv.DictReader(handle))
    if not rows or not required<=set(rows[0]):
        raise ValueError("Expected pinned source CSV missing fields")
    return rows


def audit_samples(source125,chosen111,copy118):
    if len(source125)!=SOURCE_RUNS or len({x["run"] for x in source125})!=SOURCE_RUNS:
        raise ValueError("The original125 S1 public Moreyra run universe changed")
    if len(chosen111)!=SAMPLE_COUNT or len({x["run"] for x in chosen111})!=SAMPLE_COUNT:
        raise ValueError("The 111 real 120k selected source accessions changed")
    if len(copy118)!=TARGET_COUNT or len({x["locus"] for x in copy118})!=TARGET_COUNT:
        raise ValueError("The 118 predeclared original reference-copy loci changed")
    by_run={x["run"]:x for x in source125}
    selected={}
    for row in chosen111:
        run=row["run"]
        if run not in by_run:
            raise ValueError("Selected source run has no original EAzami/S1 crosswalk")
        original=by_run[run]
        if any(original[k]!=row[k] for k in RUN_ORDER):
            raise ValueError("Original public BioSample/run/accepted-key identity changed")
        if " ".join(row["taxon_name"].replace("_"," ").split()) != \
           " ".join(original["S1_source_taxon"].replace("_"," ").split()):
            raise ValueError("Original S1 source taxon name drift")
        if row["phylogenetic_tip"]!="WCVP"+row["accepted_key"] or not row["accepted_key"].isdigit():
            raise ValueError("One accepted WCVP key not preserved in topology label")
        selected[run]=row
    copies={}
    for r in copy118:
        locus,ref=r["locus"],r["original_compositae_reference"]
        if (not ref.startswith(("lett-","saff-","sunf-")) or
            not ref.endswith("-"+locus) or
            int(r["empirically_supported_source_taxa"])<89):
            raise ValueError("Failed preregistered all111 >=80%-copy occupancy gate")
        if ref in copies.values():
            raise ValueError("Duplicated original same-reference nuclear locus")
        copies[locus]=ref
    return by_run,selected,copies


def map_and_call(original,selected,copies,reference,scratch,*,acquire=reads.acquire_larger_prefix,
                 site_caller=call_sites_from_pysam):
    run=selected["run"]
    outdir=Path(scratch)/run
    outdir.mkdir(parents=True,exist_ok=True)
    try:
        mates,_ena=acquire(original)
        if len(mates)!=2 or any(len(m)!=4*READ_PAIRS for m in mates):
            raise ValueError("Original ENA public read sample not exactly 120000 pairs")
        fqs=[]
        for k,lines in enumerate(mates,1):
            p=outdir/f"sample_{k}.fastq"
            p.write_bytes(b"".join(lines))
            fqs.append(p)
        bam=outdir/"real_source.sorted.bam"
        bwa=subprocess.Popen(["bwa","mem","-t","2",str(reference),*(str(f) for f in fqs)],
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            sorter=subprocess.Popen(["samtools","sort","-@","2","-m","128M",
                                     "-o",str(bam),"-"],
                                    stdin=bwa.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            bwa.stdout.close()
            _,sort_err=sorter.communicate(timeout=300)
            _,bwa_err=bwa.communicate(timeout=300)
            if sorter.returncode or bwa.returncode:
                raise RuntimeError("Authentic reads BWA/samtools mapping failed: "+
                                   (sort_err+bwa_err).decode("utf-8","replace")[-350:])
        finally:
            if bwa.poll() is None:
                bwa.kill()
                bwa.wait()
        source.run_commands(["samtools","index",str(bam)])
        calls=site_caller(bam,copies)
        if set(calls)!=set(copies):
            raise ValueError("Missing predeclared original locus in nucleotide calls")
        total=sum(len(x["q20_majority_bases_1based"]) for x in calls.values())
        return {
            "status":"real_source111_120k_Q20_identical_copy_base_calls",
            "run":run,"biosample":selected["biosample"],
            "accepted_key":selected["accepted_key"],
            "taxon_name":selected["taxon_name"],
            "phylogenetic_tip":selected["phylogenetic_tip"],
            "original_tree_code":original["original_tree_code"],
            "source_run_label_conflict_unresolved":run in AMBIGUOUS_ORIGINAL_SOURCE_RUNS,
            "source_public_ENA_paired_read_count":READ_PAIRS,
            "n_predeclared_80percent_reference_loci":TARGET_COUNT,
            "n_real_Q20_majority_callable_positions":total,
            "observed_original_reference_q20_site_calls":calls,
            "uncalled_bases_not_filled_by_reference":True,
            "heterozygous_or_low_quality_sites_masked":True,
            "original_single_copy_orthology_verified":False,
            "biological_nodding_from_EXIF_inferred":False,
            "full_genus_nuclear_species_tree_estimated":False,
            "FASTQ_BAM_exported":False,
        }
    finally:
        shutil.rmtree(outdir,ignore_errors=True)


def validate_one(item,selected,copies):
    if any(item.get(k)!=selected[k] for k in ("run","biosample","accepted_key","taxon_name","phylogenetic_tip")):
        raise ValueError("Original source physical sample and Q20 site receipt differ")
    if item["source_public_ENA_paired_read_count"]!=READ_PAIRS or \
       item["n_predeclared_80percent_reference_loci"]!=TARGET_COUNT:
        raise ValueError("Incorrect original120k or 118-locus source scan")
    if set(item["observed_original_reference_q20_site_calls"])!=set(copies):
        raise ValueError("One or more source-authentic nuclear loci missing")
    bases_total=0
    for locus,detail in item["observed_original_reference_q20_site_calls"].items():
        if detail["reference"]!=copies[locus] or not detail["not_a_validated_allelic_genotype"]:
            raise ValueError("Another reference copy entered the source sample")
        last=0
        for pos,base,depth,frequency in detail["q20_majority_bases_1based"]:
            if (not isinstance(pos,int) or pos<=last or pos>int(detail["reference_length_nt"])
                    or base not in "ACGT" or depth<3 or not .9<=frequency<=1):
                raise ValueError("Invalid positionwise source DNA call or non-Q20 majority")
            last=pos
            bases_total+=1
    if bases_total!=item["n_real_Q20_majority_callable_positions"]:
        raise ValueError("Real DNA site count disagrees with exact original calls")
    if not item["uncalled_bases_not_filled_by_reference"] or \
            item["original_single_copy_orthology_verified"] or item["FASTQ_BAM_exported"]:
        raise ValueError("Biological source provenance and missingness claims invalid")
    return True


def shard(source125,chosen111,copies,reference,scratch,n,n_shards=SHARDS,runner=map_and_call):
    by_run,sel,copies=audit_samples(source125,chosen111,copies)
    if not(0<=n<n_shards==SHARDS):
        raise ValueError("Required original 10 sample batches")
    group=[x for i,x in enumerate(sorted(sel.values(),key=lambda r:r["run"]))
           if i%n_shards==n]
    records=[]
    for row in group:
        t=time.monotonic()
        try:
            item=runner(by_run[row["run"]],row,copies,reference,scratch)
            validate_one(item,row,copies)
            item["elapsed_seconds"]=round(time.monotonic()-t,2)
            records.append(item)
            print("REAL111_Q20_BASE_CALL_PASS="+json.dumps({
                "run":row["run"],"called_sites":item["n_real_Q20_majority_callable_positions"]}),
                flush=True)
        except Exception as e:
            records.append({"status":"failed","run":row["run"],
                 "biosample":row["biosample"],"accepted_key":row["accepted_key"],
                 "error_type":type(e).__name__,"error":str(e)[:350]})
            print("REAL111_Q20_BASE_CALL_FAIL="+json.dumps(records[-1]),flush=True)
    if len(records)!=len(group) or {r["run"] for r in records}!={r["run"] for r in group}:
        raise ValueError("Missing or duplicated original SRA accession in nucleotides")
    return {
        "schema":"azami_v4_real111_BioSamples_118_reference_Q20_site_calls_v1",
        "input_original_S1_master_taxa":244,"input_source_names":251,
        "original_public_moreyra_candidate_taxa":113,"selected_real120k_samples":111,
        "original_source_same_copy_loci":118,
        "shard":n,"n_shards":n_shards,
        "requested":len(group),"passed":sum(x["status"].startswith("real_") for x in records),
        "failed":sum(x["status"]=="failed" for x in records),
        "no_original_FASTQ_or_BAM_preserved":True,
        "not_phylogenetic_species_tree":True,
        "outcomes":records,
    }


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--source",type=Path,required=True)
    a.add_argument("--samples",type=Path,required=True)
    a.add_argument("--targets",type=Path,required=True)
    a.add_argument("--reference",type=Path,required=True)
    a.add_argument("--scratch",type=Path,required=True)
    a.add_argument("--shard",type=int,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    source125=source.load_source(args.source)
    selected=load_csv(args.samples,SOURCE_TAXA_FIELDS)
    chosen=load_csv(args.targets,TARGET_FIELDS)
    if not args.reference.is_file():
        raise ValueError("Original public source Compositae1061 reference not recovered")
    data=shard(source125,selected,chosen,args.reference,args.scratch,args.shard)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n",
                        encoding="utf-8")
    print("T2_REAL111_ALL118_SITE_CALL_SHARD="+json.dumps({
        "shard":data["shard"],"passed":data["passed"],"failed":data["failed"],
        "original_n":data["requested"]}),flush=True)


if __name__=="__main__":
    main()
