"""Audit-preserving repair: stitch one re-obtained real SRA into 110 successful DNA site records.

Source GitHub Actions #39888260140 ran 111 source-qualified SRA accessions:
110 authentic Q20 site calls were successful, and ONLY SRR30887180 failed
because the 64MiB HTTP gzip byte prefix was too short (not gene absence).
The exact same first120k public reads are reacquired with bounded 128MiB
compressed source transfer, retaining all former 110 outcomes unchanged.

This materializes a NEW separately sourced ten-shard receipt, not rewriting
the original archived failed run and not manufacturing a completed species tree.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from analysis.v4 import t2_real111_fixed118_nucleotides as nuclear

SHARDS=10
RUN="SRR30887180"
ORIGINAL_BIO="SAMN44017849"
ORIGINAL_SOURCE_ACCEPTED="206743475"
ORIGINAL_FAILED_SHARD=3
EXPECTED_FAILURE="Official ENA bounded 64MiB prefix cannot yield 120k complete reads"
TOTAL_TAXA=111
SOURCE_SUCCESS=110


def load_stage_files(folder):
    files=sorted(Path(folder).rglob("shard_*.json"))
    if len(files)!=SHARDS:
        raise ValueError("Must retain 10 original 111-source DNA sample receipts")
    reports=[]
    for f in files:
        raw=f.read_bytes()
        parsed=json.loads(raw)
        reports.append((parsed,{"filename":str(f.name),"sha256":hashlib.sha256(raw).hexdigest()}))
    if {x["shard"] for x,_ in reports}!=set(range(SHARDS)):
        raise ValueError("Original 10 site base-call shards duplicated/missing")
    return reports


def repair(original_stages,rescue,selected,refs):
    if rescue.get("run")!=RUN or rescue.get("biosample")!=ORIGINAL_BIO or \
       rescue.get("accepted_key")!=ORIGINAL_SOURCE_ACCEPTED or \
       rescue.get("rescue_protocol")!="same_real_first120000_pairs_larger_compressed_HTTP_prefix_only" or \
       rescue.get("rescue_max_compressed_bytes_per_mate")!=128*1024*1024:
        raise ValueError("Recovery is not the exact original source run/unchanged biological read experiment")
    wanted={x["run"]:x for x in selected}
    if len(selected)!=TOTAL_TAXA or RUN not in wanted:
        raise ValueError("Frozen 111 original S1 public BioSamples absent")
    copies={x["locus"]:x["original_compositae_reference"] for x in refs}
    nuclear.validate_one(rescue,wanted[RUN],copies)
    records=[]
    repaired=[]
    initial_pass=0
    original_fail=[]
    for stage,digest in original_stages:
        if stage["n_shards"]!=SHARDS or not stage["no_original_FASTQ_or_BAM_preserved"]:
            raise ValueError("Original shard provenance changed")
        if stage["requested"]!=len(stage["outcomes"]):
            raise ValueError("One original authentic source run lost before rescue")
        updated=deepcopy(stage)
        initial_pass+=stage["passed"]
        if stage["shard"]==ORIGINAL_FAILED_SHARD:
            failures=[r for r in stage["outcomes"] if r["status"]=="failed"]
            if len(failures)!=1 or failures[0].get("run")!=RUN or \
                failures[0].get("biosample")!=ORIGINAL_BIO or \
                failures[0].get("accepted_key")!=ORIGINAL_SOURCE_ACCEPTED or \
                failures[0].get("error")!=EXPECTED_FAILURE:
                raise ValueError("Original source failure not the one narrowly allowed retry")
            original_fail=failures
            updated["outcomes"]=[
                rescue if z["run"]==RUN else z for z in updated["outcomes"]]
            updated["failed"]=0
            updated["passed"]+=1
            updated["original_failed_source_run_sra_prefix_proof_sha256"]=digest["sha256"]
            updated["only_single_targeted_source_Q20_repair_applied"]=True
        elif stage["failed"]:
            raise ValueError("An unexpected other original source accession failed")
        records.extend(updated["outcomes"])
        repaired.append(updated)
    if len(original_fail)!=1 or initial_pass!=SOURCE_SUCCESS or \
       sum(x["passed"] for x in repaired)!=TOTAL_TAXA or \
       sum(x["failed"] for x in repaired)!=0:
        raise ValueError("Cannot claim all111 genuine Q20 site calls recovered")
    names=[x["run"] for x in records]
    if len(names)!=TOTAL_TAXA or len(set(names))!=TOTAL_TAXA or \
       set(names)!={r["run"] for r in selected}:
        raise ValueError("Original 111 real source accession inventory not conserved")
    original_success={
       (x["run"],json.dumps(x,sort_keys=True))
       for stage,_ in original_stages for x in stage["outcomes"]
       if x["status"].startswith("real_")}
    retained_success={
       (x["run"],json.dumps(x,sort_keys=True))
       for x in records if x["run"]!=RUN}
    if original_success!=retained_success:
        raise ValueError("One of previously successful original 110 Q20 observations was altered")
    audit={
       "schema":"azami_v4_111_120k_sra_original_site_receipt_one_run_audited_recovery_20261023",
       "old_action_run":39888260140,
       "old_original_nucleotide_source_successes":SOURCE_SUCCESS,
       "old_original_failed_SRA_runs":[RUN],
       "new_real_source_sra_successes":TOTAL_TAXA,
       "successfully_reused_unchanged_110_original_Q20_call_receipts":True,
       "new_original_public_Q20_source_call_receipts":1,
       "larger_gzip_prefix_but_same120k_read_spots":True,
       "same_118_source_original_nuclear_reference_copies":True,
       "no_reference_base_imputation":True,
       "original_failed_shard_was_not_overwritten_in_GitHub_actions":True,
       "all_244_WCVP_keys_and_251_source_name_units_preserved_outside_sequence_panel":True,
       "accepted_nuclear_species_tree_constructed":False,
       "original_10_shard_source_json_sha256":[x for _,x in original_stages],
       "source_original_failure":original_fail[0],
    }
    return repaired,audit


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--failed-original-shards",type=Path,required=True)
    ap.add_argument("--recovered-run",type=Path,required=True)
    ap.add_argument("--samples",type=Path,required=True)
    ap.add_argument("--targets",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args()
    original=load_stage_files(args.failed_original_shards)
    rescue=json.loads(args.recovered_run.read_text(encoding="utf-8"))
    samples=nuclear.load_csv(args.samples,nuclear.SOURCE_TAXA_FIELDS)
    refs=nuclear.load_csv(args.targets,nuclear.TARGET_FIELDS)
    repaired,proof=repair(original,rescue,samples,refs)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    for stage in repaired:
        (args.out_dir/f"shard_{stage['shard']}.json").write_text(
           json.dumps(stage,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n",
           encoding="utf-8")
    (args.out_dir/"REAL111_ORIGINAL_SOURCE_ONE_SRA_RECOVERY_AUDIT.json").write_text(
       json.dumps(proof,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
       encoding="utf-8")
    print("T2_EXACT_ONE_RUN_RESCUE_111_TRUE_DNA_SOURCES="+json.dumps({
        "original_unmodified_samples":SOURCE_SUCCESS,
        "recovered_run":RUN,"restored_source_samples":TOTAL_TAXA,
        "no_imputation":proof["no_reference_base_imputation"]},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
