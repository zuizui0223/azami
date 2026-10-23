"""Merge ONLY 3 genuine network-error corrections into 10 original111 Q20 shards.

Never recompute or alter the original 108 already-successful public nuclear
DNA samples. Preserve original failed shard records, exact original evidence
digests, and all run/BioSample/WCVP/locus source checks. Repaired output is a
distinct derived artifact, not an overwrite of the original GitHub Actions
#39888192062 receipts.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

from analysis.v4 import t2_real111_original118_Q20_nucleotide_calls as nuc
from analysis.v4 import t2_real111_original118_shared_Q20_aggregate as agg
from analysis.v4 import t2_real111_three_transient_ENA_retry as fix


def patch(original,repairs,samples,targets):
    shards=agg.load_shards(original)
    samples=nuc.read_csv(samples)
    targets=nuc.read_csv(targets)
    allowed={x["run"]:x for x in samples}
    if len(allowed)!=111 or len(targets)!=118:
        raise ValueError("One of 111 original source nuclear samples or 118 loci missing")
    byfailed={}
    untouched=[]
    for shard in shards:
        for x in shard["outcomes"]:
            if x["status"]=="failed":
                byfailed[x["source_original_run"]]=x
            else:
                untouched.append(json.dumps(x,sort_keys=True))
    if set(byfailed)!=fix.TRANSIENT_RUNS or len(untouched)!=108:
        raise ValueError("Original 108 success / 3 transient ENA error receipts changed")
    paths=sorted(Path(repairs).rglob("SRR*.json"))
    if len(paths)!=3:
        raise ValueError("Require exactly three re-run actual SRA source JSON receipts")
    replacements={}
    for path in paths:
        r=json.loads(path.read_text(encoding="utf-8"))
        run=r["original_failed_run"]
        if (run not in byfailed or run in replacements or
            r["status"]!="remediated_external_ENA_transport_failure" or
            not r["predecessor_108_success_shards_not_recomputed"] or
            r["original_workflow"]!=fix.ORIGINAL_FAILED_WORKFLOW):
            raise ValueError("Replacement attempts to promote unfailed original BioSample")
        x=r["result"]
        original_case=byfailed[run]
        author=allowed[run]
        if (r["original_biosample"]!=author["biosample"] or
            original_case["source_BioSample"]!=author["biosample"] or
            x["source_BioSample"]!=author["biosample"] or
            x["source_S1_accepted_key"]!=author["accepted_key"] or
            x["tip"]!=author["phylogenetic_tip"] or
            x["status"]!="real_Q20_majority_reference_site_calls" or
            set(x["loci"])!={g["locus"] for g in targets} or
            x["n_real_read_pairs"]!=120000):
            raise ValueError("Changed run/BioSample/WCVP/118-locus Q20 source identity")
        replacements[run]=x
    if set(replacements)!=fix.TRANSIENT_RUNS:
        raise ValueError("Not all three exact source ENA transient corrections recovered")
    result=[]
    for old in shards:
        clone=json.loads(json.dumps(old))
        for i,x in enumerate(clone["outcomes"]):
            if x["source_original_run"] in replacements:
                clone["outcomes"][i]=replacements[x["source_original_run"]]
        clone["failure"]=0
        clone["success"]=len(clone["outcomes"])
        clone["repair_only_known_external_ENA_source_runs"]=sorted(
            set(x["source_original_run"] for x in old["outcomes"])&fix.TRANSIENT_RUNS)
        result.append(clone)
    retained=[json.dumps(x,sort_keys=True) for shard in result
              for x in shard["outcomes"] if x["source_original_run"] not in fix.TRANSIENT_RUNS]
    if retained!=untouched:
        raise ValueError("Original 108 successful measured basecalls were unexpectedly modified")
    audit={
        "status":"ONLY_THREE_REAL_EXTERNAL_ENA_ERRORS_REPAIRED_IN_ORIGINAL111",
        "original_source_workflow":fix.ORIGINAL_FAILED_WORKFLOW,
        "original_nuclear_BioSample_Q20_success_preserved":len(untouched),
        "transient_public_ENA_run_failures_fixed":sorted(replacements),
        "final_original_source_BioSamples":sum(len(s["outcomes"]) for s in result),
        "final_source_run_failures":0,
        "original_S1_WCVP_taxa_master":244,
        "original_S1_source_name_units":251,
        "original_author_no_warning_target_universe":241,
        "original_same_reference_target_candidates":118,
        "site_quality_cutoffs_or_original_run_changed":False,
        "n_fastq_bam_distributed":0,
        "orthology_or_species_tree_confirmed":False,
    }
    return result,audit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original",type=Path,required=True)
    p.add_argument("--repairs",type=Path,required=True)
    p.add_argument("--samples",type=Path,required=True)
    p.add_argument("--targets",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    repaired,summary=patch(args.original,args.repairs,args.samples,args.targets)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    for r in repaired:
        (args.out_dir/f"shard_{r['shard']}.json").write_text(
            json.dumps(r,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n",
            encoding="utf-8")
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(args.original).rglob("shard_*.json"))}
    summary["original_10_unchanged_SRA_shard_json_SHA256"]=hashes
    (args.out_dir/"REMEDIATED_ONLY3_ENA_TRANSPORT_AUDIT.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n",
        encoding="utf-8")
    print("T2_ORIGINAL111_ONLY3_NETWORK_RUN_REMEDIATION="+json.dumps({
        k:summary[k] for k in ("original_nuclear_BioSample_Q20_success_preserved",
            "transient_public_ENA_run_failures_fixed","final_original_source_BioSamples",
            "final_source_run_failures")},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
