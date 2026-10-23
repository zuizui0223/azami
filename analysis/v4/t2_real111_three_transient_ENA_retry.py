"""Retry ONLY 3 source-authenticated 111-run Q20 failures caused by public ENA transport.

The original 108 completed source-run nucleotides remain archived, unchanged.
Only run IDs which actually failed in GitHub Actions #39888192062 may be
reanalyzed, preserving original EAzami BioSample, author-241 and predeclared
118 single-reference copies, 120k real public paired FASTQ source.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as original
from analysis.v4 import t2_moreyra_all125_author241_site_depth as depth
from analysis.v4 import t2_real111_original118_Q20_nucleotide_calls as nuc

TRANSIENT_RUNS={"SRR25265719","SRR30887150","SRR30887204"}
ORIGINAL_FAILED_WORKFLOW=39888192062
MAX_ATTEMPTS=3
NETWORK_EXCEPTIONS=("RuntimeError","RemoteDisconnected","URLError",
                    "TimeoutError","ConnectionResetError")


def repair(one,raw,targets,loci,reference,scratch,runner=nuc.source_call,
           sleeper=time.sleep):
    if one["run"] not in TRANSIENT_RUNS:
        raise ValueError("Do not rerun already-successful 108 public BioSample read cases")
    caught=[]
    for attempt in range(1,MAX_ATTEMPTS+1):
        try:
            r=runner(one,raw,targets,loci,reference,scratch)
            if (r.get("source_original_run")!=one["run"] or
                r.get("source_BioSample")!=one["biosample"] or
                r.get("tip")!=one["phylogenetic_tip"] or
                len(r["loci"])!=len(targets)):
                raise ValueError("Q20 real source nucleotide retry changed public specimen")
            r["status"]="real_Q20_majority_reference_site_calls"
            r["elapsed_seconds"]=None
            return {
                "status":"remediated_external_ENA_transport_failure",
                "original_workflow":ORIGINAL_FAILED_WORKFLOW,
                "original_failed_run":one["run"],
                "original_biosample":one["biosample"],
                "retry_attempts_used":attempt,
                "predecessor_108_success_shards_not_recomputed":True,
                "no_change_to_original_Q20_or_120k_sampling_gate":True,
                "not_a_species_tree":True,
                "result":r,
            }
        except Exception as err:
            # Only transient transfer/network exceptions may be retried.
            kind=type(err).__name__
            if kind not in NETWORK_EXCEPTIONS:
                raise
            if kind=="RuntimeError" and not any(
                x in str(err) for x in ("ENA public source access","Remote",
                                       "HTTP","Connection","Timeout")):
                raise
            caught.append({"try":attempt,"exception":kind,"reason":str(err)[:180]})
            if attempt<MAX_ATTEMPTS:
                sleeper(5*attempt)
    raise RuntimeError(
        "Source ENA transient retry exhausted for "+one["run"]+
        " with exact run/BioSample source frozen: "+json.dumps(caught))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run",required=True)
    for x in ("--source","--samples","--targets","--author241","--reference",
              "--scratch","--out"):
        p.add_argument(x,type=Path,required=True)
    a=p.parse_args()
    if a.run not in TRANSIENT_RUNS:
        raise ValueError("Not one of the exact three observed ENA transport failures")
    sources=original.load_source(a.source)
    samples=nuc.read_csv(a.samples)
    targets=nuc.read_csv(a.targets)
    selected,refs=nuc.source_inputs(sources,samples,targets)
    found=[(one,raw) for one,raw in selected if one["run"]==a.run]
    if len(found)!=1:
        raise ValueError("Observed external failure run missing from original 111 source panel")
    loci=depth.loci_from_pinned_author_list(a.author241)
    if not a.reference.is_file():
        raise ValueError("Original source nuclear Compositae1061 reference absent")
    result=repair(*found[0],refs,loci,a.reference,a.scratch)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,
                                indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_REPAIRED_ONLY_ORIGINAL_ENA_RUN="+json.dumps({
        "run":a.run,"attempts":result["retry_attempts_used"],
        "n_Q20_sites":result["result"]["n_Q20_majority_called_sites"]},
        sort_keys=True),flush=True)


if __name__=="__main__":
    main()
