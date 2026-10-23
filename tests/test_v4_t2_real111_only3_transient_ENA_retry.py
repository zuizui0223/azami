"""Rehabilitate only three actual public ENA network failures; never redo 108 runs."""
import json
from pathlib import Path
import pytest
from analysis.v4 import t2_real111_merge_only_three_ENA_retries as merge
from analysis.v4 import t2_real111_three_transient_ENA_retry as retry


def fixtures(tmp_path):
    old=tmp_path/"original";new=tmp_path/"repairs"
    old.mkdir();new.mkdir()
    failed=sorted(retry.TRANSIENT_RUNS)
    runs=failed+[f"SRR55000{i:04d}" for i in range(108)]
    assert len(set(runs))==111
    samples=[]
    for i,run in enumerate(runs):
        samples.append({"run":run,"biosample":f"SAMN66000{i:04d}",
                        "accepted_key":str(208800000+i),
                        "phylogenetic_tip":f"WCVP{208800000+i}"})
    allrefs={f"L{i}":{"reference":f"saff-L{i}"} for i in range(118)}
    for shard in range(10):
        batch=[]
        for i,row in enumerate(samples):
            if i%10!=shard:continue
            r={"source_original_run":row["run"],"source_BioSample":row["biosample"],
               "source_S1_accepted_key":row["accepted_key"],
               "tip":row["phylogenetic_tip"],
               "status":"failed" if row["run"] in retry.TRANSIENT_RUNS else
                        "real_Q20_majority_reference_site_calls",
               "n_real_read_pairs":120000,"loci":allrefs,"old_valid":True}
            batch.append(r)
        record={"shard":shard,"n_shards":10,"outcomes":batch,
                "requested":len(batch),
                "success":sum(x["status"]!="failed" for x in batch),
                "failure":sum(x["status"]=="failed" for x in batch),
                "raw_reads_or_BAM_publication":False,
                "single_original_reference_not_confirmed_ortholog":True}
        (old/f"shard_{shard}.json").write_text(json.dumps(record))
    for i,run in enumerate(failed):
        row=next(x for x in samples if x["run"]==run)
        payload={"status":"remediated_external_ENA_transport_failure",
            "original_workflow":retry.ORIGINAL_FAILED_WORKFLOW,
            "original_failed_run":run,"original_biosample":row["biosample"],
            "predecessor_108_success_shards_not_recomputed":True,
            "result":{"source_original_run":run,"source_BioSample":row["biosample"],
                "source_S1_accepted_key":row["accepted_key"],"tip":row["phylogenetic_tip"],
                "status":"real_Q20_majority_reference_site_calls",
                "n_real_read_pairs":120000,"loci":allrefs,"old_valid":False}}
        (new/f"{run}.json").write_text(json.dumps(payload))
    targets=[{"locus":f"L{i}"} for i in range(118)]
    return old,new,samples,targets


def test_exact_three_RETRIED_real_run_merger_never_modifies_108_original(tmp_path):
    old,fix,samples,targets=fixtures(tmp_path)
    prior={f.name:f.read_bytes() for f in old.glob("shard_*.json")}
    shards,summary=merge.patch(old,fix,samples,targets)
    assert len(shards)==10
    assert summary["original_nuclear_BioSample_Q20_success_preserved"]==108
    assert summary["final_original_source_BioSamples"]==111
    assert set(summary["transient_public_ENA_run_failures_fixed"])==retry.TRANSIENT_RUNS
    assert sum(z["failure"] for z in shards)==0
    assert all(f.read_bytes()==prior[f.name] for f in old.glob("shard_*.json"))
    assert not summary["orthology_or_species_tree_confirmed"]


def test_unrequested_run_retry_is_fatal(tmp_path):
    old,fix,samples,targets=fixtures(tmp_path)
    src=sorted(fix.glob("SRR*.json"))[0]
    d=json.loads(src.read_text())
    d["original_failed_run"]="SRR999999999"
    src.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="unfailed"):
        merge.patch(old,fix,samples,targets)


def test_only_actual_public_network_failures_allowed():
    assert retry.MAX_ATTEMPTS==3
    assert retry.ORIGINAL_FAILED_WORKFLOW==39888192062
    assert retry.TRANSIENT_RUNS=={"SRR25265719","SRR30887150","SRR30887204"}


def test_successful_old_original_SRA_run_must_not_be_rerun():
    with pytest.raises(ValueError,match="already-successful"):
        retry.repair({"run":"SRR25265555"},None,None,None,None,None)
