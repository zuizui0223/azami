"""Test one-run exact public-source repair without silently rerunning/rewriting 110."""
from copy import deepcopy
import pytest

from analysis.v4 import t2_real111_audited_one_run_repair as r
from analysis.v4 import t2_real111_fixed118_nucleotides as q


def synthetic(monkeypatch):
    monkeypatch.setattr(r,"SHARDS",2)
    monkeypatch.setattr(r,"TOTAL_TAXA",3)
    monkeypatch.setattr(r,"SOURCE_SUCCESS",2)
    monkeypatch.setattr(r,"ORIGINAL_FAILED_SHARD",1)
    monkeypatch.setattr(q,"TARGET_COUNT",2)
    run=r.RUN
    samples=[{
        "run":run,"biosample":r.ORIGINAL_BIO,"accepted_key":r.ORIGINAL_SOURCE_ACCEPTED,
        "taxon_name":r.TAXON,
        "phylogenetic_tip":"WCVP"+r.ORIGINAL_SOURCE_ACCEPTED},
        {"run":"SRR100","biosample":"SAMN100","accepted_key":"100",
         "taxon_name":"Cirsium example100","phylogenetic_tip":"WCVP100"},
        {"run":"SRR200","biosample":"SAMN200","accepted_key":"200",
         "taxon_name":"Cirsium example200","phylogenetic_tip":"WCVP200"}]
    targets=[{"locus":l,"original_compositae_reference":"saff-"+l}
             for l in ("L1","L2")]
    data={
        **samples[0],"status":"real_source111_120k_Q20_identical_copy_base_calls",
        "source_public_ENA_paired_read_count":120000,
        "n_predeclared_80percent_reference_loci":2,
        "n_real_Q20_majority_callable_positions":2,
        "observed_original_reference_q20_site_calls":{
           l:{"reference":"saff-"+l,"reference_length_nt":100,
              "not_a_validated_allelic_genotype":True,
              "q20_majority_bases_1based":[[1,"A",3,1.0]]}
           for l in ("L1","L2")},
        "uncalled_bases_not_filled_by_reference":True,
        "original_single_copy_orthology_verified":False,
        "FASTQ_BAM_exported":False,
        "rescue_protocol":"same_real_first120000_pairs_larger_compressed_HTTP_prefix_only",
        "rescue_max_compressed_bytes_per_mate":128*1024*1024}
    successes=[
        {"run":"SRR100","biosample":"SAMN100","accepted_key":"100",
         "status":"real_source111_120k_Q20_identical_copy_base_calls"},
        {"run":"SRR200","biosample":"SAMN200","accepted_key":"200",
         "status":"real_source111_120k_Q20_identical_copy_base_calls"},
    ]
    failure={"run":run,"biosample":r.ORIGINAL_BIO,
              "accepted_key":r.ORIGINAL_SOURCE_ACCEPTED,"status":"failed",
              "error":r.EXPECTED_FAILURE}
    stages=[
        ({"shard":0,"n_shards":2,"requested":2,"passed":2,"failed":0,
          "no_original_FASTQ_or_BAM_preserved":True,"outcomes":successes},
         {"filename":"shard_0.json","sha256":"abc"}),
        ({"shard":1,"n_shards":2,"requested":1,"passed":0,"failed":1,
          "no_original_FASTQ_or_BAM_preserved":True,"outcomes":[failure]},
         {"filename":"shard_1.json","sha256":"xyz"}),
    ]
    return stages,data,samples,targets


def test_real_public_110_prior_successes_preserved_bit_identically(monkeypatch):
    original,rescue,samples,refs=synthetic(monkeypatch)
    pristine=deepcopy(original)
    out,manifest=r.repair(original,rescue,samples,refs)
    assert sum(stage["passed"] for stage in out)==3
    assert all(stage["failed"]==0 for stage in out)
    assert out[0]["outcomes"]==pristine[0][0]["outcomes"]
    assert original==pristine
    assert out[1]["outcomes"]==[rescue]
    assert manifest["successfully_reused_unchanged_110_original_Q20_call_receipts"]
    assert manifest["new_original_public_Q20_source_call_receipts"]==1
    assert not manifest["accepted_nuclear_species_tree_constructed"]


def test_original_fail_exact_identity_cannot_be_swapped(monkeypatch):
    original,rescue,samples,refs=synthetic(monkeypatch)
    original[1][0]["outcomes"][0]["run"]="SRRFAKE"
    with pytest.raises(ValueError,match="narrowly allowed retry"):
        r.repair(original,rescue,samples,refs)


def test_rescue_must_be_real_Q20_observations_and_same_biosample(monkeypatch):
    original,rescue,samples,refs=synthetic(monkeypatch)
    rescue["biosample"]="SAMNFAKE"
    with pytest.raises(ValueError,match="exact original source"):
        r.repair(original,rescue,samples,refs)
    original,rescue,samples,refs=synthetic(monkeypatch)
    rescue["observed_original_reference_q20_site_calls"]["L1"][
        "q20_majority_bases_1based"][0][2]=2
    with pytest.raises(ValueError,match="Invalid positionwise"):
        r.repair(original,rescue,samples,refs)


def test_repair_public_evidence_constants_immutable():
    assert r.RUN=="SRR30887180"
    assert r.ORIGINAL_BIO=="SAMN44017849"
    assert r.ORIGINAL_SOURCE_ACCEPTED=="206743475"
    assert r.SHARDS==10 and r.SOURCE_SUCCESS==110 and r.TOTAL_TAXA==111
