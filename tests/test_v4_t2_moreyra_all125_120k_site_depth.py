"""Tests for complete source-preserving 125-run depth 4x scaling."""
import pytest
from analysis.v4 import t2_moreyra_all125_120k_site_depth as mod
from analysis.v4 import t2_moreyra_all125_120k_depth_aggregate as agg


def original(monkeypatch):
    monkeypatch.setattr(mod,"RUNS",4)
    monkeypatch.setattr(mod,"TAXA",3)
    monkeypatch.setattr(mod,"SHARDS",2)
    return [{"run":f"SRR00000{i}","biosample":f"SAMN00000{i}",
             "accepted_key":f"K{i%3}","S1_source_taxon":f"Cirsium_taxon{i%3}"}
            for i in range(4)]


def test_four_sample_scaling_treats_unavailable_as_not_zero_genes(monkeypatch,tmp_path):
    source=original(monkeypatch)
    monkeypatch.setattr(mod.depth,"LOCUS_COUNT",1)
    def acquire(row):
        if row["run"]=="SRR000002":
            raise ValueError("ENA source prefix cannot provide 120k")
        mate=[b"@x\n",b"A\n",b"+\n",b"I\n"]*120000
        return [mate,mate],[{},{}]
    def mapreads(row,ref,loci,reads,path):
        return {"input_read_pairs":120000,
                **{key:100 for key in mod.METRICS}}
    out=mod.execute(source,None,["L1"],tmp_path,0,nshards=2,
                    acquire=acquire,mapper=mapreads)
    assert out["requested"]==2
    assert out["successes"]==1 and out["failures"]==1
    assert not out["nuclear_tree_generated"]
    assert not out["raw_FASTQ_or_BAM_saved_in_artifact"]
    assert out["outcomes"][1]["status"]=="failed"


def test_120k_scoring_accepts_same_original_125_biosamples_not_arbitrary_data(monkeypatch):
    monkeypatch.setattr(agg,"RUNS",2)
    monkeypatch.setattr(agg,"TAXA",1)
    original=[{
         "run":run,"biosample":f"SAMN{i}","accepted_key":"K0",
         "S1_source_taxon":"Cirsium_vulgare","real_ena_read_pairs":"30000",
         **{key:("10" if key.startswith("n_loci_") else "100") for key in agg.DEPTH_KEYS}}
         for i,run in enumerate(("SRR1","SRR2"),1)]
    reports=[{"outcomes":[{
             "run":run,"biosample":f"SAMN{i}","accepted_key":"K0",
             "status":"real_ENA_120k_quality20_depth","read_pairs":120000,
             **{key:(22 if key.startswith("n_loci_") else 200) for key in agg.DEPTH_KEYS}}]}
         for i,run in enumerate(("SRR1","SRR2"),1)]
    res=agg.reconcile(original,reports)
    taxa,summary=agg.summarize(res)
    assert len(res)==2 and len(taxa)==1
    assert summary["n_runs_with_recovered_120k_read_depth"]==2
    assert summary["n_paired_successes_with_higher_depth3_half_locus_coverage"]==2
    assert taxa[0]["best_120k_single_run_dp3_ge50pct"]==22
    assert not summary["biological_species_tree_completed"]
    reports[0]["outcomes"][0]["biosample"]="SAMN_WRONG"
    with pytest.raises(ValueError,match="Different voucher"):
        agg.reconcile(original,reports)


def test_120k_less_coverage_than_30k_is_not_silent(monkeypatch):
    monkeypatch.setattr(agg,"RUNS",1)
    monkeypatch.setattr(agg,"TAXA",1)
    old=[{"run":"SRR1","biosample":"SAMN1","accepted_key":"K",
           "S1_source_taxon":"Cirsium_x","real_ena_read_pairs":"30000",
           **{key:"23" for key in agg.DEPTH_KEYS}}]
    new=[{"outcomes":[{"run":"SRR1","biosample":"SAMN1","accepted_key":"K",
          "status":"real_ENA_120k_quality20_depth","read_pairs":120000,
          **{key:15 for key in agg.DEPTH_KEYS}}]}]
    with pytest.raises(ValueError,match="lost Q20"):
        agg.reconcile(old,new)
