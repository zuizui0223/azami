"""Aggregated numerical depth-screen tests with synthetic sample outputs."""
import json

import pytest

from analysis.v4 import t2_moreyra_author241_depth_aggregate as m


def row(run,key,genes=3):
    locus=[dict(locus=f"L{i}",reference_choice_maximizes_coverage_not_orthology=True,
         dp1_breadth=.8,dp3_breadth=.65,dp5_breadth=.3,
         dp3_at_least100bp=True) for i in range(genes)]
    counts={
      "n_loci_dp1_any":3,"n_loci_dp3_any":3,"n_loci_dp5_any":3,
      "n_loci_dp3_100bp":3,"n_loci_dp3_ge10pct":3,
      "n_loci_dp3_ge25pct":3,"n_loci_dp3_ge50pct":3,
      "n_loci_dp3_ge75pct":0,
      "total_best_reference_sites_dp1":600,
      "total_best_reference_sites_dp3":480,
      "total_best_reference_sites_dp5":260,
    }
    return {"status":"real_ENA_reference_site_depth_q20_Q20",
       "run":run,"biosample":f"SAMN{run[3:]}",
       "accepted_key":key,"S1_source_taxon":"Cirsium_sample",
       "real_ena_read_pairs":30000,"loci":locus,**counts}


def test_real_shard_source_provenance_and_241_ortholog_claim_gate(monkeypatch):
    monkeypatch.setattr(m,"RUNS",3)
    monkeypatch.setattr(m,"N_WCVP",2)
    monkeypatch.setattr(m,"LOCI",3)
    r=[row("SRR1","x"),row("SRR2","x"),row("SRR3","y")]
    source=[dict(run=z["run"],biosample=z["biosample"],
                 accepted_key=z["accepted_key"]) for z in r]
    reports=[{"requested":3,"outcomes":r}]
    full=m.validate(source,reports)
    per_run,taxa,loci,summary=m.summarize(source,full)
    assert len(full)==3 and len(per_run)==3
    assert len(taxa)==2 and len(loci)==3
    assert sorted(x["n_independent_public_runs"] for x in taxa)==[1,2]
    assert summary["strict_241_ortholog_alignments_completed"] is False
    assert summary["eazami294_astral_species_tree_completed"] is False
    assert summary["run_level_metrics_not_independent_species"]


def test_missing_source_run_blocks_success(monkeypatch):
    monkeypatch.setattr(m,"RUNS",3)
    monkeypatch.setattr(m,"N_WCVP",2)
    monkeypatch.setattr(m,"LOCI",3)
    src=[{"run":f"SRR{i}","biosample":f"SAMN{i}",
          "accepted_key":"a" if i<3 else "b"} for i in (1,2,3)]
    receipts=[{"requested":2,"outcomes":[row("SRR1","a"),row("SRR2","a")]}]
    with pytest.raises(ValueError,match="all 125"):
        m.validate(src,receipts)


def test_missing_ref_gene_or_false_full_sequence_metadata_blocks(monkeypatch):
    monkeypatch.setattr(m,"RUNS",3)
    monkeypatch.setattr(m,"N_WCVP",2)
    monkeypatch.setattr(m,"LOCI",3)
    sample=[row("SRR1","x"),row("SRR2","x"),row("SRR3","y")]
    src=[dict(run=z["run"],biosample=z["biosample"],accepted_key=z["accepted_key"])
         for z in sample]
    sample[0]["loci"][0]["reference_choice_maximizes_coverage_not_orthology"]=False
    with pytest.raises(ValueError,match="best-reference"):
        m.validate(src,[{"requested":3,"outcomes":sample}])


def test_requires_exact_ten_shard_set(tmp_path):
    (tmp_path/"shard_0.json").write_text(json.dumps({
        "shard_id":0,"n_shards":10,"read_data_in_artifact":False,
        "BAM_in_artifact":False}))
    with pytest.raises(ValueError,match="ten original"):
        m.load_receipts(tmp_path)
