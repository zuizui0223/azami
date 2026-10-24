"""Source-evidence scheduling checks. No SRA reads or fake locus alignment."""
from collections import Counter
import pytest

from analysis.v4 import t2_moreyra_capacity_shards as mod


def tiny_input():
    sizes=[("x","Cirsium_a","SRR001",4200),
           ("x","Cirsium_a","SRR002",2800),
           ("y","Cirsium_b","SRR003",800),
           ("z","Cirsium_c","SRR004",900),
           ("w","Cirsium_d","SRR005",2200)]
    return [{"S1_accepted_key":key,"S1_accepted_key_unused":"ignore",
             "source_taxon":name,"SRA_run":acc,"BioSample":"SAMN"+acc[3:],
             "archive_MB":size,"source_taxon_name_disagreement":"False"}
            for key,name,acc,size in sizes]


def test_all_runs_and_species_preserved_no_accidental_oversize_omission(monkeypatch):
    monkeypatch.setattr(mod,"RUN_N",5)
    monkeypatch.setattr(mod,"TAXON_N",4)
    groups,queue,report=mod.schedule(tiny_input(),cap=4000)
    assert len(groups)==4
    assert len(queue)==5
    assert len({x["SRA_run"] for x in queue})==5
    assert report["compressed_SRA_size_total_MB"]==10900
    assert report["n_oversized_species_groups"]==1
    assert report["oversized_groups"][0]["compressed_MB"]==7000
    assert len({r["shard"] for r in queue if r["accepted_key"]=="x"})==1
    assert all(not x["FASTQ_processed"] for x in queue)
    assert report["original_public_350_locus_tree_reproduced"] is False


def test_missing_species_is_a_hard_error(monkeypatch):
    monkeypatch.setattr(mod,"TAXON_N",5)
    monkeypatch.setattr(mod,"RUN_N",5)
    with pytest.raises(ValueError,match="113 S1 candidate taxa"):
        mod.schedule(tiny_input(),cap=4000)


def test_original_target_is_all_125_runs_not_random_small_subset():
    assert mod.RUN_N==125
    assert mod.TAXON_N==113
    assert mod.SUM_ARCHIVE_MB==57630.0
    assert mod.DEFAULT_GROUP_CAP_MB==4000
