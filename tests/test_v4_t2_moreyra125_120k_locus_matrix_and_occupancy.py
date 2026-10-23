"""Source-contract unit tests for real 120k x 125 run x 241 locus recovery.

These tests use synthetic fixtures ONLY to validate bookkeeping. Empirical
numbers must come from the pinned public-SRA GitHub Actions run.
"""
import json
import pytest

from analysis.v4 import t2_moreyra125_120k_source241_locus_matrix as matrix
from analysis.v4 import t2_moreyra125_120k_shared_locus_occupancy as occupancy


def source_fixture():
    source = [
        {"run": "SRR100", "biosample": "SAMN100", "accepted_key": "x",
         "S1_source_taxon": "Cirsium testone"},
        {"run": "SRR101", "biosample": "SAMN101", "accepted_key": "y",
         "S1_source_taxon": "Cirsium testtwo"},
        {"run": "SRR102", "biosample": "SAMN102", "accepted_key": "z",
         "S1_source_taxon": "Cirsium testthree"},
    ]
    first = {"run": "SRR100", "biosample": "SAMN100", "accepted_key": "x",
             "full_120k_comparison_status": "real_ENA_120k_quality20_depth",
             "source_30k_to_120k_valid": "True"}
    for metric in matrix.METRICS:
        first["first120k_" + metric] = "2" if metric.startswith("n_loci") else "750"
    return source, [first, {
        "run": "SRR101", "biosample": "SAMN101", "accepted_key": "y",
        "full_120k_comparison_status": "failed",
        "source_run_unavailability_reason": "archive has fewer than 120k pairs",
    }, {
        "run": "SRR102", "biosample": "SAMN102", "accepted_key": "z",
        "full_120k_comparison_status": "failed",
        "source_run_unavailability_reason": "archive has fewer than 120k pairs",
    }]


def test_validate_original_120k_identity_and_keep_short_runs(monkeypatch):
    monkeypatch.setattr(matrix, "RUNS", 3)
    monkeypatch.setattr(matrix, "TAXA", 3)
    monkeypatch.setattr(matrix, "LOCI", 2)
    monkeypatch.setattr(matrix, "PRIOR_UNAVAILABLE", {"SRR101", "SRR102"})
    src, old = source_fixture()
    matched = matrix.certify_source(src, old, ["L1", "L2"])
    assert set(matched) == {x["run"] for x in src}
    old[0]["biosample"] = "SAMN_STALE"
    with pytest.raises(ValueError, match="BioSample"):
        matrix.certify_source(src, old, ["L1", "L2"])


def test_observed_source_reads_not_equivalent_to_orthologs(monkeypatch, tmp_path):
    monkeypatch.setattr(matrix, "PAIRS", 2)
    monkeypatch.setattr(matrix, "LOCI", 2)
    src, old = source_fixture()
    def acquisition(_row):
        return ([[b"@a\n",b"A\n",b"+\n",b"I\n"] * 2] * 2, {})
    def mapper(_row, _reference, _loci, _reads, _scratch, include_per_locus=False):
        assert include_per_locus
        data = [
            {"locus": name, "best_original_reference": "lett-" + name,
             "dp1_breadth": .9, "dp3_breadth": .5, "dp5_breadth": .2,
             "reference_choice_maximizes_coverage_not_orthology": True}
            for name in ("L1", "L2")
        ]
        return {
            "input_read_pairs":2,
            "n_loci_dp3_ge25pct":2,"n_loci_dp3_ge50pct":2,
            "n_loci_dp3_100bp":2,"total_best_reference_sites_dp3":750,
            "per_locus_241":data,
        }
    out = matrix.measure_one(src[0],old[0],["L1","L2"],None,tmp_path,
                             acquisition=acquisition,mapper=mapper)
    assert out["run"] == "SRR100"
    assert len(out["per_locus_241"]) == 2
    assert out["true_phylogenetic_orthology_certified"] is False
    assert out["individual_gene_consensus_reconstructed"] is False


def test_source_120k_error_not_silently_promoted(monkeypatch,tmp_path):
    monkeypatch.setattr(matrix,"RUNS",3)
    monkeypatch.setattr(matrix,"TAXA",3)
    monkeypatch.setattr(matrix,"LOCI",2)
    monkeypatch.setattr(matrix,"SHARDS",1)
    monkeypatch.setattr(matrix,"PRIOR_UNAVAILABLE",{"SRR101","SRR102"})
    src,old=source_fixture()
    def empirical_simulation(row, old, loci, ref, workdir):
        return {"run":row["run"],"biosample":row["biosample"],
                "accepted_key":row["accepted_key"],"status":"real_ENA_120k_same_run_author241_locus_Q20",
                "n_loci_dp3_ge50pct":2}
    report=matrix.execute(src,old,["L1","L2"],None,tmp_path,0,1,
                          runner=empirical_simulation)
    assert len(report["outcomes"])==3
    assert report["mapped"]==1 and report["original_short_runs"]==2
    assert report["unexpected_failures"]==0
    assert report["raw_FASTQ_or_BAM_in_artifact"] is False


def test_shared_occupancy_requires_ten_source_shards(tmp_path):
    (tmp_path/"shard_0.json").write_text(json.dumps({
        "shard":0,"total_shards":10,"raw_FASTQ_or_BAM_in_artifact":False,
        "requested":1,"outcomes":[]}))
    with pytest.raises(ValueError,match="Ten 120k"):
        occupancy.read_shards(tmp_path)


def test_reference_copy_gate_not_same_as_verified_orthology():
    assert occupancy.PRIMARY_COVERAGE == .8
    assert occupancy.SHORT_RUNS=={"SRR30887185","SRR30887248"}
    assert len(occupancy.NODDING_SOURCE_NAMES)==2
    assert matrix.PAIRS==120000
