"""Whole-locus IQ-TREE ML bootstrap on real seven-source Q20 bases: numerical gates."""
from io import StringIO

from Bio import Phylo
import pytest

from analysis.v4 import t2_real7_source101_locus_block_ml_bootstrap as m

NAMES=[
 "Cirsium_flodmanii","Cirsium_tracyi","Cirsium_scopulorum","Cirsium_vulgare",
 "Cirsium_italicum","Cirsium_latifolium","Cirsium_pulcherrimum"
]
TREE=(
 "(((Cirsium_flodmanii:0.1,Cirsium_tracyi:0.1):0.01,"
 "(Cirsium_scopulorum:0.1,Cirsium_vulgare:0.1):0.01):0.02,"
 "(Cirsium_italicum:0.1,(Cirsium_latifolium:0.1,"
 "Cirsium_pulcherrimum:0.1):0.01):0.02);"
)

def fake_loci(n=101):
    return [{"locus":f"L{i}","alignment":["ACGT"*(i%3+1) for _ in NAMES]}
            for i in range(n)]

def test_real_source_blocks_are_resampled_as_whole_loci(monkeypatch):
    monkeypatch.setattr(m,"CANDIDATE_LOCI",3)
    loci=fake_loci(3)
    first,chosen=m.bootstrap_alignment(NAMES,loci,20261023)
    second,chosen2=m.bootstrap_alignment(NAMES,loci,20261023)
    assert first==second
    assert [x["locus"] for x in chosen]==[x["locus"] for x in chosen2]
    assert len(chosen)==3
    assert all(len(first[name])==len(next(iter(first.values()))) for name in NAMES)
    assert all(x["locus"] in {"L0","L1","L2"} for x in chosen)

def test_exact_2_shard_resampling_keeps_all_300_draws_without_duplication():
    loci=fake_loci()
    def runner(i,names,locs,workdir):
        return {
          "bootstrap_locus_draw_id":i,
          "unrooted_ML_split_labels":[f"S{k}" for k in range(4)],
          "original_7_SRA_source_SNPs_only":True,
          "full_locus_sequences_reconstructed":False,
        }
    shard0=m.run_shard(NAMES,loci,0,2,None,runner=runner)
    shard1=m.run_shard(NAMES,loci,1,2,None,runner=runner)
    assert shard0["n_results"]==150
    assert shard1["n_results"]==150
    ids={x["bootstrap_locus_draw_id"] for s in (shard0,shard1) for x in s["replicates"]}
    assert ids==set(range(300))

def test_published_4taxon_exact_quartet_matched_in_every_synthetic_bootstrap():
    original=Phylo.read(StringIO(TREE),"newick")
    sides=m.split_labels(original)
    sh=[]
    for shard in (0,1):
        reps=[{
            "bootstrap_locus_draw_id":k,"unrooted_ML_split_labels":sides,
            "original_7_SRA_source_SNPs_only":True,
            "full_locus_sequences_reconstructed":False,
        } for k in range(shard,300,2)]
        sh.append({
            "shard_id":shard,"n_original_real_SRA_samples":7,
            "n_original_101_candidate_loci":101,"source_only_no_simulated_bases":True,
            "replicates":reps,
        })
    result=m.aggregate(sh,original,original)
    assert result["n_total_independently_reestimated_ML_replicates"]==300
    assert result["ML_locus_block_published_NA63_quartet_fraction"]==1.
    assert all(x["ML_whole_locus_bootstrap_fraction"]==1
               for x in result["whole_locus_ML_branch_frequency"])
    assert not result["uncertainty"]["species_phylogeny_rooted_or_dated"]
    assert not result["uncertainty"]["global_EAzami294_species_phylogeny_completed"]

def test_block_bootstrap_requires_two_valid_shards():
    with pytest.raises(ValueError,match="Missing paired"):
        m.aggregate([],None,None)

def test_no_random_short_sequence_accepted_as_101_loci():
    with pytest.raises(ValueError,match="101 source loci"):
        m.bootstrap_alignment(NAMES,[],20261023)
