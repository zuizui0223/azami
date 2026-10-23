"""Synthetic source-Q20 masks test whole-locus NJ block resampling mechanics."""
import csv
import pytest

from analysis.v4 import t2_real110_source108_locus_block_nj_bootstrap as m


def small(monkeypatch,tmp_path):
    monkeypatch.setattr(m,"RECOVERED_LOCI",3)
    monkeypatch.setattr(m,"AFTER_SAMPLES",4)
    monkeypatch.setattr(m,"ORIGINAL_SHARED_SITES",12)
    monkeypatch.setattr(m,"MIN_PAIR_Q20_SHARED",2)
    names=[f"WCVP2080{i}" for i in range(4)]
    seq={
        names[0]:"AAAACCCCGGGG",
        names[1]:"AAAACCCCGGGA",
        names[2]:"TTTTGGGGAAAA",
        names[3]:"TTTTGGGGAAAT",
    }
    partitions=[{
      "original_locus":f"L{i}",
      "start_1based":str(4*i+1),
      "end_1based":str(4*i+4),
      "not_single_copy_Cirsium_orthology":"True"} for i in range(3)]
    truth=[]
    for i in range(4):
        for j in range(i):
            one=seq[names[j]];two=seq[names[i]]
            diff=sum(x!=y for x,y in zip(one,two))
            truth.append({"tip_A":names[j],"tip_B":names[i],
                "n_shared_Q20_positions":"12",
                "n_differing_Q20_positions":str(diff),
                "pairwise_deletion_p_distance":str(diff/12)})
    ml=f"(({names[0]}:.01,{names[1]}:.01):.01,({names[2]}:.01,{names[3]}:.01):.01);\n"
    nj=ml
    (tmp_path/"ml.nwk").write_text(ml)
    (tmp_path/"nj.nwk").write_text(nj)
    return seq,partitions,truth,tmp_path/"ml.nwk",tmp_path/"nj.nwk"


def test_exact_original_pairwise_Q20_dna_counts_retained(monkeypatch,tmp_path):
    dna,rows,original,ml,nj=small(monkeypatch,tmp_path)
    intervals=m.partitions_check(rows)
    names,cov,alt,pairs,err=m.source_pair_counts(dna,intervals,original)
    assert len(pairs)==6
    assert cov.shape==(3,6)
    assert alt.shape==(3,6)
    assert err==pytest.approx(0)
    assert int(cov.sum(axis=0).min())==12


def test_full_locus_100NJ_resamples_not_bases_and_not_orthology(monkeypatch,tmp_path):
    dna,rows,original,ml,nj=small(monkeypatch,tmp_path)
    names,cov,alt,pairs,_=m.source_pair_counts(dna,m.partitions_check(rows),original)
    split,report=m.bootstrap_analysis(names,cov,alt,pairs,ml,nj,n=100,seed=23)
    assert len(split)==1
    assert report["n_whole_locus_resampled_NJ_trees"]==100
    assert report["original_same_copy_Q20_loci"]==3
    assert report["original_masked_Q20_sites"]==12
    assert report["site_only_bootstrap_or_IQTREE_ML_bootstrap_claimed"] is False
    assert report["U_to_D_events_or_rainfall_adaptive_radiation_tested"] is False


def test_locus_partition_overlaps_are_unacceptable(monkeypatch,tmp_path):
    dna,rows,original,ml,nj=small(monkeypatch,tmp_path)
    rows[1]["start_1based"]="4"
    with pytest.raises(ValueError,match="whole-locus boundaries"):
        m.partitions_check(rows)


def test_changed_original_distance_pair_rejected(monkeypatch,tmp_path):
    dna,rows,original,ml,nj=small(monkeypatch,tmp_path)
    original[0]["n_differing_Q20_positions"]="99"
    with pytest.raises(ValueError,match="original Q20 distances"):
        m.source_pair_counts(dna,m.partitions_check(rows),original)
