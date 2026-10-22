"""Synthetic numerical parser tests only; real samples run in source-verified CI."""
import numpy as np
import pytest

from analysis.v4 import t2_moreyra_all125_author241_site_depth as m


def source_idxstats():
    return "\n".join([f"lett-L{i:04d}\t300\t10\t2" for i in range(1061)])+"\n*\t0\t0\t1000\n"


def test_1061_contig_identity_and_depth_position_order():
    lengths,loci=m.reference_lengths(source_idxstats())
    assert len(lengths)==1061 and len(loci)==1061
    base=["lett-L0000\t1\t5\n","lett-L0000\t2\t3\n",
          "lett-L0000\t3\t1\n","lett-L0001\t3\t4\n"]
    depth=m.count_depth_stream(base,lengths)
    assert depth["lett-L0000"]==[3,2,1]
    assert depth["lett-L0001"]==[1,1,0]
    with pytest.raises(ValueError,match="Duplicate/out-of-order"):
        m.count_depth_stream([base[0],base[0]],lengths)
    with pytest.raises(ValueError,match="Invalid aligned"):
        m.count_depth_stream(["lett-L0001\t999\t5\n"],lengths)


def test_best_one_ref_per_locus_does_not_add_multiple_copy_depth(monkeypatch):
    monkeypatch.setattr(m,"LOCUS_COUNT",3)
    lengths={f"{ref}-L{i}":200 for ref in ("lett","saff") for i in range(3)}
    refs={f"L{i}":[f"lett-L{i}",f"saff-L{i}"] for i in range(3)}
    depth={
      "lett-L0":[150,40,20],
      "saff-L0":[160,130,100],
      "lett-L1":[75,60,40],
      "saff-L1":[100,50,50],
      "saff-L2":[30,12,10],
    }
    rows=m.summarize_loci(["L0","L1","L2"],lengths,refs,depth)
    assert rows[0]["best_original_reference"]=="saff-L0"
    assert rows[0]["dp3_callable_sites"]==130
    assert rows[0]["dp3_breadth"]==.65
    assert rows[1]["best_original_reference"]=="lett-L1"
    assert rows[1]["dp3_at_least100bp"] is False
    result=m.per_run_metrics(rows)
    assert result["n_loci_dp3_ge50pct"]==1
    assert result["n_loci_dp3_100bp"]==1
    assert result["n_loci_dp3_any"]==3
    assert result["total_best_reference_sites_dp3"]==202
    assert result["n_loci_dp3_ge25pct"]==2
    assert result["n_loci_dp3_ge10pct"]==2
    assert not result.get("phylogeny_created",False)


def test_different_original_locus_set_not_allowed(tmp_path):
    f=tmp_path/"fake_241.txt"
    f.write_text("FAKE\n")
    with pytest.raises(ValueError,match="digest drift"):
        m.loci_from_pinned_author_list(f)


def test_unsampled_source_run_remains_in_125_source_queue(monkeypatch):
    monkeypatch.setattr(m,"RUN_COUNT",125)
    monkeypatch.setattr(m,"TAXA_COUNT",113)
    src=[{"run":f"SRR{i:08d}","biosample":f"SAMN{i:08d}",
          "accepted_key":f"K{i%113:03d}"} for i in range(125)]
    def mock(row,ref,loci,scratch):
        return {"run":row["run"],"biosample":row["biosample"],
                "status":"real_ENA_reference_site_depth_q20_Q20",
                "n_loci_dp3_ge25pct":3,"n_loci_dp3_ge50pct":1}
    output=m.collect(src,["L1"],None,None,4,10,runner=mock)
    assert output["requested"]==13
    assert output["successful"]==13
    assert output["failures"]==0
    assert output["read_data_in_artifact"] is False


def test_malformed_quality_threshold_coverage_stops(monkeypatch):
    monkeypatch.setattr(m,"LOCUS_COUNT",1)
    with pytest.raises(ValueError,match="Reference coverage exceeds"):
        m.summarize_loci(["L0"],{"lett-L0":150},{"L0":["lett-L0"]},
                         {"lett-L0":[151,50,40]})
