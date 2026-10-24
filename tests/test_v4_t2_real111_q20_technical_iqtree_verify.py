"""Synthetic source-constrained ML validation: never accept fake global species tree."""
import csv
import json
import pytest
from analysis.v4 import t2_real111_q20_technical_iqtree_verify as m


def source_files(tmp_path):
    src=tmp_path/"manifest.json"
    src.write_text(json.dumps({
        "source_111_independent_public_SRA_BioSamples":4,
        "preliminary_GTRG4_IQTREE_permitted":True,
        "single_copy_and_homeolog_orthology_verified":False,
        "global_EAzami294_species_tree_accepted":False,
        "adaptive_radiation_or_rain_transition_supported":False,
        "all111_selected_original_reference_locus_partitions":2,
        "all111_source_Q20_80pct_reference_positions":400,
    }))
    samples=tmp_path/"samples.csv"
    with samples.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["phylogenetic_tip","taxon_name"])
        w.writeheader()
        w.writerows([{"phylogenetic_tip":f"WCVP{i}","taxon_name":f"Cirsium sample{i}"}
                     for i in range(4)])
    tree=tmp_path/"tree.nwk"
    tree.write_text("((WCVP0:0.1,WCVP1:0.2):0.1,(WCVP2:0.2,WCVP3:0.3):0.1);\n")
    return src,samples,tree


def test_four_reconstructed_real_test_tips_preserved(tmp_path,monkeypatch):
    monkeypatch.setattr(m,"N_TIPS",4)
    source,samples,tree=source_files(tmp_path)
    result=m.verify(source,samples,tree,tmp_path/"out.json")
    assert result["real_original_SRA_BioSamples"]==4
    assert result["real_source_GTRG4_tree_tips"]==4
    assert not result["nuclear_copy_identity_and_paralog_homeolog_review_completed"]
    assert not result["EAzami294_independent_species_tree_completed"]
    assert not result["BIO12_BIO15_ancient_rain_adaptation_inferred"]


def test_incomplete_actual_accepted_key_tip_set_fails(tmp_path,monkeypatch):
    monkeypatch.setattr(m,"N_TIPS",4)
    source,samples,tree=source_files(tmp_path)
    tree.write_text("((WCVP0:0.1,WCVP1:0.2):0.1,(WCVP2:0.2,WCVP999:0.3):0.1);\n")
    with pytest.raises(ValueError,match="not identical"):
        m.verify(source,samples,tree,tmp_path/"out.json")


def test_no_invented_orthology_or_global_tree(tmp_path,monkeypatch):
    monkeypatch.setattr(m,"N_TIPS",4)
    source,samples,tree=source_files(tmp_path)
    a=json.loads(source.read_text())
    a["single_copy_and_homeolog_orthology_verified"]=True
    source.write_text(json.dumps(a))
    with pytest.raises(ValueError,match="not accepted"):
        m.verify(source,samples,tree,tmp_path/"out.json")


def test_full_nuclear_111_not_294_source_scope():
    assert m.N_TIPS==111
    assert m.N_ORIGINAL==244
    assert m.N_SOURCE_NAMES==251
