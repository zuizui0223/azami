"""Tests: empirical PGLS placement screen never replaces independent T2 phylogeny."""
from io import StringIO
from pathlib import Path

from Bio import Phylo
import pytest

from analysis.v4 import t2_pgls_placement_screen as screen


def author_fixture():
    rows=["Species,HeadPosition"]
    for i in range(63):
        rows.append(f"Cirsium_{i:02d},{int(i >= 58)}")
    return ("\n".join(rows)+"\n").encode()


def test_author_state_git_blob_pin_and_codes(monkeypatch):
    raw=author_fixture()
    monkeypatch.setattr(screen, "SOURCE_GIT_BLOB", screen.git_blob_sha(raw))
    states=screen.read_author_states(raw)
    assert len(states)==63
    assert sum(states.values())==5
    bad=raw.replace(b"Cirsium_00",b"Cirsium_01")
    with pytest.raises(ValueError, match="Git blob"):
        screen.read_author_states(bad)


def test_author_state_duplicates_and_wrong_counts_guard(monkeypatch):
    raw=author_fixture()
    bad=raw.replace(b"Cirsium_00",b"Cirsium_01")
    monkeypatch.setattr(screen, "SOURCE_GIT_BLOB", screen.git_blob_sha(bad))
    with pytest.raises(ValueError,match="Duplicate taxon"):
        screen.read_author_states(bad)
    fewer=raw.replace(b"Cirsium_62,1",b"Cirsium_62,0")
    monkeypatch.setattr(screen, "SOURCE_GIT_BLOB", screen.git_blob_sha(fewer))
    with pytest.raises(ValueError,match="nodding-state"):
        screen.read_author_states(fewer)


def tree(names):
    return Phylo.read(StringIO("("+",".join(f"{n}:0.2" for n in names)+");"),"newick")


def test_minimum_transitions_is_not_an_evolutionary_origin_count():
    names=["Cirsium_00","Cirsium_01","Cirsium_58","Cirsium_59"]
    tr=tree(names)
    scores={name:int(i>=2) for i,name in enumerate(names)}
    result=screen.parsimony(tr,scores)
    assert result["unknown_root_min_steps"]==2
    assert result["root_erect_min_steps"]==2
    assert result["root_nodding_min_steps"]==2


def test_52_scenario_conditional_reuse_synthetic(tmp_path):
    states={f"Cirsium_{i:02d}":int(i>=58) for i in range(63)}
    names=list(states)[:10]+list(states)[-5:]
    trees=[(f"test_scenario_{i}",tree(names)) for i in range(52)]
    output=tmp_path/"out.json"
    report=screen.evaluate_placements(trees,states,"synthetic_sha",output,draws=3,min_scored=12)
    result=report["summary"]
    assert result["n_scenario_results"]==52
    assert result["matched_taxa_n_range"]==[15,15]
    assert result["matched_nodding_n_range"]==[5,5]
    assert result["status_counts"]["conditional_placement_Mk_sensitivity_NOT_independent_history"]==52
    assert result["rain_temperature_evolutionary_test_executed"] is False
    assert result["empirical_transition_dates_estimated"] is False
    assert "unknown_root_min_steps_scenario_range" in result
    assert output.exists()


def test_no_promotion_when_too_few_published_nodding_matches(tmp_path):
    states={f"Cirsium_{i:02d}":int(i>=58) for i in range(63)}
    names=list(states)[:13]+["Cirsium_62"]
    trees=[(f"test_scenario_{i}",tree(names)) for i in range(52)]
    result=screen.evaluate_placements(trees,states,"synthetic_sha",tmp_path/"blocked.json",
                                     min_scored=12)
    assert result["summary"]["matched_nodding_n_range"] == [1,1]
    assert result["summary"]["status_counts"]["descriptive_parsimony_only_insufficient_binary_coverage"]==52
    assert "conditional_placement_Mk_sensitivity_NOT_independent_history" not in result["summary"]["status_counts"]


def test_tree_validation_catches_branch_and_tip_errors():
    tr=tree(["Cirsium_00","Cirsium_01"])
    assert screen.check_placement_tree(tr)=={"Cirsium_00","Cirsium_01"}
    tr.get_terminals()[1].name="Cirsium_00"
    with pytest.raises(ValueError,match="not unique"):
        screen.check_placement_tree(tr)
    tr.get_terminals()[1].name="Cirsium_01"
    tr.get_terminals()[0].branch_length=-1
    with pytest.raises(ValueError,match="negative"):
        screen.check_placement_tree(tr)


def test_published_nuclear_tree_blob_and_coverage_guard(monkeypatch):
    states={f"Cirsium_{i:02d}":int(i>=58) for i in range(63)}
    raw=("("+",".join(f"{n}:0.2" for n in states)+");\n").encode()
    monkeypatch.setattr(screen,"AUTHOR_NUCLEAR_TREE_GIT_BLOB",screen.git_blob_sha(raw))
    phy=screen.read_author_nuclear_tree(raw,states)
    assert len(phy.get_terminals())==63
    with pytest.raises(ValueError,match="Git blob"):
        screen.read_author_nuclear_tree(raw.replace(b"Cirsium_00",b"Cirsium_xx"),states)
    missing=states.copy()
    missing.pop("Cirsium_00")
    with pytest.raises(ValueError,match="tip set"):
        screen.read_author_nuclear_tree(raw,missing)


def test_independent_nuclear_same_44_tip_comparison_synthetic(tmp_path):
    states={f"Cirsium_{i:02d}":int(i>=58) for i in range(63)}
    original=tree(list(states))
    names=list(states)[:10]+list(states)[-5:]
    placements=[(f"test_scenario_{i}",tree(names)) for i in range(52)]
    out=screen.evaluate_placements(placements,states,"test_sha",tmp_path/"x.json",
                                   published_nuclear_tree=original,nuclear_tree_sha="source-tree-sha",
                                   draws=3,min_scored=12)
    reference=out["summary"]["independent_nuclear_reference"]
    assert reference["common_tip_comparison"]["n_taxa"]==15
    assert reference["common_tip_comparison"]["n_nodding"]==5
    assert reference["complete_published_63_tip_context"]["n_taxa"]==63
    assert reference["not_global_cirsium_tree"] is True
    assert reference["common_tip_comparison"]["chronogram_absolute_dates_inferred"] is False
    with pytest.raises(ValueError,match="original SHA256"):
        screen.evaluate_placements(placements,states,"test_sha",tmp_path/"y.json",
                                   published_nuclear_tree=original,draws=2,min_scored=12)
