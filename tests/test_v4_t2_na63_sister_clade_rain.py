"""Offline tests of the source-pinned T2 sister-clade comparison contract."""
from __future__ import annotations

from io import StringIO

import pytest
from Bio import Phylo
import analysis.v4.t2_na63_sister_clade_rain as mod


def fixture():
    tree = Phylo.read(StringIO(
        "(((D1:1,E1:1):1,(D2:1,E2:1):1):1,"
        "((D3:1,E3:1):1,(D4:1,E4:1):1):1);"), "newick")
    states = {f"D{i}": 1 for i in range(1, 5)}
    states.update({f"E{i}": 0 for i in range(1, 5)})
    rain = {"D1":600,"E1":400,"D2":300,"E2":700,
            "D3":800,"E3":700,"D4":500,"E4":900}
    climate = {taxon: {"BIO12":value,"BIO1":0.0}
               for taxon,value in rain.items()}
    return tree,states,climate


def test_four_clade_signs_and_no_inferred_rain_adaptation(monkeypatch):
    monkeypatch.setattr(mod, "EXPECTED_N", 8)
    monkeypatch.setattr(mod, "EXPECTED_DOWN", 4)
    monkeypatch.setattr(mod, "parsimony",
                        lambda tree,states: {
                            "root_erect_min_steps":4,
                            "root_nodding_min_steps":6,
                            "unknown_root_min_steps":4})
    tree,states,climate = fixture()
    out = mod.sister_clade_contrasts(tree,states,climate)
    assert out["n_maximal_pure_nodding_tip_clades"] == 4
    assert out["n_nodding_clades_wetter_than_immediate_erect_sister"] == 2
    assert out["n_nodding_clades_drier_than_immediate_erect_sister"] == 2
    assert out["repeated_independent_gain_count_established"] is False
    assert out["rain_dependent_evolutionary_rate_test_performed"] is False
    assert sorted(int(row["nodding_minus_sister_erect_BIO12_mm"])
                  for row in out["comparisons"]) == [-400,-400,100,200]


def test_missing_source_species_blocks_inference():
    tree,states,climate = fixture()
    del climate["E4"]
    with pytest.raises(ValueError, match="contract changed"):
        mod.sister_clade_contrasts(tree,states,climate)


def test_source_pins_are_not_grafted_tree_or_imputed_photo_labels():
    assert mod.SOURCE_DOI == "10.1086/724310"
    assert mod.EXPECTED_N == 63
    assert mod.EXPECTED_DOWN == 5
    assert mod.EXPECTED_PURE_DOWN_CLADES == 4
    assert mod.AUTHOR_NUCLEAR_TREE_GIT_BLOB == "6417a38793427bc7c8562dc0e2cab59b74d8c316"
