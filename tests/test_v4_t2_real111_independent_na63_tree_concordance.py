"""Synthetic independent exact-name nuclear topology and quartet comparison tests."""
from io import StringIO
import pytest
from Bio import Phylo

from analysis.v4 import t2_real111_independent_na63_tree_concordance as m


def original(monkeypatch):
    monkeypatch.setattr(m,"NA63",6)
    monkeypatch.setattr(m,"SAMPLE_UNIVERSE",6)
    tips=[f"Cirsium_{c}" for c in "abcdef"]
    sample=[{"accepted_key":str(1000+i),
             "taxon_name":"Cirsium "+letter,
             "phylogenetic_tip":"WCVP"+str(1000+i),
             "run":f"SRR200{i}"}
            for i,letter in enumerate("abcdef")]
    # Rooting of an unrooted six-tip topology is incidental, all edges positive.
    pub="(((Cirsium_a:1,Cirsium_b:1):2,Cirsium_c:3):1,((Cirsium_d:1,Cirsium_e:1):2,Cirsium_f:3):1);"
    tech="(((WCVP1000:1,WCVP1001:1):2,WCVP1002:3):1,((WCVP1003:1,WCVP1004:1):2,WCVP1005:3):1);"
    return sample,Phylo.read(StringIO(pub),"newick"),Phylo.read(StringIO(tech),"newick")


def test_independent_matching_topologies_no_bm_ou_or_adaptation_claim(monkeypatch):
    s,p,t=original(monkeypatch)
    result=m.compare(t,p,s)
    assert result["exact_original_taxon_overlap"]==6
    assert result["normalized_unrooted_RF"]==0.0
    assert result["n_resolved_discordant_quartets"]==0
    assert result["resolved_quartet_concordance_fraction"]==1.0
    assert result["n_possible_quartets"]==15
    assert not result["adaptive_radiation_or_environmental_causality_tested"]
    assert result["real111_technical_tree_is_not_single_copy_orthology_species_tree"]


def test_true_topology_conflict_must_not_be_labelled_independent_match(monkeypatch):
    s,p,t=original(monkeypatch)
    new=Phylo.read(StringIO(
        "(((WCVP1000:1,WCVP1002:1):2,WCVP1001:3):1,"
        "((WCVP1003:1,WCVP1004:1):2,WCVP1005:3):1);"
    ),"newick")
    result=m.compare(new,p,s)
    assert result["normalized_unrooted_RF"]>0
    assert result["n_resolved_discordant_quartets"]>0
    assert result["resolved_quartet_concordance_fraction"]<1.0


def test_taxon_name_is_not_auto_collapsed_to_accepted_species(monkeypatch):
    s,p,t=original(monkeypatch)
    s[2]["taxon_name"]="Cirsium a"
    with pytest.raises(ValueError,match="Duplicate original"):
        m.compare(t,p,s)


def test_insufficient_original_taxon_overlap_returns_stop_not_adaptation(monkeypatch):
    s,p,t=original(monkeypatch)
    m.MIN_COMMON=5
    limited="((Cirsium_a:1,Cirsium_b:1):1,(Cirsium_c:1,(Cirsium_d:1,Cirsium_Z:1):1):1);"
    published=Phylo.read(StringIO(limited),"newick")
    m.NA63=5
    result=m.compare(t,published,s)
    assert result["exact_original_taxon_overlap"]==4
    assert not result["technical_unrooted_NA63_validation_completed"]
    assert not result["exact_taxonomic_voucher_identity_verified"]


def test_quartet_rejects_exact_polytomy_not_new_history():
    distances={
        ("a","b"):2,("a","c"):2,("a","d"):2,
        ("b","c"):2,("b","d"):2,("c","d"):2
    }
    assert m.quartet_choice(distances,("a","b","c","d")) is None
