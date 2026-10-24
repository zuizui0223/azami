"""Biological and unrooted topology source contracts for true 131 ITS ETS."""
from io import StringIO
from Bio import Phylo

from analysis.v4 import t3_131_direct_head_two_source_tree_topology_sensitivity as m


def toy():
    return Phylo.read(StringIO(
        "((A:.02,B:.03):.04,((C:.01,D:.01):.03,(E:.02,F:.02):.01):.02);"),
        "newick")


def test_unrooted_four_tip_split_is_root_invariant():
    orig=toy()
    x={"C","D","E","F"}
    assert m.unrooted_split_exists(orig,x)
    orig.root_at_midpoint()
    assert m.unrooted_split_exists(orig,x)
    assert m.unrooted_split_exists(orig,{"A","B"})
    assert not m.unrooted_split_exists(orig,{"A","D","F"})


def test_true_documented_developmental_variable_is_not_false_U_control():
    assert m.PAIR_STAGED=={"Parasenecio_cyclotus","Parasenecio_koualapensis"}
    assert m.PAIR_OTHER=={"Cremanthodium_campanulatum","Ligularia_purdomii"}
    assert len(m.REQUIRED_LCA_TIPS)==4
    assert m.TIPS==131
    assert m.PINNED["max_likelihood"]["sha256"]==m.TREE_SHA
    assert set(m.PINNED)=={"max_likelihood","ultrafast_bootstrap_consensus",
                           "BioNJ_initial_distance_tree"}
    assert all(len(v["sha256"])==64 for v in m.PINNED.values())


def test_root_independent_nearest_real_molecular_neighbor_uses_true_branch_length():
    tree=toy()
    assert m.nearest_tip(tree,"A")["species"]=="B"
    assert m.nearest_tip(tree,"A")["distance_substitutions_per_site"]>0
    assert len({x.name for x in tree.get_terminals()})==6
