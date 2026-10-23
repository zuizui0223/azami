"""Real 131 GenBank linked-rDNA T3 tree node audit and biological stop rules."""
import io
from Bio import Phylo
from analysis.v4 import t3_2025_real_nuclear_rDNA_phylogeny_posture_topology as m


def test_root_independent_original_genus_exclusivity_from_actual_unrooted_splits():
    tree=Phylo.read(io.StringIO("((Cremanthodium_A:1,Cremanthodium_B:1):1,"
                "((Parasenecio_X:1,Parasenecio_Y:1):1,(Ligularia_M:1,Ligularia_N:1):1):1);"),"newick")
    assert m.root_independent_split_is_exclusive(tree,"Cremanthodium")
    assert m.root_independent_split_is_exclusive(tree,"Parasenecio")
    assert m.root_independent_split_is_exclusive(tree,"Ligularia")
    tree.root_at_midpoint()
    assert m.root_independent_split_is_exclusive(tree,"Cremanthodium")
    assert m.root_independent_split_is_exclusive(tree,"Parasenecio")


def test_genus_mixing_breaks_root_invariant_genus_split_without_falsifying_morphology():
    tree=Phylo.read(io.StringIO("((Cremanthodium_A:1,Parasenecio_X:1):1,"
                "((Parasenecio_Y:1,Ligularia_M:1):1,(Ligularia_N:1,Cremanthodium_B:1):1):1);"),"newick")
    assert not m.root_independent_split_is_exclusive(tree,"Cremanthodium")
    assert not m.root_independent_split_is_exclusive(tree,"Parasenecio")


def test_actual_published_genealogy_provenance_boundaries():
    assert m.N_GENUINE_GB_ITS_ETS==131
    assert m.D_N==6 and m.UNKNOWN_N==125
    assert m.STRONG_ALRT==80.0 and m.STRONG_UFBOOT==95.0
    assert m.TREE_SHA=="fb1c2019b06e3691fac13e19d2f250ad7a34287fbd29abd209b9405d30748310"
    assert m.TREE_RUN==39961377067
    assert "NOT a published whole-nuclear coalescent species tree" in m.__doc__
    assert "NO independent" in m.__doc__
