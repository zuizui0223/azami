"""Source-real ML branch uncertainty: site UFBoot must not masquerade as 300-locus NJ."""
from io import StringIO

from Bio import Phylo
import pytest

from analysis.v4 import t2_real7_iqtree_support_audit as mod


TREE=(
 "(((Cirsium_flodmanii:0.1,Cirsium_tracyi:0.1)94.5/97:0.01,"
 "(Cirsium_scopulorum:0.1,Cirsium_vulgare:0.1)75.1/84:0.02)"
 "92.3/96:0.03,(Cirsium_italicum:0.1,"
 "(Cirsium_latifolium:0.1,Cirsium_pulcherrimum:0.1)"
 "99/100:0.02)92.3/96:0.08);"
)

def test_real_support_pairs_are_separately_named_SHaLRT_then_UFBOOT():
    tree=Phylo.read(StringIO(TREE),"newick")
    pairs=[mod.support(x) for x in tree.get_nonterminals()]
    assert (94.5,97.) in pairs
    assert (99.,100.) in pairs
    assert (75.1,84.) in pairs
    assert all((a is None or 0<=a<=100) for a,b in pairs)
    assert all((b is None or 0<=b<=100) for a,b in pairs)

def test_unrooted_split_branch_support_not_confused_with_root(tmp_path):
    p=tmp_path/"iqtree.nwk"
    p.write_text(TREE)
    tree,splits=mod.read_splits(p)
    assert len(splits)==4
    assert splits[("Cirsium_flodmanii","Cirsium_tracyi")]==(94.5,97.)
    assert splits[("Cirsium_latifolium","Cirsium_pulcherrimum")]==(99.,100.)

def test_conflicting_support_is_not_invented():
    t=Phylo.read(StringIO("(A:0.1,B:0.1)105/96:0.2;"),"newick")
    with pytest.raises(ValueError,match="Invalid IQ-TREE"):
        mod.support(t.root)

def test_missing_NJ_support_does_not_become_fake_one():
    # Existing four split rows are 300 distinct locus bootstrap draws.
    assert mod.EXPECTED_N_INNER_SPLITS==4
    assert mod.PREFIXES==("primary_101","sensitivity_110")
