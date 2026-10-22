"""IQ-TREE observed-only seven-source ML: no imputation or species-tree promotion."""
from io import StringIO

import pytest
from Bio import Phylo

from analysis.v4 import t2_real7_author101_iqtree_ml as m


NAMES=[
    "Cirsium_flodmanii","Cirsium_scopulorum","Cirsium_tracyi",
    "Cirsium_vulgare","Cirsium_italicum","Cirsium_latifolium",
    "Cirsium_pulcherrimum"
]

def test_multiple_observed_nuclear_loci_as_seven_complete_case_only():
    loci=[
        {"locus":"L1","alignment":["ACGT" for _ in NAMES]},
        {"locus":"L2","alignment":["ATGC" if k%2 else "AGGC" for k in range(7)]},
    ]
    s=m.alignment_by_samples(NAMES,loci)
    assert len(s)==7 and all(len(x)==8 for x in s.values())
    assert s[NAMES[0]]=="ACGTAGGC"
    assert s[NAMES[1]]=="ACGTATGC"
    payload=m.fasta_payload(s)
    assert payload.count(">")==7
    assert payload.endswith("\n")

def test_bad_base_without_source_Q20_is_not_used_as_ancestral_imputation():
    with pytest.raises(ValueError,match="unobserved or reference-imputed"):
        m.alignment_by_samples(NAMES,[{
            "locus":"L1","alignment":["AAAA"]*6+["AANA"]
        }])

def test_reference_source_names_are_never_inferred_taxa():
    bad=NAMES.copy()
    bad[0]="OUTGROUP_saff"
    with pytest.raises(ValueError,match="Unverified author"):
        m.alignment_by_samples(bad,[{
            "locus":"L1","alignment":["AAAA"]*7
        }])

def test_ml_unrooted_internal_splits_and_known_quartet_preserved():
    # Explicit known F+T split against S+V; other 3 tips arbitrary.
    tree=Phylo.read(StringIO(
        "(((Cirsium_flodmanii:0.1,Cirsium_tracyi:0.1):0.1,"
        "(Cirsium_scopulorum:0.2,Cirsium_vulgare:0.2):0.1):0.1,"
        "(Cirsium_italicum:0.1,(Cirsium_latifolium:0.1,"
        "Cirsium_pulcherrimum:0.1):0.1):0.2);"),"newick")
    splits=m.splits_with_support(tree,7)
    assert len(splits)==4
    assert tuple(sorted(("Cirsium_flodmanii","Cirsium_tracyi"))) in splits
    assert m.quartet_split_four(tree)==[
        "Cirsium_flodmanii","Cirsium_tracyi"]

def test_absent_or_unresolved_technical_species_split_fails():
    tree=Phylo.read(StringIO(
       "(" +",".join(x+":0.1" for x in NAMES)+");"),"newick")
    with pytest.raises(ValueError,match="fully binary"):
        m.splits_with_support(tree,7)

def test_technical_IQtree_model_and_dating_are_prespecified():
    assert m.MODEL=="GTR+G4"
    assert m.REPS==1000
    assert m.SEED==20261023
    assert m.EXPECTED_101==101
    assert m.EXPECTED_101_SITES==28726
    assert m.EXPECTED_110_SITES==29270
