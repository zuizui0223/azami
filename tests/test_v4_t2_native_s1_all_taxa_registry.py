"""Taxon concept guards for the complete (not preselected) S1 244-key roster."""
from analysis.v4.t2_native_s1_all_taxa_registry import generic_taxon_label


def test_author_suffix_is_not_a_new_taxon_and_infraspecific_rank_survives():
    assert generic_taxon_label("Cirsium alpicola Nakai")=="Cirsium alpicola"
    assert generic_taxon_label("Cirsium japonicum Fisch. ex DC. var. horridum Nakai")=="Cirsium japonicum var. horridum"
    assert generic_taxon_label("Cirsium_japonicum_var._albescens")=="Cirsium japonicum var. albescens"
    assert generic_taxon_label("Cirsium japonicum var. horridum")!="Cirsium japonicum"


def test_unknown_or_non_cirsium_source_never_fabricates_a_tip():
    assert generic_taxon_label("Cirsium")==""
    assert generic_taxon_label("Cirsium SP")==""
    assert generic_taxon_label("Carduus nutans")==""
