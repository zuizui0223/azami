"""Fixture tests for empirical seven-taxon NJ vs original NA63 nuclear quartet."""
import json
from io import StringIO
import pytest
from Bio import Phylo

from analysis.v4 import t2_real7_NA63_quartet_crosscheck as m


def tree(nested=True,published=False):
    if nested:
        t="((Cirsium_flodmanii:1,Cirsium_tracyi:1):1,(Cirsium_scopulorum:1,Other1:1):1,(Cirsium_vulgare:1,(Other2:1,Other3:1):1):1);"
    else:
        t="((Cirsium_flodmanii:1,Cirsium_scopulorum:1):1,(Cirsium_tracyi:1,Other1:1):1,(Cirsium_vulgare:1,(Other2:1,Other3:1):1):1);"
    if published:
        for label in ('Other1','Other2','Other3'):
            t=t.replace(label,'NA_'+label)
    return Phylo.read(StringIO(t),"newick")


def receipt():
    return dict(n_original_public_SRA_biosamples=7,
                n_primary_loci_shared_atleast_100_Q20_sites=101,
                source_shared_Q20_positions_all_targets=29270,
                n_complete_case_invariant_plus_variable_sites_primary=28726,
                nuclear_tree_preliminary_only=True,
                EAzami_294_tip_species_tree_constructed=False)


def bootstrap():
    return [
        {"unrooted_split_side":"Cirsium_flodmanii|Cirsium_tracyi",
         "technical_NJ_locus_resampling_fraction":".75",
         "n_locus_bootstrap_resamples":"300"},
    ]


def test_quartet_concordance_is_exact_and_not_species_tree(monkeypatch):
    monkeypatch.setattr(m,"ORIGINAL_NA63",7)
    same=tree()
    published=tree(published=True)
    result=m.evaluate(receipt(),same,bootstrap(),published)
    assert result["unrooted_four_tip_topology_match"] is True
    assert result["independent_author_nuclear_quartet_pair"]==[
        "Cirsium_flodmanii","Cirsium_tracyi"]
    assert result["seven_tip_branch_bootstrap_not_quartet_bootstrap"]
    assert result["seven_taxa_NJ_not_accepted_as_genus_species_tree"]
    assert result["no_BM_OU_or_rain_adaptive_radiation_result"]


def test_discordant_pilot_quartet_is_reported_without_false_positive():
    published=tree(nested=True,published=True)
    discordant=tree(nested=False)
    assert m.quartet_split(published)!=m.quartet_split(discordant)
    result=m.evaluate(receipt(),discordant,bootstrap(),published)
    assert result["unrooted_four_tip_topology_match"] is False


def test_exact_name_collision_and_taxon_synonym_not_conflated():
    t=tree()
    with pytest.raises(ValueError,match="Four exact"):
        m.quartet_split(t,("Cirsium_flodmanii","Cirsium_scopulorum",
                            "Cirsium_tracyi","Cirsium_pulcherrimum"))


def test_pinned_author_nuclear_tree_source_digest_required():
    with pytest.raises(ValueError,match="Git blob"):
        m.load_na63(fetch=lambda url:b"(Cirsium_flodmanii,Cirsium_tracyi);")


def test_wrong_original_site_count_blocks_independent_quartet():
    r=receipt();r["source_shared_Q20_positions_all_targets"]=12345
    with pytest.raises(ValueError,match="authentic"):
        m.evaluate(r,tree(),bootstrap(),tree())
