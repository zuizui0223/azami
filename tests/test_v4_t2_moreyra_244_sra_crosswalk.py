"""Offline checks: full original S1 source and Moreyra public SRA evidence gates."""
from collections import defaultdict
import pytest

import analysis.v4.t2_moreyra_244_sra_crosswalk as r


def small_census():
    return [
       dict(accepted_key="a",normalized_exact_source_name_candidate="Cirsium arvense",
            source_names="Cirsium_arvense",source_concept_collision="False",
            eligible_photo_observations="20",Moreyra2025_SRA_taxon_name_candidate="True"),
       dict(accepted_key="b",normalized_exact_source_name_candidate="Cirsium vulgare",
            source_names="Cirsium_vulgare",source_concept_collision="False",
            eligible_photo_observations="30",Moreyra2025_SRA_taxon_name_candidate="True"),
       dict(accepted_key="c",normalized_exact_source_name_candidate="Cirsium funkiae",
            source_names="Cirsium_funkiae",source_concept_collision="False",
            eligible_photo_observations="40",Moreyra2025_SRA_taxon_name_candidate="True"),
       dict(accepted_key="d",normalized_exact_source_name_candidate="",
            source_names="Cirsium_brevicaule|Cirsium_irumtiense",
            source_concept_collision="True",eligible_photo_observations="80",
            Moreyra2025_SRA_taxon_name_candidate="False"),
    ]


def public_sample(name,run,match_status="exact"):
    return {
        "tree_code":name, "biosample":"SAMN"+run[3:],
        "run":run,"experiment":"SRX000",
        "voucher_and_herbarium":"public_author_voucher",
        "scope_class":"outside_scope",
        "name_reconciliation_priority":"high" if match_status!="exact" else "low",
        "tree_code_vs_sra_name":match_status,
    }


def test_public_name_matching_requires_same_244_flag_and_two_s2(monkeypatch):
    monkeypatch.setattr(r,"EXPECTED_MATCH_N",3)
    source=defaultdict(list,{
       "Cirsium arvense":[public_sample("Cirsium arvense","SRR111")],
       "Cirsium vulgare":[public_sample("Cirsium vulgare 1","SRR222")],
       "Cirsium funkiae":[public_sample("Cirsium funkiae","SRR333",
                                       "different_submitted_or_published_name")],
       "Cirsium brevicaule":[public_sample("Cirsium brevicaule","SRR444")],
    })
    out,per_run,report=r.reconcile(small_census(),source)
    assert len(out)==4 and len(per_run)==3
    assert report["S1_with_Moreyra_public_SRA_name_candidate"]==3
    assert report["candidate_records_where_submitted_name_differs"]==1
    assert report["S2_focal_both_have_SRA_sample_candidates"] is True
    assert report["original_350_locus_Moreyra_tree_recovered"] is False
    assert report["accepted_taxonomic_concept_and_same_voucher_verified"] is False
    assert out[3]["reconstruction_status"]=="WCVP_source_name_collision"
    assert out[3]["n_public_SRA_samples"]==0
    assert all(row["accepted_concept_verified"]=="no" for row in per_run)


def test_unverified_name_without_run_does_not_silently_pass(monkeypatch):
    monkeypatch.setattr(r,"EXPECTED_MATCH_N",3)
    census=small_census()
    census[2]["Moreyra2025_SRA_taxon_name_candidate"]="False"
    sra={
       "Cirsium arvense":[public_sample("Cirsium arvense","SRR111")],
       "Cirsium vulgare":[public_sample("Cirsium vulgare","SRR222")],
       "Cirsium funkiae":[public_sample("Cirsium funkiae","SRR333")],
    }
    with pytest.raises(ValueError,match="matching drifted"):
        r.reconcile(census,sra)


def test_taxon_rank_is_not_collapsed_into_species_name():
    from analysis.v4.t2_native_s1_all_taxa_registry import generic_taxon_label
    assert generic_taxon_label("Cirsium vulgare 3")=="Cirsium vulgare"
    assert generic_taxon_label("Cirsium japonicum var. albescens")!=\
           generic_taxon_label("Cirsium japonicum var. fukienense")
