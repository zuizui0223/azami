"""Synthetic source identity tests. Empirical workflow also builds exact EAzami294."""
from copy import deepcopy
import pytest
from analysis.v4 import t2_eazami294_s1_run_bridge as mod

def example(monkeypatch):
    monkeypatch.setattr(mod,"EXPECTED",{
       "s1_taxa":3,"s1_photos":7,"eazami_samples":3,
       "eazami_runs":4,"eazami_moreyra":2,
       "eazami_chang":1,"s1_sra_rows":2,"s1_sra_taxa":1})
    s1=[
        {"accepted_key":"a","source_names":"Cirsium_a",
         "source_concept_collision":"False",
         "normalized_exact_source_name_candidate":"Cirsium a",
         "eligible_photo_observations":"2"},
        {"accepted_key":"b","source_names":"Cirsium_b",
         "source_concept_collision":"False",
         "normalized_exact_source_name_candidate":"Cirsium b",
         "eligible_photo_observations":"3"},
        {"accepted_key":"c","source_names":"Cirsium_c|Cirsium_d",
         "source_concept_collision":"True",
         "normalized_exact_source_name_candidate":"",
         "eligible_photo_observations":"2"}
    ]
    sra=[{
        "accepted_key":"a","S1_source_taxon":"Cirsium_a",
        "original_tree_code":"Cirsium a",
        "biosample":"SAMN1","run":run,"voucher_and_herbarium":"Original voucher"
    } for run in ("SRR1","SRR2")]
    panel=[
       {"panel_id":"MRY_SAMN1","biosample":"SAMN1",
        "run_accessions":"SRR1|SRR2","run_count":"2",
        "source_studies":"Moreyra2025","source_taxon_label":"Cirsium a",
        "analysis_taxon_label":"Cirsium a","assay":"Compositae1061_target_capture",
        "name_review_required":"false"},
       {"panel_id":"MRY_SAMN2","biosample":"SAMN2",
        "run_accessions":"SRR3","run_count":"1",
        "source_studies":"Moreyra2025","source_taxon_label":"Cirsium b",
        "analysis_taxon_label":"Cirsium b","assay":"Compositae1061_target_capture",
        "name_review_required":"false"},
       {"panel_id":"CH_SAMN3","biosample":"SAMN3",
        "run_accessions":"SRR4","run_count":"1",
        "source_studies":"Chang2026","source_taxon_label":"Cirsium b",
        "analysis_taxon_label":"Cirsium b","assay":"leaf_RNAseq_transcriptome",
        "name_review_required":"false"},
    ]
    return s1,sra,panel

def test_two_runs_of_one_biosample_not_duplicate_nuclear_taxa(monkeypatch):
    s1,sra,panel=example(monkeypatch)
    full,perrun,report=mod.link(s1,sra,panel)
    assert len(full)==3 and len(perrun)==2
    assert report["source_run_identity_status"]=={
       "exact_run_BioSample_tree_code":2}
    assert report["taxon_level_status"]["confirmed_EAzami294_same_run_taxa"]==1
    assert report["taxon_level_status"]["Chang_label_only_taxa"]==1
    assert report["EAzami_original_biological_samples"]==3
    assert report["EAzami_original_run_accessions"]==4
    assert report["EAzami_294_global_nuclear_tree_completed"] is False
    assert all(x["accepted_concept_voucher_identity_verified"] is False for x in perrun)
    assert full[1]["n_EAzami294_Chang_label_only_samples"]==1
    assert full[2]["WCVP_source_alias_collision"] is True
    assert full[2]["n_EAzami294_Moreyra_label_only_samples"]==0

def test_biosample_same_run_disagreement_fails(monkeypatch):
    s1,sra,panel=example(monkeypatch)
    sra[0]["biosample"]="SAMN_OTHER"
    with pytest.raises(ValueError,match="BioSample"):
        mod.link(s1,sra,panel)

def test_tree_code_name_disagreement_fails(monkeypatch):
    s1,sra,panel=example(monkeypatch)
    sra[0]["original_tree_code"]="Cirsium incorrectly_labeled"
    with pytest.raises(ValueError,match="tree_code"):
        mod.link(s1,sra,panel)

def test_chang_duplicate_run_in_294_is_unacceptable(monkeypatch):
    s1,sra,panel=example(monkeypatch)
    panel[2]["run_accessions"]="SRR3"
    with pytest.raises(ValueError,match="duplicates source run"):
        mod.link(s1,sra,panel)

def test_real_contract_targets_frozen_census_and_not_global_tree():
    assert mod.EA_PANEL_FIELDS>={
        "biosample","run_accessions","analysis_taxon_label"}
    assert mod.normalize_name("Cirsium_vulgare")=="Cirsium vulgare"
