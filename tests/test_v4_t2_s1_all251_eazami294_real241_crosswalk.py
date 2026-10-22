"""No arbitrary collapsing of the original 251 S1 source names into 244 species."""
import pytest
from analysis.v4 import t2_s1_all251_eazami294_real241_crosswalk as m


def fixture(monkeypatch):
    for key,value in (("N_SOURCE",3),("N_ACCEPTED",2),("N_PHOTOS",10),
                      ("N_MERGED_KEYS",1),("N_Q20_TAXA",1),
                      ("N_EAZAMI_BIOSAMPLES",3),("N_EAZAMI_RUNS",3)):
        monkeypatch.setattr(m,key,value)
    rows=[
        {"accepted_key":"K0","original_source_name":"Cirsium_other",
         "original_photo_count":"4","mean_image_angle_deg":"60",
         "chelsa_bio12":"1000"},
        {"accepted_key":"206745753","original_source_name":"Cirsium_brevicaule",
         "original_photo_count":"5","mean_image_angle_deg":"75",
         "chelsa_bio12":"2500"},
        {"accepted_key":"206745753","original_source_name":"Cirsium_irumtiense",
         "original_photo_count":"1","mean_image_angle_deg":"80",
         "chelsa_bio12":"3000"},
    ]
    depth=[{"accepted_key":"K0","source_taxon":"Cirsium_other",
            "maximum_per_run_loci_dp3_ge50pct":"150",
            "maximum_per_run_loci_dp3_ge25pct":"210"}]
    panel=[
       {"tip_id":"JOG0001","biosample":"SAMN1","run_accessions":"SRR1",
        "source_study":"Moreyra2025","analysis_taxon_label":"Cirsium other"},
       {"tip_id":"JOG0002","biosample":"SAMN2","run_accessions":"SRR2",
        "source_study":"Chang2026","analysis_taxon_label":"Cirsium brevicaule"},
       {"tip_id":"JOG0003","biosample":"SAMN3","run_accessions":"SRR3",
        "source_study":"Chang2026","analysis_taxon_label":"Cirsium irumtiense"},
    ]
    return rows,depth,panel


def test_exact_original_source_names_keep_brevicaule_versus_irumtiense(monkeypatch):
    rows,depth,panel=fixture(monkeypatch)
    result,report=m.build(rows,depth,panel)
    assert len(result)==3 and report["total_original_S1_WCVP_accepted_keys"]==2
    assert report["original_S1_units_with_real_241_depth"]==1
    assert result[0]["Q20_dp3_ge50pct_loci_best_single_run"]==150
    assert result[1]["Q20_dp3_ge50pct_loci_best_single_run"] is None
    assert result[2]["Q20_dp3_ge50pct_loci_best_single_run"] is None
    assert report["C_brevicaule_and_C_irumtiense_separate_source_name_units"]=={
        "Cirsium_brevicaule":{"n_photo":5,"n_eazami294_name_matched_samples":1,
                              "n_eazami294_Chang_samples":1},
        "Cirsium_irumtiense":{"n_photo":1,"n_eazami294_name_matched_samples":1,
                              "n_eazami294_Chang_samples":1},
    }
    assert not report["EAzami_global_nuclear_tree_reconstructed"]


def test_collapsed_wcvp_species_not_pooled_with_single_other_specimen(monkeypatch):
    rows,depth,panel=fixture(monkeypatch)
    depth[0]["accepted_key"]="206745753"
    with pytest.raises(ValueError,match="multi-name collision"):
        m.build(rows,depth,panel)


def test_wrong_biosample_count_fails(monkeypatch):
    rows,depth,panel=fixture(monkeypatch)
    with pytest.raises(ValueError,match="294 deduplicated"):
        m.build(rows,depth,panel[:2])
