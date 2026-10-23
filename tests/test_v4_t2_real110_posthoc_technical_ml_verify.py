"""Synthetic provenance guards for IQ-TREE on the one-source-excluded 110 panel."""
import csv
import json
from pathlib import Path

import pytest

from analysis.v4 import t2_real110_posthoc_technical_ml_verify as mod


def fixture(tmp_path):
    root=Path(tmp_path)
    receipt={
      "posthoc_sensitivity_taxa":110,
      "original_all111_Q20_successful_real_SRA_BioSamples":111,
      "original111_technical_gate_passed":False,
      "original_111_source_108_retained_loci":108,
      "observed_Q20_reference_positions_reused_unchanged":29310,
      "original_111_predeclared_89_of111_quorum_unmodified":89,
      "original_111_reference_position_selection_recomputed_on110":False,
      "original_S1_WCVP_master_accepted_keys":244,
      "original_S1_source_name_operational_units":251,
      "source111_single_species_Q20_outlier":{"tip":"WCVP206740970","n_real_Q20_positions":119},
      "copy_number_paralogy_polyploidy_review_completed":False,
      "global_EAzami294_true_species_tree_accepted":False,
      "true_Mk_nodding_transition_or_rainfall_adaptation_supported":False,
      "technical_GTRG4_ML_source_sensitivity_permitted":True,
    }
    (root/"audit.json").write_text(json.dumps(receipt))
    rows=[
        {"source_tip":f"WCVP{208000000+i}",
         "source_BioSample":f"SAMN{i+8000000}",
         "original_run":f"SRR{i+8000000}",
         "n_actual_Q20_called":"25200",
         "taxon_name":("Cirsium conspicuum" if i==0 else
                       "Cirsium scopulorum" if i==1 else f"Cirsium synthetic{i:03d}")}
        for i in range(110)
    ]
    with (root/"samples.csv").open("w",newline="") as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (root/"tree.nwk").write_text("(" + ",".join(
        f"{r['source_tip']}:0.001" for r in rows)+");\n")
    (root/"iqtree.report").write_text("Model of substitution: GTR+G4\n")
    return root,rows


def test_110_validated_dna_tree_not_original111_success(tmp_path):
    path,_=fixture(tmp_path)
    report=mod.verify(path/"tree.nwk",path/"audit.json",path/"samples.csv",path/"iqtree.report")
    assert report["recovered_technical_ML_tree_tip_count"]==110
    assert report["ML_tree_positive_branch_count"]==110
    assert report["original_locus_selection_108_no_posthoc_refit"]
    assert not report["per_locus_single_copy_ortholog_and_homeolog_verified"]
    assert not report["precipitation_driven_adaptive_radiation_tested"]
    assert not report["EAzami294_full_nuclear_species_tree_finished"]


def test_single_sparse_original_tip_is_never_reintroduced(tmp_path):
    path,_=fixture(tmp_path)
    text=(path/"tree.nwk").read_text().replace("WCVP208000000","WCVP206740970")
    (path/"tree.nwk").write_text(text)
    with pytest.raises(ValueError,match="authentic 110"):
        mod.verify(path/"tree.nwk",path/"audit.json",path/"samples.csv",path/"iqtree.report")


def test_a_priori_111_failure_cannot_be_relabelled_as_success(tmp_path):
    path,_=fixture(tmp_path)
    data=json.loads((path/"audit.json").read_text())
    data["original111_technical_gate_passed"]=True
    (path/"audit.json").write_text(json.dumps(data))
    with pytest.raises(ValueError,match="strict posthoc110"):
        mod.verify(path/"tree.nwk",path/"audit.json",path/"samples.csv",path/"iqtree.report")
