"""Small exact-tree tests, no simulated phylogenies passed off as empirical output."""
import hashlib
import json

import pytest

from analysis.v4 import t2_real110_ml_vs_nj_same_q20_topology as mod


def setup_case(tmp_path):
    ids=[f"WCVP20800000{i}" for i in range(6)]
    r=tmp_path
    ml="(("+ids[0]+":.01,"+ids[1]+":.01):.01,("+ids[2]+":.01,"+ids[3]+":.01):.01,("+ids[4]+":.01,"+ids[5]+":.01):.01);\n"
    nj="(("+ids[0]+":.01,"+ids[2]+":.01):.01,("+ids[1]+":.01,"+ids[3]+":.01):.01,("+ids[4]+":.01,"+ids[5]+":.01):.01);\n"
    (r/"ml.nwk").write_text(ml)
    (r/"nj.nwk").write_text(nj)
    source={
      "posthoc_sensitivity_taxa":6,
      "original_all111_Q20_successful_real_SRA_BioSamples":111,
      "original111_technical_gate_passed":False,
      "original_111_source_108_retained_loci":108,
      "observed_Q20_reference_positions_reused_unchanged":29310,
      "original_111_predeclared_89_of111_quorum_unmodified":89,
      "original_111_reference_position_selection_recomputed_on110":False,
      "original_S1_WCVP_master_accepted_keys":244,
      "original_S1_source_name_operational_units":251,
      "source111_single_species_Q20_outlier":{"tip":"WCVP206740970"}
    }
    (r/"source.json").write_text(json.dumps(source))
    ml_info={
      "recovered_technical_ML_tree_tip_count":6,
      "original_observed_Q20_nucleotide_sites":29310,
      "precipitation_driven_adaptive_radiation_tested":False,
      "unrooted_substitutions_per_site_not_clock":True,
      "treefile_sha256":hashlib.sha256((r/"ml.nwk").read_bytes()).hexdigest(),
    }
    (r/"ml.json").write_text(json.dumps(ml_info))
    nj_info={
      "source_posthoc_Q20_BioSamples":6,
      "original_Q20_reference_positions":29310,
      "original_failed111_quality_gate_relabelled_as_success":False,
      "quartets_are_not_independent_replicates":True,
    }
    (r/"nj.json").write_text(json.dumps(nj_info))
    comparison={
      "exact_original_taxon_overlap":4,
      "real111_technical_tree_is_not_single_copy_orthology_species_tree":True,
      "resolved_quartet_concordance_fraction":.75,
      "normalized_unrooted_RF":.50,
    }
    (r/"na_ml.json").write_text(json.dumps(comparison))
    comparison["resolved_quartet_concordance_fraction"]=.5
    (r/"na_nj.json").write_text(json.dumps(comparison))
    return r


def run(r):
    return mod.source_checks(
       r/"ml.nwk",r/"nj.nwk",r/"source.json",
       r/"na_ml.json",r/"na_nj.json",r/"ml.json",r/"nj.json")


def test_source_sameDNA_tree_split_RF_is_not_na63_precision(monkeypatch,tmp_path):
    monkeypatch.setattr(mod,"TIPS",6)
    monkeypatch.setattr(mod,"INDEPENDENT_NA63_OVERLAP",4)
    s=run(setup_case(tmp_path))
    assert s["ML_unrooted_internal_splits"]==3
    assert s["NJ_unrooted_internal_splits"]==3
    assert s["ML_and_NJ_shared_unrooted_internal_splits"]==1
    assert s["ML_vs_NJ_symmetric_Robinson_Foulds"]==4
    assert s["ML_vs_NJ_normalized_unrooted_RF"]==pytest.approx(4/6)
    assert s["ML_vs_NA63_quartet_fraction"]==.75
    assert s["NJ_vs_NA63_quartet_fraction"]==.5
    assert s["original111_quality_source_passed"] is False
    assert not s["botanical_U_to_D_change_number_or_rain_causality_inferred"]


def test_unpinned_ml_source_artifact_or_changed_tree_must_fail(monkeypatch,tmp_path):
    monkeypatch.setattr(mod,"TIPS",6)
    monkeypatch.setattr(mod,"INDEPENDENT_NA63_OVERLAP",4)
    p=setup_case(tmp_path)
    (p/"ml.nwk").write_text((p/"ml.nwk").read_text().replace(":.01",":.03",1))
    with pytest.raises(ValueError,match="SHA differs"):
        run(p)


def test_wrong_original111_scientific_quality_gate_is_not_acceptance(monkeypatch,tmp_path):
    monkeypatch.setattr(mod,"TIPS",6)
    monkeypatch.setattr(mod,"INDEPENDENT_NA63_OVERLAP",4)
    p=setup_case(tmp_path)
    s=json.loads((p/"source.json").read_text())
    s["original111_technical_gate_passed"]=True
    (p/"source.json").write_text(json.dumps(s))
    with pytest.raises(ValueError,match="failed Q20"):
        run(p)
