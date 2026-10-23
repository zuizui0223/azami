"""Validate source-level literature additions cannot silently inflate primary T2."""
from copy import deepcopy
import pytest
from analysis.v4 import t2_s1_251_posture_literature_expansion_guard as audit


def real():
    return audit.load(audit.CENSUS),audit.load(audit.CANDIDATES),audit.load(audit.QUEUE)


def test_real_original_251_and_literature_72_counts():
    all_rows,external,queue=real()
    v=audit.checked_summary(all_rows,external,queue)
    base=v["source_authority_cohorts"]["original_independent_authorities_66"]
    new=v["source_authority_cohorts"]["expanded_single_authority_sensitivity_72"]
    assert v["original_total_taxon_source_name_units"]==251
    assert v["original_accepted_WCVP_keys"]==244
    assert v["original_photographic_observations"]==170066
    assert (base["n_U"],base["n_D"],base["n_original_source_images"])==(54,12,60751)
    assert (new["n_U"],new["n_D"],new["n_original_source_images"])==(57,15,121071)
    assert (base["n_original_real111_source_nuclear_samples"],new["n_original_real111_source_nuclear_samples"])==(38,42)
    assert 900 < new["descriptive_difference_D_minus_U_mm"] < 970
    assert v["remaining_original_source_name_units_unscored"]==177
    assert v["remaining_unscored_original_photo_observations"]==41779
    assert v["remaining_unscored_real111_source_nuclear_candidates"]==67
    assert v["priority_one_real111_nuclear_with_10photos"]==41
    assert not new["phylogenetic_independent_tip_effect_fitted"]
    assert v["global_independent_molecular_species_tree_not_constructed"]
    assert v["rain_induced_evolutionary_transition_not_inferred"]


def test_two_variable_ontogenetic_species_are_kept_unknown():
    all_rows,external,queue=real()
    v=audit.checked_summary(all_rows,external,queue)
    assert set(v["source_documented_variable_or_ontogenetic_taxa"])=={
        "Cirsium heterophyllum","Cirsium japonicum"}
    rows={x["original_source_name"]:x for x in all_rows}
    for taxon in v["source_documented_variable_or_ontogenetic_taxa"]:
        assert rows[taxon]["expanded_posture_sensitivity_U_D_or_unknown"]=="?"
        assert rows[taxon]["new_single_source_literature_posture"]=="V"
        assert rows[taxon]["developmental_or_taxon_polymorphism_flag"]=="True"


def test_cannot_promote_new_single_authority_to_blind_consensus():
    all_rows,external,queue=real()
    bad=deepcopy(external)
    bad[0]["strict_dual_blind_taxon_concept_checked"]="yes"
    with pytest.raises(ValueError,match="single-authority"):
        audit.checked_summary(all_rows,bad,queue)


def test_cannot_convert_cam_relative_angle_or_variable_head_into_nodding():
    all_rows,external,queue=real()
    x=deepcopy(all_rows)
    next(r for r in x if r["original_source_name"]=="Cirsium heterophyllum")[
        "expanded_posture_sensitivity_U_D_or_unknown"]="D"
    with pytest.raises(ValueError,match="Stage-varying"):
        audit.checked_summary(x,external,queue)
    x=deepcopy(all_rows)
    x[0]["botanical_state_not_exif_angle_threshold"]="False"
    with pytest.raises(ValueError,match="Original botanical source"):
        audit.checked_summary(x,external,queue)


def test_mutation_of_frozen_66_or_177_pending_is_rejected():
    all_rows,external,queue=real()
    x=deepcopy(all_rows)
    target=next(r for r in x if r["botanical_U_D_state"]=="D")
    target["botanical_U_D_state"]="U"
    with pytest.raises(ValueError,match="baseline"):
        audit.checked_summary(x,external,queue)
    x=deepcopy(queue)
    x.pop()
    with pytest.raises(ValueError,match="177 still-unscored"):
        audit.checked_summary(all_rows,external,x)


def test_no_fabricated_exact_species_accepted_concept_for_new_literature():
    all_rows,external,queue=real()
    bad=deepcopy(external)
    bad[0]["original_source_name"]="Cirsium never_recorded"
    with pytest.raises(ValueError,match="Published sources not matched"):
        audit.checked_summary(all_rows,bad,queue)
