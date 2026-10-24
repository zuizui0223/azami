"""All 244 S1 key preservation tests using synthetic source table only."""
import pytest
from analysis.v4 import t2_s1_all244_real241_site_depth_coverage as mod


def mock_original():
    original=[];depth=[]
    for i in range(244):
        key=str(206700000+i)
        matched=i<113
        collided=i>=237
        name=f"Cirsium_test{i:03d}"
        row={
            "accepted_key":key,"source_names":name,
            "source_concept_collision":str(collided),
            "normalized_exact_source_name_candidate":name.replace("_"," "),
            "eligible_photo_observations":"695" if i==0 else "697",
            "Moreyra2025_SRA_taxon_name_candidate":str(matched),
            "S1_mean_EXIF_angle_deg":str(67.+i*.1),
            "chelsa_bio12_photo_eligible_mean":str(1200.+i*.2),
        }
        original.append(row)
        if matched:
            val=120 if i%2 else 60
            depth.append({
                "accepted_key":key,"source_taxon":name,
                "n_independent_public_runs":"1",
                "maximum_per_run_loci_dp3_ge25pct":"180",
                "maximum_per_run_loci_dp3_ge50pct":str(val),
                "maximum_per_run_loci_dp3_100bp":"140",
            })
    return original,depth


def test_unmeasured_nuclear_gene_coverage_is_null_not_fake_zero():
    original,depth=mock_original()
    result,summary=mod.build(original,depth)
    assert len(result)==244
    assert sum(x["source_original_photo_count"] for x in result)==170066
    assert summary["original_clean_accepted_taxa"]==237
    assert summary["actual_public_read_name_candidate_taxa"]==113
    assert result[112]["nuclear_read_coverage_observed"]
    assert result[113]["best_single_run_241_loci_dp3_ge50pct"] is None
    assert result[113]["nuclear_missing_reason"]=="no_exact_name_matched_Moreyra_public_SRA_run"
    assert result[243]["nuclear_missing_reason"]=="WCVP_source_name_collision"
    assert summary["all244_on_independent_nuclear_tree"] is False
    assert summary["actual_read_supported_taxa_ge100_loci_dp3_ge50pct"]==56
    assert summary["photo_depth_sensitivity"]["10"]["n_real_nuclear_run_sample_taxa"]==113


def test_name_collision_cannot_silently_enter_113_nuclear_cohort():
    original,depth=mock_original()
    original[0]["source_concept_collision"]="True"
    with pytest.raises(ValueError,match="WCVP source-name collision"):
        mod.build(original,depth)


def test_taxonomic_source_drift_blocks_promotion():
    original,depth=mock_original()
    depth[0]["source_taxon"]="Cirsium invented"
    with pytest.raises(ValueError,match="Wrong nuclear source taxon"):
        mod.build(original,depth)


def test_missing_one_of_113_source_taxa_fails_closed():
    original,depth=mock_original()
    with pytest.raises(ValueError,match="Incomplete actual"):
        mod.build(original,depth[:-1])
