import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "reproducibility" / "v3_environment_design_upstream_20260907.json"


def load_receipt():
    return json.loads(RECEIPT.read_text(encoding="utf-8"))


def test_full_source_and_native_denominators_are_locked():
    r = load_receipt()
    assert r["status"] == "V3_ENVIRONMENT_DESIGN_UPSTREAM_COMPLETE_BEFORE_CHELSA_DIAGNOSTICS"
    source = r["source_observations"]
    assert source["unique_observations"] == 665139
    assert source["metadata_rows"] == 1122855
    assert source["api_rows"] == 640141
    assert source["exact_date_rows"] == 663255
    assert source["public_coordinate_rows"] == 632927

    native = r["native_range_join"]
    counts = native["status_counts"]
    assert sum(counts.values()) == 665139
    assert counts["native"] == 384072
    assert counts["introduced"] == 157716
    assert counts["unmapped_or_unusable_location"] == 64283
    assert counts["unresolved_taxon"] == 49281
    assert counts["unlisted"] == 9787
    assert native["primary_native_wild_public_exact_date_resolved_rows"] == 382258


def test_taxon_balancing_is_upstream_and_phenotype_blind():
    r = load_receipt()
    cohort = r["environment_design_cohort"]
    assert cohort["eligible_rows_before_taxon_support"] == 319244
    assert cohort["eligible_taxa_before_support"] == 355
    assert cohort["taxa_meeting_support"] == 243
    assert cohort["selected_rows"] == 9503
    assert cohort["selected_taxa"] == 243
    assert cohort["taxon_cap"] == 60
    assert cohort["minimum_taxon_support"] == 5
    assert abs(cohort["equal_taxon_weight_fraction"] - 1 / 243) < 1e-15
    assert cohort["largest_selected_taxon_fraction"] < cohort["largest_raw_eligible_taxon_fraction"]

    separation = r["separation_checks"]
    assert separation == {
        "trait_files_read": 0,
        "image_files_read": 0,
        "environment_values_used_to_select_cohort": 0,
        "ecological_models_executed": 0,
        "source_rows_deleted": 0,
    }


def test_cancelled_chelsa_step_is_explicitly_not_evidence():
    r = load_receipt()
    workflow = r["source_workflow"]
    assert workflow["upstream_steps_6_to_9"] == "success"
    assert workflow["chelsa_step_10"] == "cancelled_and_not_used_from_this_run"
    assert "cancelled_after_step_9" in workflow["workflow_terminal_state"]
    assert len(workflow["artifact_sha256"]) == 64
    assert set(r["report_sha256"]) == {
        "environment_source_observations_report.json",
        "native_range_join_report.json",
        "environment_design_cohort_report.json",
    }
