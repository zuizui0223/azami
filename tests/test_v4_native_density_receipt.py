"""Verify the exact scientific claim boundaries of the 30-fold native output.

These tests preserve already-executed aggregate estimates. They do not
reconstruct private user/source/photo observation identifiers in the repo.
"""
from pathlib import Path
import json
import math

ROOT=Path(__file__).resolve().parents[1]
EVIDENCE=ROOT/"analysis/v4/results/20261021_unified_abcde"


def test_density_robustness_receipt_matches_frozen_source_and_a_to_e_cohort():
    robust=json.loads((EVIDENCE/"native_density_robustness_receipt.json").read_text())
    original=json.loads((EVIDENCE/"verified_aggregate_receipt.json").read_text())
    assert robust["schema_version"]=="native_abcde_density_robustness_20261021_v1"
    source=robust["source"]
    assert (source["workflow_run"],source["artifact_id"])==(39528835793,12314092145)
    assert source["json_member_sha256"] == (
        "0ddf0339e994157ccb75fbc1309917c94265cadaf1dd4fbe80d9baf944483642"
    )
    assert source["source_native_sha256"]==original["provenance"]["source_decoded_native_sha256"]
    assert source["previous_unified_result_sha256"] == (
        "8dd590f286961285372e0879ca9b9c9033ea38b10d16b510ef890c696c1789a3"
    )
    d=robust["design"]
    assert d["eligible_same_nine_outputs"]==157779
    assert d["accepted_taxa_with_measured_core9"]==242
    assert d["qualifying_taxa"]==40
    assert d["repeated_rarefaction"]==d["same_species_null_reps"]==30
    assert d["taxon_cluster_bootstrap"]==1000
    assert d["KDE_bandwidth_multipliers"]==[0.75,1,1.25]
    assert not d["photo_calibration_verified"]


def test_colour_and_outline_share_identical_taxon_pairs_and_taxon_dependent_uncertainty():
    robust=json.loads((EVIDENCE/"native_density_robustness_receipt.json").read_text())
    x=robust["results"]
    n=x["same_species_baseline_taxa"]
    pairs=x["pairs_with_both_same_species_baselines"]
    assert n==38 and pairs==n*(n-1)//2==703
    assert x["colour"]["centred_w2_above_null_pairs"]==498
    assert x["outline"]["centred_w2_above_null_pairs"]==2
    assert x["presentation"]["centred_w2_above_null_pairs"]==37
    for group in ("presentation","colour","outline"):
        assert x[group]["centred_w2_above_null_fraction"]== (
            x[group]["centred_w2_above_null_pairs"]/pairs
        )
    p=x["paired_colour_outline"]
    assert (p["both_above"],p["colour_only_above"],p["outline_only_above"],
            p["neither_above"])==(2,496,0,205)
    assert sum((p["both_above"],p["colour_only_above"],p["outline_only_above"],
                p["neither_above"]))==pairs
    assert p["difference_centred_exceedance_rate"]==(
        p["colour_only_above"]-p["outline_only_above"])/pairs
    lo,med,hi=p["bootstrap_q025_q50_q975"]
    assert 0<lo<med<hi<1
    assert p["leave_one_taxon_out_min_max"][0]>0
    assert p["exclude_top_two_taxa"]["pairs"]==630
    assert p["exclude_top_two_taxa"]["difference"]>0


def test_bandwidth_and_observer_quadrant_sensitivities_remain_only_conditional():
    robust=json.loads((EVIDENCE/"native_density_robustness_receipt.json").read_text())
    bw=robust["results"]["KDE"]
    for module in ("colour","outline"):
        x=bw[module]["overlap_075_1_125"]
        assert len(x)==3 and all(0<=v<=1 for v in x)
        assert x[0]<x[1]<x[2]
        assert all(.9<r<=1 for r in bw[module]["rank_rho_vs_1"])
    assert all(a<b for a,b in zip(
        bw["colour"]["overlap_075_1_125"],
        bw["outline"]["overlap_075_1_125"]
    ))
    deletions=robust["results"]["fixed_deletions"]
    assert deletions["observer_top1pct"]["observers_removed"]==521
    assert deletions["observer_top1pct"]["difference_colour_minus_outline"]>0
    assert deletions["reference_replicates"]==10
    assert deletions["between_species_replicates"]==6
    assert len(deletions["longitude_regions"])==4
    assert all(x["difference"]>0 for x in deletions["longitude_regions"])
    assert all(x["taxa"] <38 for x in deletions["longitude_regions"])
    assert "not a multiple-comparisons-corrected" in robust["interpretation"]["conditional_uncertainty"]
    assert robust["interpretation"]["not_established"]


def test_frozen_v3_contract_preserved_while_new_A_to_E_is_exploratory():
    text=(ROOT/"README.md").read_text()
    assert "analysis/v4/UNIFIED_NATIVE_GATE_A_TO_E_JA.md" in text
    assert "analysis/v3/README.md" in text
    assert "native_density_robustness" in (
        ROOT/".github/workflows/v4-native-integration.yml"
    ).read_text()
    assert (ROOT/"reproducibility/current_reference/manifest.json").is_file()
