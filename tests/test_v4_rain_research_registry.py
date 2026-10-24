"""S1 spatial evidence and S2/T2/T3 inference guardrails for current V4.

Tests ensure new rain-shelter research design cannot turn unexecuted invasion or
evolution into a claimed numerical result or treat all non-native as introduced.
"""
from pathlib import Path
import json

import pandas as pd
import pytest

from analysis.v4.native_range_status import distribution_status, inclusion_reason

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"
CONTRACT=ROOT/"analysis/v4/RAIN_ORIENTATION_SPATIAL_INTRO_EVOLUTION_20261021_JA.md"


def _evidence():
    return json.loads(SOURCE.read_text(encoding="utf-8"))


def test_induced_WCVP_introduced_status_not_generic_foreign_or_uncertain():
    assert distribution_status({
        "locationId":"TDWG:ABC","occurrenceStatus":"PRESENT",
        "establishmentMeans":"INTRODUCED"})=="introduced"
    assert distribution_status({
        "locationId":"TDWG:ABC","occurrenceStatus":"PRESENT",
        "establishmentMeans":"NATIVE"})=="native"
    assert distribution_status({
        "locationId":"TDWG:ABC","occurrenceStatus":"PRESENT",
        "establishmentMeans":""})=="native"
    assert distribution_status({
        "locationId":"TDWG:ABC","occurrenceStatus":"ABSENT",
        "establishmentMeans":"INTRODUCED"})=="uncertain"
    assert distribution_status({"locationId":"TDWG:ABC","occurrenceStatus":"PRESENT",
                                "establishmentMeans":"UNKNOWN"})=="uncertain"
    assert distribution_status({"locationId":"country-name","occurrenceStatus":"PRESENT",
                                "establishmentMeans":"INTRODUCED"})=="uncertain"
    assert inclusion_reason("introduced",True,False,False)=="range_introduced"
    assert inclusion_reason("native",True,False,False)=="include"
    assert inclusion_reason("native",True,False,True)=="captive"


def test_S1_existing_BIO12_effect_matches_saved_preexisting_actual_model():
    b=_evidence()
    y=b["S1"]
    df=pd.read_csv(ROOT/"analysis/v4/results/20261020_native_integration/native_trait_effects.csv")
    row=df.loc[(df.module=="orientation") &
                (df.endpoint=="orientation_image_vertical_angle") &
                (df.term=="chelsa_bio12__between")]
    assert len(row)==1
    row=row.iloc[0]
    assert row.n==170066 and row.species==244
    assert float(row.trait_change)==pytest.approx(y["annual_bio12_between_deg_per_100mm"])
    assert float(row.lower95)==pytest.approx(y["annual_bio12_between_95pct_interval"][0])
    assert float(row.upper95)==pytest.approx(y["annual_bio12_between_95pct_interval"][1])
    assert y["bio12_within_interval_crosses_zero"]
    assert y["phylogeny_scenarios_are_independent_evolutionary_events"] is False
    assert b["source_frame"]["native_primary_observations"]==313350


def test_no_imagined_introduction_or_ancestral_origin_dates():
    b=_evidence()
    assert b["study_structure"]==[
        "S1 global native spatial pattern",
        "S2 native-to-introduced within-species transfer",
        "T2 Cirsium transitions",
        "T3 Asteraceae cross-tribe convergence"
    ]
    # Original registry first described merely checking establishment
    # status; S2 has since executed and authenticated all-introduced
    # comparisons. Refuse a return to an unexecuted-only status while still
    # rejecting fabricated introduction dates or adaptation claims.
    assert "SOURCE-VERIFIED MAIN" in b["S2"]["status"]
    assert "WCVP native vs introduced" in b["S2"]["status"]
    assert "adaptive niche shift not established" in b["S2"]["status"]
    assert b["S2"]["introduction_age_inferred_from_photograph_date"] is False
    assert b["T2"]["existing_52_trees_sufficient_for_dated_origin_counts"] is False
    assert b["T3"]["cross_tribe_convergence_inferred"] is False
    assert b["quality_controls"]["same_observation_is_same_head_verified"] is False
    assert b["quality_controls"]["photographer_random_effect_fitted"] is False
    assert not b["quality_controls"]["monthly_climate_precipitation_12_rasters_joined"]
    assert b["frozen_v3_reference"]=={
        "stages":8, "outputs":16, "unchanged":True,
        "submission_replacement_decision":"conditional_after_new_v4_validation"
    }


def test_readme_uses_spatial_introduction_evolution_as_primary_question():
    readme=(ROOT/"README.md").read_text(encoding="utf-8")
    design=CONTRACT.read_text(encoding="utf-8")
    assert "RAIN_ORIENTATION_SPATIAL_INTRO_EVOLUTION_20261021_JA.md" in readme
    assert "RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json" in readme
    assert "雨はアザミの頭花の向きを形づくってきたか" in design
    for stage in ("S1", "S2", "T2", "T3"):
        assert stage in design
    assert "逆" not in design or "符号" in design
    assert b"".decode() is not None
