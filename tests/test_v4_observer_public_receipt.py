"""Confirm only aggregate, internally consistent observer-bias evidence is versioned."""
import json
from pathlib import Path
import math
from collections.abc import Mapping

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"analysis/v4/observer_effort_aggregate_20261020.json"

def test_public_receipt_contains_no_identifiers_or_locations():
    r=json.loads(RECEIPT.read_text())
    assert r["source"]["native_observations"]==313350
    assert r["source"]["native_accepted_taxa"]==306
    assert r["source"]["native_observers"]==92839
    forbidden={"obs_id","photo_id","user_id","latitude","longitude",
               "photo_url","observation_id","url","original_path"}
    def walk(d):
        if isinstance(d,Mapping):
            assert not (set(d)&forbidden),set(d)&forbidden
            for v in d.values():walk(v)
        elif isinstance(d,list):
            for v in d:walk(v)
    walk(r)

def test_pooled_between_and_within_contributions_close_exactly():
    r=json.loads(RECEIPT.read_text())
    p=r["pooled_color_radiation"]
    for x in p.values():
        assert math.isclose(x["pooled_beta"],
                            x["between_contribution"]+x["within_contribution"],
                            abs_tol=2e-5)
    assert p["observation"]["pooled_beta"]<0
    assert p["equal_taxon_region_observer"]["pooled_beta"]<p["observation"]["pooled_beta"]
    assert p["observation_omit_top_two_taxa"]["pooled_beta"]<p["observation"]["pooled_beta"]
    assert p["observation_omit_top_five_taxa"]["pooled_beta"]<p["observation_omit_top_two_taxa"]["pooled_beta"]

def test_aggregate_claim_ceiling_keeps_spatial_uncertainty():
    r=json.loads(RECEIPT.read_text())
    s=r["joint_environment_plus_space"]
    assert s["color_radiation_beta"]<0
    assert s["presentation_rain_beta"]<0
    assert s["independent_environmental_causation_supported"] is False
    for interval in (s["color_radiation_beta"], s["presentation_rain_beta"]):
        assert math.isfinite(interval)
    p=r["orientation_color_integration_observer_balanced"]
    assert p["bootstrap_matched_n_interval_95"][0]<0<p["bootstrap_matched_n_interval_95"][1]
