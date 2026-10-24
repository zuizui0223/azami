"""Source-verified S2 equivalence feasibility without retroactive SESOI selection."""
import json

import numpy as np
import pytest

from analysis.v4 import s2_ecdf_equivalence_frontier as eq


def fake_curve(sample,unit,shift=0.35,radius=.10):
    x=np.linspace(0,1,301)
    n=np.clip(x,0,1)
    i=np.clip(x-shift,0,1)
    i[-1]=1.0
    delta=i-n
    low=np.maximum(-1,delta-radius)
    high=np.minimum(1,delta+radius)
    return {
       "unit":unit,
       "native_original_observations":sample[0],
       "all_introduced_original_observations":sample[1],
       "native_and_introduced_geographic_cells_equal_weight":True,
       "display_only_unselected_common_knots":x.tolist(),
       "native_empirical_cdf":n.tolist(),
       "all_introduced_empirical_cdf":i.tolist(),
       "whole_curve_spatial_uncertainty":{
         "status":"EXPLORATORY_FULL_ECDF_GEOGRAPHIC_CELL_UNCERTAINTY",
         "resamples_requested":399,
         "resamples_valid":399,
         "simultaneous_confidence_envelope_exploratory_conditional_on_source":True,
         "introduced_minus_native_ECDF":delta.tolist(),
         "simultaneous_95pct_spatial_sensitivity_lower":low.tolist(),
         "simultaneous_95pct_spatial_sensitivity_upper":high.tolist(),
         "native_5degree_cells":45,
         "introduced_5degree_cells":54,
       }
    }


def source_example():
    return {
      "status":"ORIGINAL_NATIVE_ALL_INTRODUCED_FULL_EMPIRICAL_CURVES",
      "no_selected_outcome_threshold_or_effect_score":True,
      "complete_ECDF_geographic_uncertainty_without_selected_metric":True,
      "privacy":{"observation_id_exported":False,"photographer_id_exported":False,
                 "source_location_exported":False},
      "species":{
        k:{
          "rainfall_environment_full_distribution":fake_curve(v["rain"],"mm_per_year",.17),
          "image_angle_full_distribution":fake_curve(v["photo"],"image_degrees",.05),
        } for k,v in eq.SPECIES.items()
      }
    }


def test_equivalence_margin_upper_is_NOT_a_formal_equivalence_test():
    r=eq.inspect_curve(fake_curve((21870,25657),"image_degrees",shift=.05,radius=.11),
         n_native=21870,n_introduced=25657,unit="image_degrees")
    assert 0<r["observed_sup_absolute_CDF_difference"]<.1
    assert r["simultaneous_95pct_envelope_upper_bound_on_sup_CDF_difference"]>.1
    assert r["minimum_delta_for_95pct_simultaneous_band_containment"]>=r[
      "observed_sup_absolute_CDF_difference"]
    assert r["NOT_a_prospectively_chosen_SESOI"] is True
    assert r["NOT_a_calibrated_formal_equivalence_TOST"] is True
    assert r["NOT_a_proof_of_unchanged_biological_head_posture"] is True


def test_report_cannot_claim_equivalence_when_band_zero_crosses():
    source=source_example()
    for sp in source["species"].values():
        band=sp["image_angle_full_distribution"]
        v=eq.inspect_curve(band,
              n_native=band["native_original_observations"],
              n_introduced=band["all_introduced_original_observations"],
              unit="image_degrees")
        assert v["all_95pct_CDF_difference_bands_include_zero_at_reported_knots"]


def test_original_SHA_prevents_artifact_replacement(tmp_path,monkeypatch):
    data=source_example()
    raw=(json.dumps(data,sort_keys=True)+"\n").encode()
    p=tmp_path/"S2_unselected_full_rain_and_head_angle_distributions.json"
    p.write_bytes(raw)
    monkeypatch.setattr(eq,"ORIGINAL_SOURCE_SHA256",__import__("hashlib").sha256(raw).hexdigest())
    checked=eq.input_check(p)
    assert checked["species"].keys()==data["species"].keys()
    out=tmp_path/"result.json"
    result=eq.evaluate(p,out)
    assert result["status"].startswith("ACTUAL_SOURCE_GRID")
    assert out.exists()
    assert set(result["species"])=={"Cirsium arvense","Cirsium vulgare"}
    with pytest.raises(ValueError,match="unverified"):
        p.write_bytes(raw+b"tampered")
        eq.input_check(p)


def test_bad_original_ecdf_denominator_uncertainty_or_source_precision_fails():
    source=fake_curve((21870,25657),"image_degrees")
    with pytest.raises(ValueError,match="denoms"):
        eq.inspect_curve(source,n_native=21871,n_introduced=25657,unit="image_degrees")
    source["whole_curve_spatial_uncertainty"]["resamples_valid"]=100
    with pytest.raises(ValueError,match="resamples"):
        eq.inspect_curve(source,n_native=21870,n_introduced=25657,unit="image_degrees")
    source["whole_curve_spatial_uncertainty"]["resamples_valid"]=399
    source["whole_curve_spatial_uncertainty"]["simultaneous_95pct_spatial_sensitivity_lower"]=[1]*301
    with pytest.raises(ValueError,match="inconsistent"):
        eq.inspect_curve(source,n_native=21870,n_introduced=25657,unit="image_degrees")
