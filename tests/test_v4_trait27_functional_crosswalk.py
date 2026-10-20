"""Keep the v4 function-led translation in exact lockstep with all 27 source endpoints.

These gates do not certify that photograph-derived quantities are calibrated
biological functions or that new anatomical keypoints have been measured.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"analysis/v4/trait27_functional_crosswalk_v1.json"
SOURCE=ROOT/"legacy/ch1_global/v2/ontology/ch1_continuous_trait_contract.csv"
INVENTORY=ROOT/"analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_endpoint_inventory.csv"

def load():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))

def test_every_frozen_endpoint_retained_exactly_once_with_true_frozen_status():
    v=load()
    original=list(csv.DictReader(SOURCE.open(encoding="utf-8-sig",newline="")))
    inventory=list(csv.DictReader(INVENTORY.open(encoding="utf-8-sig",newline="")))
    endpoints=v["endpoints"]
    assert len(original)==len(inventory)==len(endpoints)==27
    names=[x["id"] for x in endpoints]
    assert len(names)==len(set(names))
    assert set(names)=={x["endpoint_id"] for x in original}
    original_by_id={x["endpoint_id"]:x for x in original}
    inv_by_id={x["endpoint_id"]:x for x in inventory}
    for item in endpoints:
        o=original_by_id[item["id"]]
        x=inv_by_id[item["id"]]
        assert item["source_tier"]==o["analysis_tier"]
        assert item["source_formula_id"]==o["formula_id"]
        assert item["source_v2_observations"]==int(x["n_observations_measured"])
        assert item["source_v2_status"]==x["measurement_status"]
        assert item["storage"].startswith("raw_per_head")
    assert sum(x["source_v2_observations"]==0 for x in endpoints)==5
    assert sum(x["source_v2_observations"]>0 for x in endpoints)==22

def test_modules_are_biological_roles_not_27_independent_regressions():
    v=load()
    assert v["primary_modules"]==["presentation","visible_floret_colour","involucre_structure"]
    assert set(v["contextual_modules"])=={"whole_head_geometry","display_context"}
    groups=Counter(e["mapping_group"] for e in v["endpoints"])
    expected={"presentation":1,"flower_colour":4,"colour_composition":4,"floral_display":1,
              "whole_head_shape":4,"involucre_body":3,"phyllary_contour_proxies":6,
              "involucre_surface_texture":3,"involucre_specularity":1}
    assert groups==expected
    assert sum(groups.values())==27
    for row in v["dependent_measurement_groups"]:
        assert groups[row["group"]]==row["n_endpoints"]
        assert row["needed_new_measurement"]

def test_photo_glare_is_not_sticky_exudate_and_outline_is_not_phyllary_spine():
    v=load();m={x["id"]:x for x in v["endpoints"]}
    spec=m["involucre_surface_specular_fraction"]
    assert spec["v4_role"]=="technical_validation_only_not_stickiness"
    assert spec["biological_module"]=="involucre_structure"
    assert m["bract_projection_maximum"]["v4_role"]=="candidate_contour_proxies_not_true_phyllary_anatomy"
    assert m["capitulum_outline_aspect_ratio"]["biological_module"]=="whole_head_geometry"
    required={x["id"] for x in v["bridge_measurements_not_in_legacy27"]}
    assert {"sticky_exudate_or_gland_assay","phyllary_free_tip_axis_angle",
            "gravity_calibrated_head_angle","floral_pixel_Lab_joint_distribution"}<=required
    assert all(x["unmeasured_is_null"] for x in v["bridge_measurements_not_in_legacy27"])

def test_hue_and_compositions_are_not_separate_univariate_discoveries():
    v=load();rules=v["rules"]
    assert rules["no_hue_two_independent_traits"]
    assert rules["no_colour_pixel_fractions_four_independent_predictors"]
    assert rules["no_27_independent_pvalues"]
    assert rules["no_27_endpoint_complete_case"]
    assert rules["no_image_vertical_as_gravity"]
    assert rules["no_phenotype_correlation_as_genetic_module"]
    assert rules["phylogenetic_placement_only_as_sensitivity"]
