"""Synthetic validation for continent-specific original image-head orientation.

Actual published evidence must come from authenticated 665115 source Actions;
no invented toy/photo effects are included in source receipts.
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4 import native_intro_continent_orientation as site
from analysis.v4.native_cached_climate import VARIABLES


def toy_source():
    rng=np.random.default_rng(20261021)
    data=[]
    cases=[
        ("native",1,40.,0.,0.),
        ("introduced",7,49.,7.,0.),
        ("introduced",2,36.,0.,8.),
        ("introduced",8,55.,0.,7.),
        ("introduced",5,44.,0.,0.),
    ]
    for species in sorted(site.FOCAL):
        for state,continent,angle,warm,seasonal in cases:
            for cell in range(40):
                for k in range(4):
                    row={
                        "obs_id":f"{species}-{state}-{continent}-{cell}-{k}",
                        "taxon_name":species,"range_state":state,
                        "user_id":f"person-{species}-{state}-{continent}-{cell}-{k}",
                        "tdwg1":continent,
                        "latitude":float(-25),
                        "longitude":float(-170+cell*7),
                        "observed_on":f"2024-{(cell%12+1):02d}-18",
                        "quality_grade":"research",
                        "orientation_angle_degrees_median":float(angle+rng.normal(scale=1.5)),
                        "source_photo_count":1,
                    }
                    for var in VARIABLES:
                        row[var]=float(rng.normal(scale=.25))
                    # Half introduced locations are novel warmer / more seasonal;
                    # the other half are a climate analogue within this continent.
                    if state=="introduced" and cell>=20:
                        row["chelsa_bio01"]+=warm
                        row["chelsa_bio15"]+=seasonal
                    data.append(row)
    return pd.DataFrame(data)


def test_continental_photo_angles_compare_original_native_and_multiple_invasions(monkeypatch):
    f=toy_source()
    monkeypatch.setattr(site,"EXPECTED_FOCAL_SOURCE_N",len(f))
    result=site.analyse_regions(f)
    assert set(result["source_taxa"])==site.FOCAL
    assert result["no_private_original_identifiers_in_public_result"]
    assert result["two_species_after_original_observer_cell_date_balancing"]==len(f)
    for species in site.FOCAL:
        row=result["source_taxa"][species]
        assert row["native_all_regions"]["source_photo_orientation_measurements"]==160
        regions=row["introduced_by_TDWG_level1"]
        assert set(regions)=={"North America","Africa","South America","Australasia"}
        assert regions["North America"]["vs_native_descriptive_difference_deg"]>5
        assert regions["South America"]["vs_native_descriptive_difference_deg"]>10
        assert regions["Africa"]["vs_native_descriptive_difference_deg"]<0
        assert all(region["vs_native_cell_mean_angle"]["status"]=="ORIGINAL_SOURCE_CELL_BOOTSTRAP_DESCRIPTIVE"
                   for region in regions.values())
        for region in regions.values():
            assert region["original_5degree_equal_sin_cells"]==40
            assert region["vs_native_cell_mean_angle"]["bootstrap_site_cell_replicates"]==1000
            assert len(region["photo_orientation_quantiles_q10_q50_q90"])==3
            assert region["vs_native_cell_mean_angle"]["bootstrap_intervals_do_not_identify_technical_camera_or_genetic_change"]
    public=json.dumps(result)
    assert "person-Cirsium" not in public
    assert '"obs_id"' not in public
    assert '"latitude"' not in public
    assert result["source_pixel_climate_coverage_by_region"]=={"status":"NOT_COMPUTED"}


def test_continental_original_climate_proxy_identifies_warm_and_seasonality_groups_without_angle_use(monkeypatch):
    f=toy_source()
    monkeypatch.setattr(site,"EXPECTED_FOCAL_SOURCE_N",len(f))
    arrays={p:f[p].to_numpy(float) for p in VARIABLES}
    result=site.analyse_regions(f,climate=arrays)
    key="Cirsium vulgare"
    climate=result["climate_novelty_x_photo_angle_by_introduced_region"]
    assert climate[key+"|North America|chelsa_bio01"]["status"].startswith(
        "SUPPORTED_DESCRIPTIVE_IMAGE_ANGLE")
    for continent in ("Africa","South America"):
        case=climate[key+"|"+continent+"|chelsa_bio15"]
        assert case["status"].startswith("SUPPORTED_DESCRIPTIVE_IMAGE_ANGLE")
        assert case["introduced_beyond_native_q95"]["original_5degree_equal_sin_cells"]>=15
    assert result["source_pixel_climate_coverage_by_region"][key+"|introduced|Africa"][
        "five_original_CHELSA_pixels_complete"]==160
    native=f.copy()
    native["orientation_angle_degrees_median"]=179-np.asarray(native["orientation_angle_degrees_median"])
    changed=site.analyse_regions(native,climate=arrays)
    for k,v in climate.items():
        assert v.get("reference_native_cell_equal_q95")==changed[
            "climate_novelty_x_photo_angle_by_introduced_region"][k].get(
                "reference_native_cell_equal_q95")


def test_no_source_pixels_in_introduced_Africa_means_no_wetness_response(monkeypatch):
    f=toy_source()
    monkeypatch.setattr(site,"EXPECTED_FOCAL_SOURCE_N",len(f))
    arrays={p:f[p].to_numpy(float).copy() for p in VARIABLES}
    mask=f.tdwg1.eq(2).to_numpy()
    for value in arrays.values():
        value[mask]=np.nan
    a=site.analyse_regions(f,climate=arrays)
    for taxon in site.FOCAL:
        row=a["climate_novelty_x_photo_angle_by_introduced_region"][
            taxon+"|Africa|chelsa_bio15"]
        assert row["status"]=="STOP_NO_ORIGINAL_CLIMATE_COVERAGE_OR_INDEPENDENT_CELLS"
        assert row["within_continent_angle_change"] is None
        assert a["source_taxa"][taxon]["introduced_by_TDWG_level1"]["Africa"][
            "source_photo_orientation_measurements"]==160


def test_no_cell_based_certainty_for_two_patch_introduction():
    toy=toy_source()
    one=toy.loc[toy.tdwg1.eq(7)].head(20).copy()
    one["_cell"]=np.arange(len(one))%2
    assert site.cell_weighted_metrics(one)["status"]=="INSUFFICIENT_SOURCE_OBSERVERS_OR_SPATIAL_CELLS"
    lo=site._paired_cell_bootstrap(np.arange(20,dtype=float),np.arange(2,dtype=float),"test")
    assert lo["status"]=="STOP_INSUFFICIENT_SOURCE_CELLS_FOR_INTERVAL"
    assert lo["interval"] is None


def test_invalid_camera_angles_are_not_assigned_as_true_downturned_heads():
    toy=toy_source()
    one=toy.head(30).copy()
    one["_cell"]=np.arange(30)
    one.iloc[0,one.columns.get_loc("orientation_angle_degrees_median")]=190.
    with pytest.raises(ValueError,match="Invalid original image-vertical angle"):
        site.cell_weighted_metrics(one)
