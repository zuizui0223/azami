"""Original-photo public-aggregate curve recovery and global region reassembly."""
import copy
import numpy as np
import pytest

from analysis.v4.s2_region_curve_diagnostics import reconstruct,INPUT_STATUS,SPECIES


def _fixture():
    nr=12
    na=12
    n=np.ones((nr,na),float)
    n=n/n.sum()
    # A deliberately different introduced angle distribution in two regions.
    left=np.ones((nr,na),float)
    left[:,0]*=3.
    left/=left.sum()
    right=np.ones((nr,na),float)
    right[:,-1]*=3.
    right/=right.sum()
    intro=.7*left+.3*right
    def panel(mat):
        probs=mat.sum(axis=1)
        return {
            "original_photo_observations":2500,
            "original_geographic_cells":55,
            "original_wild_source_rain_bin_mass":probs.tolist(),
            "original_photo_rain_bin_mass":probs.tolist(),
            "original_photo_joint_mass":mat.tolist(),
            "photo_angle_conditional_on_rain_or_null":
                (mat/probs[:,None]).tolist(),
        }
    def group(mat,weight):
        return {
            "photo_n":220 if weight>.5 else 35,
            "photo_cells":40 if weight>.5 else 6,
            "fraction_global_cell_weight":weight,
            "joint_mass_within_region":mat.tolist(),
        }
    case={
        "display_climate_BIO12_edges_mm_year":list(range(100,230,10)),
        "display_angle_edges_image_degrees":list(range(0,195,15)),
        "native":panel(n),
        "all_introduced":panel(intro),
        "introduced_continents_supplementary":{
            "groups":{
                "North America":group(left,.7),
                "other_low_support":group(right,.3),
            }
        }
    }
    return {
        "status":INPUT_STATUS,
        "all_introduced_remains_primary":True,
        "visualization_only_no_trait_width_or_tail_indices_selected":True,
        "privacy":{"obs_ids":False,"photographers":False},
        "species":{name:copy.deepcopy(case) for name in SPECIES},
    }


def test_exact_full_curve_recovery_preserves_low_support_in_main():
    data=_fixture()
    output=reconstruct(data)
    assert output["status"]=="DERIVED_FROM_PREVERIFIED_SOURCE_AGGREGATE_ONLY"
    assert output["no_selected_quantile_width_or_extreme_angle_indicator"]
    for case in output["species"].values():
        assert case["regional_reconstruction_exact"]
        assert len(case["all_introduced"]["same_photo_image_angle_ECDF"])==12
        assert case["all_introduced"]["photo_eligible_BIO12_ECDF"][-1]==pytest.approx(1)
        assert case["all_introduced"]["original_wild_source_BIO12_ECDF"][-1]==pytest.approx(1)
        assert len(case["same_rain_photo_angle_conditional_difference_or_null"])==12
        # Exclusion is always a separate sensitivity: the full primary curve
        # still includes low-support regional observations.
        assert case["introduced_supported_regions_only_ECDF_sensitivity"]!=(
            case["all_introduced"]["same_photo_image_angle_ECDF"])


def test_mixture_mismatch_is_not_silently_renormalized():
    data=_fixture()
    bad=data["species"][SPECIES[0]]["introduced_continents_supplementary"][
        "groups"]["North America"]
    bad["fraction_global_cell_weight"]=.8
    with pytest.raises(ValueError,match="could not be reassembled"):
        reconstruct(data)


def test_refuse_source_or_privacy_degradation():
    data=_fixture()
    data["privacy"]["photos"]=True
    with pytest.raises(ValueError,match="authentic|privacy-safe"):
        reconstruct(data)
    data=_fixture()
    data["species"][SPECIES[0]]["native"]["original_photo_joint_mass"][0][0]=-.001
    with pytest.raises(ValueError,match="two-dimensional"):
        reconstruct(data)
