"""The observed full S1 / ITS109 gap is source-taxonomy and model-scope, not
an argument for instantaneous plasticity or assumed adaptive radiation."""
import numpy as np
import pytest

from analysis.v4 import t2_s1_109_vs251_phylogenetic_sampling_explanation as m


def mock(n=26):
    output=[]
    for i in range(n):
        record={
            "original_source_name":f"Cirsium_species{i:03}",
            "mean_image_angle_deg":str(40+0.3*i+9*np.sin(i*1.61)),
            "median_image_angle_deg":str(42+0.3*i+9*np.cos(i*.83)),
            "chelsa_bio01":str(9+2*np.sin(.77*i)),
            "chelsa_bio04":str(500+24*i+57*np.sin(.31*i)),
            "chelsa_bio12":str(1000+50*i+75*np.cos(.13*i)),
            "chelsa_bio15":str(36+2*np.cos(.57*i)+.15*i),
            "chelsa_rsds_mean":str(13+1.2*np.sin(.93*i)+.08*i),
        }
        output.append(record)
    return output


def test_HC3_frozen_S1_five_axis_effect_does_not_dichotomise_angle():
    source=mock()
    fit=m.fit_equal_species_5axis_HC3(source,"mean_image_angle_deg")
    assert fit["n_equal_weighted_original_S1_source_name_units"]==26
    assert fit["continuous_S1_photo_response"]=="mean_image_angle_deg"
    assert len(fit["source_covariates"])==5
    assert len(fit["HC3_95percent_normal_ci"])==2
    assert fit["HC3_se"]>0
    assert np.isfinite(fit["BIO12_effect_original_photo_angle_degrees_per100mm"])
    assert fit["photo_angle_is_EXIF_image_vertical_not_botanical_nodding"]
    assert fit["not_original_244_accepted_key_observation_mixed_model"]


def test_nominal_genus_posture_labels_cannot_change_continuous_frozen_photo_fit():
    x=mock()
    first=m.fit_equal_species_5axis_HC3(x,"mean_image_angle_deg")
    for i,r in enumerate(x):
        r["floristic_head_posture"]="D" if i%2==0 else "U"
    second=m.fit_equal_species_5axis_HC3(x,"mean_image_angle_deg")
    assert first==second


def test_source_109_photo_coverage_and_branch_support_expectations_locked():
    assert m.EXPECTED_IMAGES==170066
    assert m.NUCLEAR_TIPS==109
    assert m.SOURCE_NAMES==251
    assert m.ACCEPTED_KEYS==244
    assert len(m.FIVE)==5
    assert "chelsa_bio12" in m.FIVE
    assert m.EXPECTED_251_SHA256=="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
    assert m.EXPECTED_109_ML_SHA256=="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"


def test_no_false_original_uncertainty_when_rank_collapse():
    x=mock()
    for r in x:
        r["chelsa_bio12"]="1000"
    with pytest.raises(ValueError,match="invariant"):
        m.fit_equal_species_5axis_HC3(x,"mean_image_angle_deg")
