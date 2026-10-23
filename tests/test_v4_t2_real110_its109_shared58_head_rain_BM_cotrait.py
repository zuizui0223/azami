"""Check actual S1 photo-mean + logBIO12 phylogenetic contrast mathematics."""
import numpy as np
import pytest
from analysis.v4 import t2_real110_its109_shared58_head_rain_BM_cotrait as m


def table(n=32):
    rng=np.random.default_rng(304)
    names=[f"Cirsium_source{i:03d}" for i in range(n)]
    data={}
    for i,name in enumerate(names):
        data[name]={
            "chelsa_bio12":str(400+27*i+37*np.sin(i*.33)),
            "mean_image_angle_deg":str(48+3*np.sin(i*.71)+rng.normal(0,2)),
            "original_photo_count":str(i+1),
            "botanical_U_D_state":"U" if i%3 else "D",
        }
    return names,data


def test_analytical_whitened_trait_covariation_has_bounded_uncertainty_and_original_photo():
    names,d=table()
    coords=np.arange(len(names))
    C=np.exp(-abs(coords[:,None]-coords[None,:])/3.)
    x=m.phylogenetic_pair_correlation(names,C,d,.05)
    assert -1<x["S1_original_own_photo_vs_logBIO12_BM_whitened_correlation"]<1
    assert 0<=x["phylogeny_conditioned_approx_two_sided_t_p"]<=1
    assert x["conditional_Fisher_z_95CI"][0]<x["conditional_Fisher_z_95CI"][1]
    assert x["S1_exact_species_count"]==32
    assert x["tip_independent_noise_fraction_of_median_real_gene_tree_depth"]==.05
    assert x["head_angle_is_original_S1_camera_relative_not_gravity"]
    assert x["present_CHELSA_annual_rain_is_not_historical_gain_rain"]


def test_original_botanical_flora_categories_can_all_change_without_changing_photo_result():
    names,source=table()
    n=len(names)
    C=.75*np.eye(n)+.25*np.ones((n,n))
    reference=m.phylogenetic_pair_correlation(names,C,source,.05)
    for row in source.values():
        row["botanical_U_D_state"]="D"
    assert m.phylogenetic_pair_correlation(names,C,source,.05)==reference


def test_root_and_technical_unmeasured_noise_are_not_called_true_evolutionary_rates():
    names,source=table()
    C=np.eye(len(names))
    z=m.phylogenetic_pair_correlation(names,C,source,0.)
    w=m.phylogenetic_pair_correlation(names,C,source,.05)
    assert z["S1_original_own_photo_vs_logBIO12_BM_whitened_correlation"]==pytest.approx(
        w["S1_original_own_photo_vs_logBIO12_BM_whitened_correlation"],abs=1e-12)
    assert z["not_BM_joint_evolutionary_rate_or_independent_gain_count"]
    assert z["spatial_and_biogeographic_confounders_not_adjusted"]
    with pytest.raises(ValueError,match="outside contract"):
        m.phylogenetic_pair_correlation(names,C,source,-.1)
    with pytest.raises(ValueError,match="outside contract"):
        m.phylogenetic_pair_correlation(names,C,source,.4)


def test_actual_t2_two_source_genealogies_and_true_shared_58_scope_are_preserved():
    assert len(m.SOURCE_SCOPES)==4
    assert set(z[2] for z in m.SOURCE_SCOPES)=={"full","shared"}
    assert m.N_EXACT_OVERLAP==58
    assert m.N_MULTILOCUS==110
    assert m.N_ITS==109
    assert m.NUISANCE_TIP_FRACTIONS==(0.,.05)
    assert len(m.EXPECTED_S1_SHA)==64
    assert len(m.EXPECTED_ML108_SHA)==64
    assert len(m.EXPECTED_ITS_SHA)==64
