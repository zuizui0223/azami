"""Validated genuine nuclear Brownian branches for S1 continuous measured angle."""
from analysis.v4 import t2_true_nuclear_S1_continuous_eazami_brownian_fourclimate as m
def test_genuine_original_frozen_data_not_botanical_U_D():
    assert m.FOUR==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.PANELS=={"NorthAmerica_true_nuclear42":42,"Japan_true_176loci_nuclear12":12}
    assert m.N_PERM==9999
