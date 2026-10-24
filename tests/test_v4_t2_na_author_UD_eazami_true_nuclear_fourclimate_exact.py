"""Genuine published NA author U/D & EAzami CTMC but NOT S1 camera angle U/D."""
import math
from analysis.v4 import t2_na_author_UD_eazami_true_nuclear_fourclimate_exact as m
def test_full_source_PINNED_real_science_contract():
    assert m.NA_TRAITS=="traits/traits.csv"
    assert m.NA_TRAIT_BLOB=="fa127f0587b676e6e004ea121424d26069586bcb"
    assert m.FOUR==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.SOURCE_S1_SHA=="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
    assert m.SOURCE_NEW4_SHA=="823ecb8442eb03598253cb33f3e16bac86538c6c2652f48ca2dc0886d0e19626"
    assert math.comb(42,3)==11480 and math.comb(40,3)==9880
