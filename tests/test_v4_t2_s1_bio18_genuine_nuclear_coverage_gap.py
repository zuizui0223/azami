"""The complete native S1 high BIO18 tail is not represented by NA nuclear tree."""
import pytest
from analysis.v4 import t2_s1_bio18_genuine_nuclear_coverage_gap as m
def test_frozen_source_and_nuclear_coverage_universe():
    assert m.SOURCE_251_SHA=="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
    assert m.NEW_S1_FOUR_SHA=="823ecb8442eb03598253cb33f3e16bac86538c6c2652f48ca2dc0886d0e19626"
    assert m.FOUR==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.EXPECTED_CUTOFF==pytest.approx(588.2739645)
    assert m.EXPECTED_WET==63
