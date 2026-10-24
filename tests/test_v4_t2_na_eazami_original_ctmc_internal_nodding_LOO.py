"""EAzami original internal-only score and all independent author nodding species LOO."""
import math
from analysis.v4 import t2_na_eazami_original_ctmc_internal_nodding_LOO as m
def test_frozen_exact_baseline_and_all_three_D_LOO_contract():
    assert m.BASELINE_SHA=="c97756b6806128aefdca7d09701474c7d21c0bd4ffd6afb9ce27f68369ecba03"
    assert math.comb(42,3)==11480
    assert math.comb(41,2)==820
    assert m.FOUR==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
