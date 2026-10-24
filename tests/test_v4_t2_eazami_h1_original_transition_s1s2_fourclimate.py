"""Ensure frozen EAzami CTMC U->D workflow is used, not previous PIC/PGLS."""
import math
import numpy as np
from analysis.v4 import t2_eazami_h1_original_transition_s1s2_fourclimate as m

def test_frozen_source_and_four_s1s2_climates():
    assert m.AXES==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.PANELS["n5_primary"]==("5",12,7,5,792)
    assert m.PANELS["n3_sensitivity"]==("3",13,7,6,1716)
    assert m.PANELS["n10_strict_sensitivity"]==("10",9,5,4,126)
    assert math.comb(12,5)==792 and math.comb(13,6)==1716
    assert m.EA_H1_SOURCE_BLOB=="09731b563907b3ec2778ea04086d729771376f7e"
    assert m.EA_SIX_TREE_BLOB=="d7f8c0b96ee9c257c815db06a5a0aa9fdae4841e"
    assert m.EA_UD_BLOB=="f655544e77f1659323c1e915567dacbb51b01cc3"

def test_exact_rank_one_sided_positive_bio18_and_both_tails():
    data=np.array([-2.,-1.,0.,1.,2.])
    plus=m.tails_null(data,1.,"positive")
    assert plus["higher_tail_count"]==2
    assert plus["one_sided_directional_fraction"]==.4
    assert plus["two_sided_finite_map_fraction"]==.8
    minus=m.tails_null(data,-1.,"two_sided")
    assert minus["higher_tail_count"]==4
    assert minus["one_sided_directional_fraction"] is None

def test_source_content_blob_sha_provenance():
    assert len(m.EA_HEAD_COMMIT)==40
    assert len(m.EA_CONTRACT_BLOB)==40
    assert len(m.EA_COVERAGE_BLOB)==40
    assert len(m.EA_REFERENCE_BLOB)==40
    assert m.MIN_COVERAGE==.9
