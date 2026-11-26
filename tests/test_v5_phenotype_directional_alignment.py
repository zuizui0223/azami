import numpy as np
from analysis.v5.phenotype_directional_alignment import existing_tails, alignment, exact_block_test

def test_native_tails_exclude_outside_observed_range():
    ref=np.arange(10.,dtype=float)
    x=np.array([-1.,0.,5.,9.,10.])
    h,th=existing_tails(x,np.zeros(len(x),int),1,ref,np.ones(len(ref)))
    np.testing.assert_allclose(h[0],[.2,.2,.4])
    assert th==(0.,1.,7.,9.)

def test_aligned_tail_swaps_when_climate_reverses():
    b=np.repeat(np.arange(8),2)
    dx=np.tile([-1.,1.],8)
    dy=np.column_stack([-dx*.1,dx*.1])
    a,c=alignment(dx,dy,b)
    np.testing.assert_allclose(a["effect"],.2)
    np.testing.assert_allclose(a["predicted_tail_change"],.1)
    np.testing.assert_allclose(a["opposite_tail_change"],-.1)
    np.testing.assert_allclose(c["effect"],.2)
    assert a["p"]==c["p"]==2/256

def test_group_shift_not_a_within_region_gradient():
    b=np.repeat(np.arange(8),2)
    dx=np.tile([1.,2.],8)
    dy=np.tile([-.1,.1],(16,1))
    a,c=alignment(dx,dy,b)
    np.testing.assert_allclose(a["effect"],.2)
    np.testing.assert_allclose(c["effect"],0.)
    assert c["p"]==1

def test_single_pair_blocks_cannot_identify_within_gradient():
    a,c=alignment(np.ones(8),np.tile([-.1,.1],(8,1)),np.arange(8))
    assert a["status"]=="estimated"
    assert c["status"]=="insufficient_support" and c["p"]==1
