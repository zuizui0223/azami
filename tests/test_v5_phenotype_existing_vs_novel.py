import numpy as np
from analysis.v5.phenotype_existing_vs_novel import codes,hist_by_cell,gains,evaluate

def test_explicit_outside_values():
    r=np.arange(10.)[:,None]
    c,k=codes(r,np.ones(10),np.array([[-1.],[0.],[9.],[10.]]))
    assert k==7 and c[0]==0 and c[-1]==6
    assert 1<=c[1]<=5 and 1<=c[2]<=5

def test_partition():
    a=np.array([.5,.5,0.]);b=np.array([.1,.7,.2])
    np.testing.assert_allclose(gains(a,b,np.array([1,1,0],bool)),[.2,.2])

def test_joint_unobserved_inside_marginal_bounds():
    r=np.array([[0.,0.],[1.,1.],[2.,2.],[3.,3.],[4.,4.]])
    x=np.array([[0.,4.],[4.,4.]])
    row,y=evaluate(r,np.ones(5),None,x,np.array([0,1]),2,np.array([1]),np.array([0]),False)
    assert row["analogue_unobserved"]==0 and row["novel_unobserved"]==1
    assert row["new_gain"]==1 and row["existing_gain"]==0

def test_circular_wrap():
    angles=np.radians([1.,359.,360.])
    r=np.column_stack([np.sin(angles),np.cos(angles)])
    c,k=codes(r,np.ones(3),r,True)
    assert k==12 and np.all((c>=0)&(c<12))
