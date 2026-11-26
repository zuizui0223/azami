import numpy as np
from analysis.v5.phenotype_native_transport import neighbors,predict_native
def test_native_prediction_and_support():
    x=np.column_stack([np.arange(40.),np.arange(40.)**2])
    y=np.column_stack([np.ones(40),np.arange(40.)])
    ca=x[[3,12,26,35]]
    query=np.vstack([x[12],x[12]+100000])
    pc,pi,ok,q=predict_native(x,y,ca,query,k=3)
    np.testing.assert_allclose(pc[:,0],1)
    np.testing.assert_allclose(pi[:,0],1)
    assert bool(ok[0]) and not bool(ok[1])
    np.testing.assert_allclose(pc[1],pi[0])
def test_k1_shape():
    x=np.array([[0.,0.],[1.,2.],[2.,4.]])
    d,i=neighbors(x,x,k=1)
    np.testing.assert_allclose(d,0)
    assert i.shape==(3,1)
