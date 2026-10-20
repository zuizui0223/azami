import importlib,importlib.util
import numpy as np

def mod():
    assert importlib.util.find_spec('analysis.v4.native_exploration') is not None,'native-exploration executor absent'
    return importlib.import_module('analysis.v4.native_exploration')

def test_moment_r2_equals_direct_standardized_weighted_fit():
    m=mod();rng=np.random.default_rng(44);X=rng.normal(size=(40,4));Y=rng.normal(size=(40,2));w=rng.uniform(.1,2,40)
    a=m.moments(X,Y,w);XZ,YZ,ww=m.standardize(X,Y,w)
    b=np.linalg.lstsq(XZ*np.sqrt(ww[:,None]),YZ*np.sqrt(ww[:,None]),rcond=None)[0]
    expected=1-np.sum(ww[:,None]*(YZ-XZ@b)**2)/np.sum(ww[:,None]*YZ**2)
    assert np.isclose(m.r2(a,list(range(4))),expected)

def test_partition_preserves_signed_shared_variation_and_sums():
    m=mod();rng=np.random.default_rng(5);X=rng.normal(size=(70,17));Y=X[:,[0]]+rng.normal(size=(70,1))
    r=m.environment_summary(X,Y)
    p=r['partition']
    assert np.isclose(p['unique_environment']+p['unique_space']+p['shared']+p['unexplained'],1)
    assert np.isclose(sum(r['shapley_given_space'].values()),p['unique_environment'])

def test_weight_scaling_does_not_change_fit():
    m=mod();rng=np.random.default_rng(8);X=rng.normal(size=(25,3));Y=rng.normal(size=(25,2));w=rng.uniform(.1,2,25)
    assert np.isclose(m.r2(m.moments(X,Y,w),[0,2]),m.r2(m.moments(X,Y,w*99),[0,2]))
