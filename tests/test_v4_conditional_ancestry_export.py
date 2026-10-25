from io import StringIO
import numpy as np
import pytest
from Bio import Phylo
from analysis.v4 import s1s2_multitrait_branch_core as core
from analysis.v4 import conditional_ancestry_export as m

def star():
    return Phylo.read(StringIO('('+','.join(f'T{i}:1' for i in range(8))+');'),'newick')

def test_star_analytic_mean_rate_uncertainty_and_tip_conservation():
    x=np.arange(8,dtype=float)[:,None]
    out=m.reconstruct(star(),[f'T{i}' for i in range(8)],x,['angle'])
    root=out['nodes'][0]
    assert root['estimate']['angle']==pytest.approx(3.5)
    assert out['diffusion_rate']['angle']==pytest.approx(6.)
    assert root['conditional_se']['angle']==pytest.approx(np.sqrt(6/8))
    assert [n['estimate']['angle'] for n in out['nodes'] if n['is_tip']]==list(range(8))
    assert all(n['conditional_se']['angle']==0 for n in out['nodes'] if n['is_tip'])
    assert out['measurement_error_included'] is False

def test_exported_edges_match_frozen_branch_operator():
    tree=Phylo.read(StringIO('(((T0:1,T1:2):1,(T2:1,T3:1):2):1,((T4:1,T5:1):1,(T6:2,T7:1):1):1);'),'newick')
    names=[f'T{i}' for i in range(8)]; x=np.column_stack((np.arange(8),np.arange(8)**2))
    out=m.reconstruct(tree,names,x,['a','b'])
    D,internal=core.branch_operator(tree,names)
    assert np.allclose([[e['delta'][k] for k in ['a','b']] for e in out['edges']],D@x)
    assert [e['internal'] for e in out['edges']]==internal.tolist()
    assert len({n['node_id'] for n in out['nodes']})==len(out['nodes'])

def test_invalid_coordinates_and_duplicate_tips_fail_closed():
    names=[f'T{i}' for i in range(8)]
    with pytest.raises(ValueError):m.reconstruct(star(),names,np.full((8,1),np.nan),['x'])
    tree=star();tree.get_terminals()[1].name='T0'
    with pytest.raises(ValueError):m.reconstruct(tree,names,np.ones((8,1)),['x'])

def test_zero_length_regularization_is_disclosed():
    tree=star();tree.get_terminals()[0].branch_length=0
    out=m.reconstruct(tree,[f'T{i}' for i in range(8)],np.arange(8)[:,None],['x'])
    assert out['zero_length_edges_regularized']==1
    assert out['branch_length_floor']==1e-10

def test_holm_resolution_detects_unreachable_alpha():
    low=core.permutation_resolution(999,84)
    assert low['minimum_attainable_first_holm']==pytest.approx(.084)
    assert low['can_reject_at_alpha'] is False
    assert core.permutation_resolution(9999,84)['can_reject_at_alpha'] is True
    assert core.permutation_resolution(999,40)['can_reject_at_alpha'] is True
    with pytest.raises(ValueError):core.permutation_resolution(0,84)

def test_runner_exports_same_cohort_raw_ancestors_on_all_edges_only():
    import pandas as pd
    from analysis.v4.run_s1s2_multitrait_nuclear_branch import _test_unit
    names=[f'T{i}' for i in range(8)]
    rng=np.random.default_rng(55)
    frame=pd.DataFrame(rng.normal(size=(8,5)),index=names,columns=['angle',*core.FOUR])
    frame['n_original_observations']=10
    info={'tree':star(),'names':names}
    out=_test_unit('test',info,frame,('angle',),'all_edges',999,1,'endpoint:angle')
    assert out['conditional_ancestry']['nodes'][0]['estimate']['angle']==pytest.approx(frame.angle.mean())
    assert set(out['conditional_ancestry']['diffusion_rate'])=={'angle',*core.FOUR}

def test_negative_length_is_rejected():
    tree=star(); tree.get_terminals()[0].branch_length=-1
    with pytest.raises(ValueError,match='branch length'):
        m.reconstruct(tree,[f'T{i}' for i in range(8)],np.arange(8)[:,None],['x'])
