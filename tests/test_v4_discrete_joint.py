import itertools
from io import StringIO
import numpy as np
from Bio import Phylo
from analysis.v4.discrete_joint import posterior

def test_three_state_joint_matches_exhaustive_enumeration():
    tree=Phylo.read(StringIO('((a:.2,b:.3):.1,c:.4);'),'newick')
    states={'a':{'red'},'b':{'white'},'c':{'red','yellow'}}
    universe=('red','white','yellow');q=.7
    result=posterior(tree,states,universe,q)
    nodes=list(tree.find_clades(order='preorder'));idx={n:i for i,n in enumerate(nodes)}
    edges=[(idx[p],idx[c],c.branch_length) for p in nodes for c in p.clades]
    free=[i for i,n in enumerate(nodes) if not n.is_terminal()]
    weights=[];assignments=[]
    for ss in itertools.product(range(3),repeat=len(nodes)):
        weight=1/3
        for i,n in enumerate(nodes):
            if n.is_terminal():weight*=int(universe[ss[i]] in states[n.name])/len(states[n.name])
        for p,c,t in edges:
            e=np.exp(-3*q*t);weight*=1/3+((1 if ss[p]==ss[c] else 0)-1/3)*e
        weights.append(weight);assignments.append(ss)
    total=sum(weights)
    assert np.isclose(np.exp(result['log_likelihood']),total)
    for row in result['edges']:
        expected=np.zeros((3,3))
        for ss,w in zip(assignments,weights):expected[ss[row['parent']],ss[row['child']]]+=w/total
        assert np.allclose(row['joint'],expected)
    for row in result['nodes']:
        expected=np.zeros(3)
        for ss,w in zip(assignments,weights):expected[ss[row['node']]]+=w/total
        assert np.allclose(row['probabilities'],expected)

def test_unresolved_tip_is_not_forced_to_a_state():
    tree=Phylo.read(StringIO('(a:.2,b:.3,c:.4);'),'newick')
    r=posterior(tree,{'a':{'0'},'b':{'1'}},('0','1'),.4)
    c=next(n for n in r['nodes'] if n['name']=='c')
    assert all(0<p<1 for p in c['probabilities'])
