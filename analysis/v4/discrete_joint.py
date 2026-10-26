"""Equal-rates Mk ancestral and edge joint probabilities at a fixed rate.

Scaled pruning extends EAzami's binary edge-joint calculation to nominal
multistate observations. Allowed tip sets encode uncertainty, not observed
within-species frequencies. Endpoint-state differences are not jump counts.
"""
import numpy as np


def posterior(tree, states, universe, q):
    universe=tuple(universe);k=len(universe)
    if k<2 or len(set(universe))!=k or not np.isfinite(q) or q<=0:
        raise ValueError('Invalid state universe or rate')
    nodes=list(tree.find_clades(order='preorder'));idx={n:i for i,n in enumerate(nodes)}
    tip_names=[n.name for n in tree.get_terminals()]
    if len(set(tip_names))!=len(tip_names) or not set(states)<=set(tip_names):
        raise ValueError('Duplicated or unmatched observed taxa')
    transitions={}
    for p in nodes:
        for c in p.clades:
            t=c.branch_length
            if t is None or not np.isfinite(t) or t<0:raise ValueError('Invalid branch length')
            e=np.exp(-k*q*t);transitions[c]=np.ones((k,k))*(1-e)/k+np.eye(k)*e
    up={};scales={};messages={}
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            allowed=set(states.get(node.name,universe))
            if not allowed or not allowed<=set(universe):raise ValueError('Invalid allowed states')
            up[node]=np.array([s in allowed for s in universe],float)/len(allowed);scales[node]=0.
        else:
            raw=np.ones(k);scale=0.
            for child in node.clades:
                messages[child]=transitions[child]@up[child]
                raw*=messages[child];scale+=scales[child]
            z=raw.sum()
            if z<=0:raise ValueError('Zero likelihood')
            up[node]=raw/z;scales[node]=scale+np.log(z)
    root=tree.root;down={root:np.ones(k)/k};edge_rows=[]
    for p in nodes:
        for c in p.clades:
            outside=down[p].copy()
            for sibling in p.clades:
                if sibling is not c:outside*=messages[sibling]
            outside/=outside.sum()
            down[c]=outside@transitions[c];down[c]/=down[c].sum()
            joint=outside[:,None]*transitions[c]*up[c][None,:];joint/=joint.sum()
            edge_rows.append(dict(parent=idx[p],child=idx[c],joint=joint.tolist(),
                descendant_taxa=sorted(t.name for t in c.get_terminals()),branch_length=float(c.branch_length)))
    node_rows=[]
    for n in nodes:
        prob=up[n]*down[n];prob/=prob.sum()
        node_rows.append(dict(node=idx[n],name=n.name,is_tip=n.is_terminal(),
            descendant_taxa=sorted(t.name for t in n.get_terminals()),probabilities=prob.tolist()))
    return dict(universe=list(universe),q=float(q),log_likelihood=float(scales[root]+np.log(np.mean(up[root]))),nodes=node_rows,edges=edge_rows)
