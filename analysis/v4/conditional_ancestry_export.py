"""Conditional Brownian node estimates in original coordinate units.

Observed tips are fixed, the root has a flat prior, and each coordinate has
a plug-in REML diffusion estimate. Marginal normal intervals condition on
that estimate and one tree; they are not simultaneous confidence intervals,
measurement-error intervals, historical climate observations or event counts.
Hue sine/cosine stay coordinates, never independent inferred pigment states.
"""
from __future__ import annotations
import hashlib
import numpy as np
from scipy.linalg import solve


def reconstruct(tree, tip_order, values, columns):
    x=np.asarray(values,dtype=float)
    tips=tree.get_terminals()
    if (len(tip_order)<8 or len(set(tip_order))!=len(tip_order)
            or len({t.name for t in tips})!=len(tips)
            or {t.name for t in tips}!=set(tip_order)):
        raise ValueError('Require eight or more unique exactly matched tips')
    if (x.shape!=(len(tip_order),len(columns)) or not columns
            or len(set(columns))!=len(columns) or not np.isfinite(x).all()):
        raise ValueError('Invalid observed coordinate matrix')
    nodes=list(tree.find_clades(order='preorder')); ix={n:i for i,n in enumerate(nodes)}
    L=np.zeros((len(nodes),len(nodes))); edges=[]; regularized=0
    for p in nodes:
        for c in p.clades:
            length=c.branch_length
            if length is None or not np.isfinite(length) or length<0:
                raise ValueError('Invalid source branch length')
            regularized+=int(length==0)
            effective=max(float(length),1e-10)
            a,b=ix[p],ix[c];w=1/effective
            L[a,a]+=w;L[b,b]+=w;L[a,b]-=w;L[b,a]-=w
            edges.append((a,b,effective,float(length)))
    by={n.name:n for n in tips}
    fixed=[ix[by[name]] for name in tip_order]
    free=[ix[n] for n in nodes if not n.is_terminal()]
    K=L[np.ix_(free,free)]
    inverse=solve(K,np.eye(len(free)),assume_a='pos')
    fitted=np.zeros((len(nodes),x.shape[1]));fitted[fixed]=x
    fitted[free]=solve(K,-L[np.ix_(free,fixed)]@x,assume_a='pos')
    delta=np.array([fitted[b]-fitted[a] for a,b,_,_ in edges])
    rate=np.sum(delta**2/np.array([e[2] for e in edges])[:,None],axis=0)/(len(tips)-1)
    variances=np.zeros_like(fitted)
    variances[free]=np.maximum(np.diag(inverse),0)[:,None]*rate
    se=np.sqrt(variances)
    conditional_covariance=np.zeros((len(nodes),len(nodes)))
    conditional_covariance[np.ix_(free,free)]=inverse
    edge_scale=np.array([conditional_covariance[a,a]+conditional_covariance[b,b]-2*conditional_covariance[a,b] for a,b,_,_ in edges])
    if np.min(edge_scale)<-1e-10:raise ValueError('Negative conditional edge variance')
    edge_se=np.sqrt(np.maximum(edge_scale,0)[:,None]*rate)
    if not all(np.isfinite(v).all() for v in (fitted,se,rate)):
        raise ValueError('Nonfinite Brownian solution')
    descendants=[sorted(t.name for t in n.get_terminals()) for n in nodes]
    ids=[f"node_{i:04d}_"+hashlib.sha256('\n'.join(d).encode()).hexdigest()[:12]
         for i,d in enumerate(descendants)]
    def record(v):return dict(zip(columns,map(float,v)))
    return {
        'model':'fixed_tip_flat_root_Brownian_plugin_REML_marginals',
        'uncertainty':'95% marginal normal intervals conditional on one tree and fitted diffusion; no simultaneous coverage',
        'measurement_error_included':False,'tree_uncertainty_included':False,
        'diffusion_estimation_uncertainty_included':False,
        'tip_zero_se_means_conditioned_not_error_free_measurement':True,
        'historical_environment_observed':False,
        'diffusion_rate':record(rate),'diffusion_degrees_of_freedom':len(tips)-1,
        'zero_diffusion_coordinates':[c for c,r in zip(columns,rate) if r==0],
        'zero_length_edges_regularized':regularized,'branch_length_floor':1e-10,
        'nodes':[{'node_id':ids[i],'is_tip':n.is_terminal(),
                  'source_taxon':n.name if n.is_terminal() else None,
                  'descendant_source_taxa':descendants[i],
                  'estimate':record(fitted[i]),'conditional_se':record(se[i]),
                  'conditional_lower_95':record(fitted[i]-1.959963984540054*se[i]),
                  'conditional_upper_95':record(fitted[i]+1.959963984540054*se[i])}
                 for i,n in enumerate(nodes)],
        'edges':[{'parent_id':ids[a],'child_id':ids[b],
                  'internal':not nodes[b].is_terminal(),'source_branch_length':length,
                  'delta':record(delta[j]),'conditional_delta_se':record(edge_se[j]),
                  'conditional_delta_lower_95':record(delta[j]-1.959963984540054*edge_se[j]),
                  'conditional_delta_upper_95':record(delta[j]+1.959963984540054*edge_se[j])}
                 for j,(a,b,_,length) in enumerate(edges)],
    }
