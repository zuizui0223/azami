"""Join exact frozen S1 aggregate cohorts to independently inferred nuclear trees."""
import copy
import math
import numpy as np
import pandas as pd
from analysis.v4 import s1s2_multitrait_branch_core as core
from analysis.v4.conditional_ancestry_export import reconstruct


def evaluate_unit(tree,means,columns,*,nperm=9999,seed=20261025):
    names=list(means.index)
    D,internal=core.branch_operator(tree,names)
    trait=means[list(columns)].to_numpy(float)
    environment=means[list(core.FOUR)].to_numpy(float)
    ancestry=reconstruct(tree,names,np.column_stack([trait,environment]),[*columns,*core.FOUR])
    edge_values=np.array([[e['delta'][c] for c in columns] for e in ancestry['edges']])
    if not np.allclose(edge_values,D@trait,atol=1e-8):
        raise ValueError('Ancestral export and branch-test operator disagree')
    scopes={}
    for scope,operator in [('all_edges',D),('internal_edges',D[internal])]:
        if len(operator)<3:
            scopes[scope]={'status':'INSUFFICIENT_INTERNAL_EDGES','edge_count':len(operator)}
        else:
            scopes[scope]={'status':'COMPUTED_CONDITIONAL_TIP_SHUFFLES','edge_count':len(operator),
                'axes':core.alignment_test(operator,environment,trait,nperm=nperm,seed=seed)}
    return {'ancestry':ancestry,'scopes':scopes}


def join_moments(moments,manifest,cohort,columns):
    """Require the same photographs for every trait/environment coordinate."""
    required={'cohort','accepted_key','coordinate','n_observations','mean','sample_sd'}
    if not required<=set(moments) or not {'run','accepted_key'}<=set(manifest):
        raise ValueError('Missing aggregate or manifest fields')
    mapping=manifest.copy();mapping['accepted_key']=mapping.accepted_key.astype(str)
    if mapping.run.duplicated().any() or mapping.accepted_key.duplicated().any():
        raise ValueError('One-to-one run and accepted-key manifest required')
    data=moments.copy();data['accepted_key']=data.accepted_key.astype(str)
    data=data.loc[(data.cohort==cohort)&data.accepted_key.isin(mapping.accepted_key)&data.coordinate.isin(columns)]
    if not columns or len(set(columns))!=len(columns) or data.duplicated(['accepted_key','coordinate']).any():
        raise ValueError('Missing or duplicated coordinate identity')
    records=[]
    for row in mapping.itertuples():
        block=data.loc[data.accepted_key==row.accepted_key].set_index('coordinate')
        if set(block.index)!=set(columns):raise ValueError('Incomplete exact matched taxon')
        numeric=block[['n_observations','mean','sample_sd']].to_numpy(float)
        counts=block.n_observations.to_numpy(float)
        if (not np.isfinite(numeric).all() or (block.sample_sd<0).any()
                or not np.equal(counts,np.floor(counts)).all() or (counts<5).any()
                or len(set(counts))!=1):
            raise ValueError('Different photo cohorts or invalid aggregate moments')
        records.append(dict(run=row.run,n_original_observations=int(counts[0]),
                            **{c:float(block.loc[c,'mean']) for c in columns}))
    return pd.DataFrame(records).set_index('run')


def rooted_ingroup(source,manifest_runs,retained_runs):
    """Root on the observed saff tip, then prune; do not midpoint-root the ingroup."""
    all_runs=list(manifest_runs);keep=list(retained_runs)
    tips=[n.name for n in source.get_terminals()]
    if (len(set(tips))!=len(tips) or len(set(all_runs))!=len(all_runs)
            or set(tips)!=set(all_runs)|{'OUTGROUP_saff'} or not keep
            or len(set(keep))!=len(keep) or not set(keep)<=set(all_runs)):
        raise ValueError('Tree and explicit run identities differ')
    for node in source.find_clades():
        if node is source.root:continue
        if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length<0:
            raise ValueError('Invalid molecular branch length')
    expected={(a,b):source.distance(a,b) for a in keep for b in keep}
    tree=copy.deepcopy(source);tree.root_with_outgroup('OUTGROUP_saff')
    tree.prune('OUTGROUP_saff')
    for tip in list(tree.get_terminals()):
        if tip.name not in keep:tree.prune(tip)
    if set(n.name for n in tree.get_terminals())!=set(keep):raise ValueError('Pruning identity changed')
    if any(not math.isclose(tree.distance(a,b),v,rel_tol=1e-10,abs_tol=1e-12) for (a,b),v in expected.items()):
        raise ValueError('Pruning changed retained pairwise distances')
    return tree
