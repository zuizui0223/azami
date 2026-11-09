"""Terminal-branch approximation to convergence, conditional on two nuclear trees."""
import itertools
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from Bio import Phylo
from statsmodels.stats.multitest import multipletests
def rooted_ingroup(tree, expected, keep):
    out="OUTGROUP_saff"
    assert {n.name for n in tree.get_terminals()}==set(expected)|{out}
    before=[tree.distance(a,b) for a,b in itertools.combinations(keep,2)]
    tree.root_with_outgroup(out)
    tree.prune(out)
    for n in list(tree.get_terminals()):
        if n.name not in keep:tree.prune(n)
    assert np.allclose(before,[tree.distance(a,b) for a,b in itertools.combinations(keep,2)])
    return tree

def maps(tree,names):
    nodes=list(tree.find_clades(order='preorder'))
    idx={node:i for i,node in enumerate(nodes)}
    tips=[idx[next(t for t in tree.get_terminals() if t.name==n)] for n in names]
    parent={idx[c]:idx[p] for p in nodes for c in p.clades}
    branch=np.array([0. if n is tree.root else max(n.branch_length,1e-10) for n in nodes])
    M=np.zeros((len(nodes),len(nodes)))
    for i in range(1,len(nodes)):
        M[i]=M[parent[i]];M[i,i]=1.
    V=(M*branch)@M.T
    K=V[np.ix_(tips,tips)]
    precision=np.linalg.inv(K)
    one=np.ones(len(tips))
    intercept=one@precision/(one@precision@one)
    U=V[:,tips]@precision
    # Estimated root mean plus Gaussian conditional mean, for every node.
    A=U+np.outer(1-U@one,intercept)
    assert np.allclose(A[tips],np.eye(len(tips)),atol=1e-5)
    return K,precision,intercept,A,np.array([parent[t] for t in tips])

def scores(tips,anc,i,j):
    return (np.linalg.norm(anc[...,i,:]-anc[...,j,:],axis=-1)-np.linalg.norm(tips[...,i,:]-tips[...,j,:],axis=-1)).mean(axis=-1)

def analyse(root):
    inp=root/'inputs';d=pd.read_csv(inp/'tip_means.csv')
    climate=['chelsa_bio12','chelsa_bio01','chelsa_rsds_mean']
    head=['orientation_image_vertical_angle','corolla_lab_chroma','capitulum_outline_aspect_ratio']
    for c in climate+head:d[c]=(d[c]-d[c].mean())/d[c].std(ddof=0)
    rows=[]
    for file in sorted(inp.glob('*.nwk')):
        for low in [True,False]:
            part=d if low else d[d.S1_source_taxon!='Cirsium esculentum']
            names=list(part.run)
            tr=rooted_ingroup(Phylo.read(file,'newick'),d.run,names)
            tr.root.branch_length=0.
            K,P,intercept,A,parent=maps(tr,names)
            pairs=np.array(list(itertools.combinations(range(len(names)),2)))
            env=part[climate].to_numpy()
            for threshold in [.9,1.1]:
                i,j=pairs.T
                keep=(parent[i]!=parent[j])&(np.linalg.norm(env[i]-env[j],axis=1)<=threshold)
                background=(parent[i]!=parent[j]) & ~keep
                bi,bj=i[background],j[background]
                i,j=i[keep],j[keep]
                for cols in [head[:2],head]:
                    if len(i)<3:continue
                    y=part[cols].to_numpy()
                    y=(y-y.mean(0))/y.std(0)
                    centered=y-intercept@y
                    rate=centered.T@P@centered/(len(y)-1)
                    rng=np.random.default_rng(20261109)
                    z=np.linalg.cholesky(K)@rng.standard_normal((4999,len(y),len(cols)))@np.linalg.cholesky(rate).T
                    z=(z-z.mean(axis=1,keepdims=True))/z.std(axis=1,keepdims=True)
                    obs=scores(y,A[parent]@y,i,j);null=scores(z,A[parent]@z,i,j)
                    background_obs=scores(y,A[parent]@y,bi,bj)
                    background_null=scores(z,A[parent]@z,bi,bj)
                    gap=obs-background_obs
                    rows.append(dict(tree='ASTRAL' if 'ASTRAL' in file.name else 'ML',n_taxa=len(y),threshold=threshold,head_dims=len(cols),pairs=len(i),
                                     taxa_in_pairs=len(np.unique(np.r_[i,j])),other_pairs=len(bi),contrast=gap,contrast_p=(1+((null-background_null)>=gap).sum())/5000,parent_minus_tip=obs,bm_mean=null.mean(),bm_p=(1+(null>=obs).sum())/5000))
    tab=pd.DataFrame(rows)
    primary=(tab.n_taxa==24)&(tab.threshold==.9)&(tab.head_dims==2)
    assert primary.sum()==2
    ids=tab.index[primary]
    corrected=multipletests(tab.loc[ids,['bm_p','contrast_p']].to_numpy().ravel(),method='holm')[1].reshape(-1,2)
    tab.loc[ids,['holm_approach','holm_contrast']]=corrected
    tab.to_csv(root/'terminal_approach.csv',index=False)
    print(tab.round(4).to_string(index=False))
if __name__=='__main__':analyse(Path(sys.argv[1]))
