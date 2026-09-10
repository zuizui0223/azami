"""Pairwise density overlap and exact empirical W2 on frozen pilot coordinates."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from run_distribution_breadth_pilot import ENV, log_kde, rng_for


def wasserstein_parts(x,y):
    # Equal empirical weights: optimal assignment gives exact finite-sample W2.
    if x.shape != y.shape or x.ndim != 2 or len(x)==0 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('W2 assignment requires equal nonempty finite sample shapes')
    cost=cdist(x,y,'sqeuclidean')
    i,j=linear_sum_assignment(cost)
    total=float(cost[i,j].mean())
    centre=float(np.sum((x.mean(0)-y.mean(0))**2))
    assert total>=centre-1e-9
    return np.sqrt(total),np.sqrt(centre),np.sqrt(max(0,total-centre))


def overlap(x,y,h,draw_x,draw_y,own_x,own_y):
    # Importance sampling from the equal KDE mixture, estimating integral min(p,q).
    a=log_kde(draw_x,y,h); b=log_kde(draw_y,x,h)
    values=np.r_[2*np.exp(np.minimum(own_x,a)-np.logaddexp(own_x,a)),
                 2*np.exp(np.minimum(own_y,b)-np.logaddexp(own_y,b))]
    return float(values.mean())


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--traits',required=True,type=Path)
    p.add_argument('--environment',required=True,type=Path)
    p.add_argument('--basis-dir',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--replicates',type=int,default=30)
    p.add_argument('--draws-per-taxon',type=int,default=128)
    p.add_argument('--seed',type=int,default=20260910)
    args=p.parse_args(); args.out.mkdir(parents=True,exist_ok=False)
    old=json.loads((args.basis_dir/'contract.json').read_text())
    for key,path in [('traits',args.traits),('environment',args.environment)]:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==old['sha256'][key]
    bases=json.loads((args.basis_dir/'shared_bases.json').read_text())
    contract=dict(status='exploratory_pilot',source_sha256=old['sha256'],
                  basis_sha256=hashlib.sha256((args.basis_dir/'shared_bases.json').read_bytes()).hexdigest(),
                  sample_size=50,replicates=args.replicates,draws_per_taxon=args.draws_per_taxon,seed=args.seed,
                  native_filter=False,new_images=False,
                  overlap='integral min(p,q), Monte Carlo from equal KDE mixture; not geometric intersection/union',
                  wasserstein='exact W2 assignment between equal-weight 50-point empirical distributions; not KDE distance',
                  decomposition='W2 squared = centroid distance squared + centred distribution W2 squared',
                  distinctiveness='equal-weight mean 1-overlap to other eligible taxa, conditional on sampled reference set',
                  limitations=['shared fixed KDE bandwidth can inflate overlap','finite-sample W2 positive even for same population',
                               'rarefaction ranges are not confidence intervals','same projection limitations as first pilot',
                               'no biological equivalence or convergence inference','no tests treating pairs as independent'])
    (args.out/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    t=pd.read_csv(args.traits,usecols=['obs_id','taxon_name','endpoint_id','value','measurement_available'],dtype={'obs_id':str})
    t=t[t.measurement_available.astype(str).str.lower().isin(['true','1']) & np.isfinite(t.value)]
    wide=t.pivot(index=['obs_id','taxon_name'],columns='endpoint_id',values='value').reset_index()
    e=pd.read_csv(args.environment,dtype={'obs_id':str}).replace([np.inf,-np.inf],np.nan).dropna(subset=[*ENV,'latitude','longitude'])
    summaries=[]; all_pairs=[]; distinct=[]
    for basis in bases:
        group=basis['group']; members=basis['members']; h=basis['bandwidth']; d=basis['dimensions']
        f=wide[['obs_id','taxon_name',*members]].dropna().merge(e[['obs_id','taxon_name']],on=['obs_id','taxon_name'],validate='one_to_one')
        scores=((f[members].to_numpy()-basis['centre'])/basis['scale'])@np.array(basis['loadings'])
        f=f.reset_index(drop=True)
        samples={taxon:scores[np.array(list(idx))] for taxon,idx in f.groupby('taxon_name').groups.items() if len(idx)>=50}
        taxa=sorted(samples); assert len(taxa)==basis['eligible_taxa_50']
        pairs=[(a,b) for i,a in enumerate(taxa) for b in taxa[i+1:]]
        metrics=np.empty((args.replicates,len(pairs),4))
        for rep in range(args.replicates):
            cache={}
            for taxon in taxa:
                r=rng_for(args.seed,group,taxon,rep,'overlap')
                x=samples[taxon][r.choice(len(samples[taxon]),50,replace=False)]
                draw=x[r.integers(50,size=args.draws_per_taxon)]+r.normal(size=(args.draws_per_taxon,d))*h
                cache[taxon]=(x,draw,log_kde(draw,x,h))
            for j,(a,b) in enumerate(pairs):
                x,dx,px=cache[a];y,dy,py=cache[b]
                ov=overlap(x,y,h,dx,dy,px,py)
                metrics[rep,j]=[ov,*wasserstein_parts(x,y)]
            if (rep+1)%10==0: print(f'{group}: {rep+1}/{args.replicates} repetitions, {len(pairs)} pairs',flush=True)
        assert np.isfinite(metrics).all() and np.all((metrics[:,:,0]>=0)&(metrics[:,:,0]<=1))
        np.savez_compressed(args.out/f'{group}_replicates.npz',metrics=metrics,taxa=np.array(taxa))
        rows=[]
        for j,(a,b) in enumerate(pairs):
            row=dict(group=group,taxon_a=a,taxon_b=b)
            for k,name in enumerate(['density_overlap','wasserstein_w2','centroid_distance','centred_w2']):
                row[name]=float(np.median(metrics[:,j,k]))
                row[name+'_rarefaction_p025']=float(np.quantile(metrics[:,j,k],.025))
                row[name+'_rarefaction_p975']=float(np.quantile(metrics[:,j,k],.975))
            # Calculate proportions per replicate before summary; medians do not obey the identity.
            ratio=np.divide(metrics[:,j,3]**2,metrics[:,j,1]**2,out=np.zeros(args.replicates),where=metrics[:,j,1]>0)
            row['centred_share_w2_squared']=float(np.median(ratio))
            rows.append(row)
        df=pd.DataFrame(rows);all_pairs.extend(rows)
        for taxon in taxa:
            selected=df[(df.taxon_a==taxon)|(df.taxon_b==taxon)]
            distinct.append(dict(group=group,taxon_name=taxon,reference_taxa=len(taxa)-1,
                                 mean_density_distinctiveness=float((1-selected.density_overlap).mean())))
        summaries.append(dict(group=group,taxa=len(taxa),pairs=len(pairs),dimensions=d,retained_variance=basis['retained_variance'],
                              median_pair_density_overlap=float(df.density_overlap.median()),
                              median_pair_w2=float(df.wasserstein_w2.median()),
                              median_pair_centred_share=float(df.centred_share_w2_squared.median()),
                              interpretation='descriptive_projection_only' if group in ['involucre','whole'] else 'exploratory'))
        print(json.dumps(summaries[-1]),flush=True)
    pd.DataFrame(all_pairs).to_csv(args.out/'pairwise_distributions.csv',index=False)
    pd.DataFrame(distinct).to_csv(args.out/'taxon_distinctiveness.csv',index=False)
    (args.out/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')


if __name__=='__main__':main()
