"""Decompose the fixed V5 climate-analogue distance test on identical draws."""
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v4.s1_climate_trait_hypervolume import load_source,geo_aggregate,standardise_cells
from analysis.v4 import s1_h3_climate_analogue_alignment as h

root=Path('work/bio18_geo')
def main():
    df,_=load_source(root/'private/decoded/observations.csv.gz',root/'private/decoded/measurement_receipt.json')
    records=[]
    for deg in [1,5]:
        cells,_=geo_aggregate(df,deg);z,_=standardise_cells(cells);taxa=h.eligible_taxa(z)
        rng=np.random.default_rng(20261109);streams=np.random.default_rng(20261110)
        estimates=[]; nulls=[]
        for rep in range(24):
            x,y,sector=h._sample_taxa(taxa,rng)
            pairs,support=h.matched_pairs(x,y,sector)
            if support['status']!='qualified':continue
            differences=y[pairs['i'],:,None,:]-y[pairs['j'],None,:,:]
            total=pairs['d']
            contributions=np.divide(differences**2,total[...,None],out=np.zeros_like(differences),where=total[...,None]>0)
            assert np.allclose(contributions.sum(axis=3),total)
            # Identical shuffles across components and metrics.
            seed=int(streams.integers(2**31));obs=[];null=[]
            for metric in [total,*np.moveaxis(contributions,-1,0),*np.moveaxis(abs(differences),-1,0)]:
                score,perm=h.observed_and_permutations(dict(pairs,d=metric),len(taxa),np.random.default_rng(seed),499)
                obs.append(score);null.append(perm)
            assert np.allclose(sum(null[1:4]),null[0])
            estimates.append(obs);nulls.append(null)
        assert len(estimates)>=12
        a=np.mean(estimates,axis=0);b=np.mean(nulls,axis=0)
        for j,term in enumerate(['joint','angle_part','chroma_part','outline_part','angle_alone','chroma_alone','outline_alone']):
            null=b[j];excess=a[j]-null.mean()
            p=(1+np.count_nonzero(abs(null-null.mean())>=abs(excess)-1e-12))/(1+len(null))
            records.append(dict(degrees=deg,taxa=len(taxa),draws=len(estimates),metric=term,observed=a[j],null_mean=null.mean(),excess=excess,p=p))
    tab=pd.DataFrame(records)
    for typ in ['part','alone','joint']:
        test=tab.metric.str.endswith(typ)
        tab.loc[test,'holm_family']=multipletests(tab.loc[test,'p'],method='holm')[1]
    tab.to_csv(root/'public/v5_trait_convergence.csv',index=False)
if __name__=='__main__':main()
