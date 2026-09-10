"""Common-cohort omnibus environment test; equal weight per biological construct."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v3 import run_biological_axis_reanalysis as axis
from analysis.v3.run_construct_scale_upgrade import common_complete_features,CORE


def statistic(beta,groups):
    return np.mean([np.mean(beta[...,idx]**2,axis=-1) for idx in groups],axis=0)


def test_gradient(y,x,groups,blocks,perms,rng):
    den=x@x;beta=x@y/den;observed=float(statistic(beta,groups));exceed=0
    for start in range(0,perms,128):
        size=min(128,perms-start);xp=np.broadcast_to(x,(size,len(x))).copy()
        for idx in blocks:
            xp[:,idx]=rng.permuted(xp[:,idx],axis=1)
        null=statistic(xp@y/den,groups);exceed+=int(np.sum(null>=observed-1e-15))
    return beta,observed,(exceed+1)/(perms+1)


def main():
    p=argparse.ArgumentParser();p.add_argument('--traits',type=Path,required=True);p.add_argument('--environment',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--permutations',type=int,default=9999);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    contract={'base_commit':'4deae0815880baf5db3e2795bd79369344e8bace','constructs':CORE,'predictors':axis.PREDICTORS,
      'cohort':'exact complete18 minimum5 cohort; same observations across all constructs and complete nine environments',
      'among':'common-cohort taxon medians of construct coordinates and environmental variables, one row per taxon',
      'within':'taxon demeaned standardized coordinates and predictor; observation weighted as original within atlas',
      'statistic':'mean across 9 constructs of mean squared standardized component slopes; joint constructs do not gain weight merely from more dimensions',
      'permutation':'permute predictor once per replicate, shared across all response coordinates; unrestricted among taxa, within-taxon blocks for observations',
      'permutations':a.permutations,'seed':20260910,'multiplicity':'one new separate BH family of 18 omnibus tests; original PR93 families untouched',
      'contribution':'component mean-square slope / sum across constructs; descriptive, not variance partition or separately tested effects',
      'limits':['omnibus association may be driven by one construct, not coordinated response','new rows have not passed spatial or phylogenetic sensitivity','within permutation assumes exchangeability within taxa, not spatial independence','complete-case cohort has restricted support; not full genus','marginal correlated gradients, not causal effects']}
    contract['input_sha256']={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in [('traits',a.traits),('environment',a.environment)]}
    assert contract['input_sha256']['traits']=='d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344'
    assert contract['input_sha256']['environment']=='e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7'
    (a.out/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    t,e=axis.load_data(a.traits,a.environment);features,_=common_complete_features(t,5);fmap=features.attrs['feature_map']
    f=features.merge(e[['obs_id','taxon_name',*axis.PREDICTORS]],on=['obs_id','taxon_name'],validate='one_to_one').dropna()
    assert len(f)==1734 and f.taxon_name.nunique()==42
    cols=sum([fmap[c] for c in CORE],[]);groups=[np.array([cols.index(col) for col in fmap[c]]) for c in CORE]
    rows=[];contributions=[]
    for scale in ['among_taxon','within_taxon']:
        if scale=='among_taxon':
            data=f.groupby('taxon_name')[cols+axis.PREDICTORS].median();y=np.column_stack([axis.z(data[col]) for col in cols]);blocks=[np.arange(len(data))]
        else:
            data=f.copy();y=np.column_stack([axis.demean_std(data,col) for col in cols]);blocks=[np.flatnonzero(data.taxon_name.to_numpy()==taxon) for taxon in sorted(data.taxon_name.unique())]
        for predictor in axis.PREDICTORS:
            x=axis.z(data[predictor]) if scale=='among_taxon' else axis.demean_std(data,predictor)
            beta,value,pv=test_gradient(y,x,groups,blocks,a.permutations,axis.stable_rng(20260910,scale,predictor,'whole_environment'))
            scores=np.array([np.mean(beta[idx]**2) for idx in groups]);shares=scores/scores.sum()
            rows.append(dict(scale=scale,predictor=predictor,n=len(data),taxa=42,statistic=value,p_value=pv,largest_contributor=CORE[int(scores.argmax())],largest_share=float(shares.max())))
            for c,idx,score,share in zip(CORE,groups,scores,shares):
                contributions.append(dict(scale=scale,predictor=predictor,construct=c,mean_squared_slope=score,statistic_share=share,component_slopes=json.dumps(beta[idx].tolist())))
            print(scale,predictor,value,pv,flush=True)
    result=pd.DataFrame(rows);result['q_bh_18']=axis.bh(result.p_value);result['fdr_0_05']=result.q_bh_18<.05
    result.to_csv(a.out/'omnibus_environment.csv',index=False);pd.DataFrame(contributions).to_csv(a.out/'construct_contributions.csv',index=False)
    (a.out/'summary.json').write_text(json.dumps({'observations':1734,'taxa':42,'tests':18,'fdr_supported':int(result.fdr_0_05.sum()),'hits':result[result.fdr_0_05].to_dict('records'),'robustness':'not_yet_spatial_or_phylogenetic_tested'},indent=2)+'\n')

if __name__=='__main__':main()
