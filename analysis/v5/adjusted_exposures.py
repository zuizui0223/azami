"""Compare biologically defined climate exposures on one fixed introduced cohort."""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

SPECS={
    "month_mean":["prec","tavg","srad"],
    "month_max":["prec","tmax","srad"],
    "previous_max":["prec","prior_tmax","srad"],
    "month_max_annual_background":["prec","tmax","srad","annual_prec","annual_tavg","annual_srad"],
    "month_mean_annual_background":["prec","tavg","srad","annual_prec","annual_tavg","annual_srad"],
}


def compare(sp,frame,monthly,month,out):
    from analysis.v5.record_density import fit_gradients
    d=frame.copy()
    idx=np.arange(len(d))
    for axis,table in monthly.items():
        d[axis]=table[idx,month]
        d["annual_"+axis]=table.mean(1)
    d["prior_tmax"]=monthly["tmax"][idx,(month-1)%12]
    needed=list(dict.fromkeys(col for axes in SPECS.values() for col in axes))
    d=d[d.novel.isin([0,1])].dropna(subset=[*needed,"year"]).copy()
    rows=[];diagnostics=[]
    for scope,g in d.groupby("novel"):
        w=1/g.groupby(["i","j"])["novel"].transform("size")
        x=g[needed].to_numpy()
        cov=np.cov(x.T,aweights=w)
        corr=cov/np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
        for j,a in enumerate(needed):
            for k in range(j):
                diagnostics.append(dict(species=sp,novel=int(scope),a=a,b=needed[k],r=corr[j,k]))
    for label,axes in SPECS.items():
        rows.extend(fit_gradients(sp,d,label,axes))
    coverage=[dict(species=sp,novel=int(k),records=len(g),cells=len(g[["i","j"]].drop_duplicates()),
                   blocks=g.block.nunique()) for k,g in d.groupby("novel")]
    for name,data in [("adjusted_slopes",rows),("adjusted_correlations",diagnostics),("adjusted_support",coverage)]:
        path=out/(name+".csv")
        pd.DataFrame(data).to_csv(path,index=False,mode="a",header=not path.exists())


def finalize(out):
    p=out/"adjusted_slopes.csv";d=pd.read_csv(p)
    d["holm_model"]=d.groupby(["climate","regional","term"]).p.transform(lambda s:multipletests(s,method="holm")[1])
    # In the three temperature substitutions, also account for inspecting all three specifications.
    primary=d.climate.isin(["month_mean","month_max","previous_max"])
    d.loc[primary,"holm_substitutions"]=d[primary].groupby(["regional","term"]).p.transform(lambda s:multipletests(s,method="holm")[1])
    d.to_csv(p,index=False)
