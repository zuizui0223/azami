"""Keep observed within-cell variation with equal total weight per occupied cell."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.ndimage import gaussian_filter
from statsmodels.stats.multitest import multipletests
from analysis.v5.novel_head_overlap import T,make_cells
from analysis.v5.seasonal_exposure import normals


def observations(frame, cell_table):
    d=frame[frame[T[0]].between(0,180)].dropna(subset=T).copy()
    d["i"]=np.minimum(179,np.floor(d.latitude+90)).astype(int)
    d["j"]=np.floor((d.longitude+180)%360).astype(int)
    d=d.merge(cell_table[["i","j","novel","block"]],on=["i","j"],suffixes=("_record",""),validate="many_to_one")
    d["w"]=1/d.groupby(["i","j"])[T[0]].transform("size")
    np.testing.assert_allclose(d.w.sum(),len(cell_table))
    assert len(d)==cell_table.photos.sum()
    return d


def density(z,w,bounds,h):
    edges=[np.linspace(lo-3*h,hi+3*h,41) for lo,hi in bounds.T]
    hist=np.histogramdd(np.clip(z,*bounds),bins=edges,weights=w)[0]
    smooth=gaussian_filter(hist,[h/np.diff(e)[0] for e in edges],mode="constant")
    return smooth/smooth.sum()


def gradients(sp,d,folder):
    date=pd.to_datetime(d.observed_on,format="mixed",errors="coerce")
    month=date.dt.month.fillna(1).astype(int).to_numpy()-1
    axes=["prec","tavg","srad"]
    monthly=normals(d,folder)
    d=d.copy()
    for axis in axes:
        d[axis]=monthly[axis][np.arange(len(d)),month]
    d["sin"]=np.sin(2*np.pi*month/12)
    d["cos"]=np.cos(2*np.pi*month/12)
    d["year"]=date.dt.year-2020
    d=d[d.novel.isin([0,1])].dropna(subset=[*axes,"year"]).copy()
    models={"head_month":axes,"niche":["bio01","bio04","bio12","bio15","rsds_mean"]}
    result=[]
    for model,axes in models.items():
        result.extend(fit_gradients(sp,d,model,axes))
    return result


def fit_gradients(sp,d,model,axes):
    w=1/d.groupby(["i","j"])[T[0]].transform("size").to_numpy()
    def standardized(v):
        mean=np.average(v,axis=0,weights=w)
        sd=np.sqrt(np.average((v-mean)**2,axis=0,weights=w))
        return (v-mean)/sd,sd
    x,climate_sd=standardized(d[axes].to_numpy())
    y,trait_sd=standardized(d[T].to_numpy())
    new=d.novel.to_numpy()
    raw=np.column_stack([x,x*new[:,None],d[["sin","cos","year"]]])
    rows=[]
    for regional in [False,True]:
        # Separate intercepts for every region x novelty category in sensitivity.
        group=(d.block*2+new if regional else new)
        group=pd.Series(group,index=d.index)
        z=pd.DataFrame(np.column_stack([raw,y]),index=d.index)
        avg=z.mul(w,axis=0).groupby(group).transform("sum").div(pd.Series(w,index=d.index).groupby(group).transform("sum"),axis=0)
        centered=(z-avg).to_numpy()
        X=centered[:,:raw.shape[1]]
        assert np.linalg.matrix_rank(X)==X.shape[1]
        loo=[]
        err=np.zeros(len(T)); baseline_error=np.zeros(len(T))
        novel=X[new==1,:len(axes)]
        time=X[new==1,-3:]
        nw=np.sqrt(w[new==1])[:,None]
        residual=novel*nw-(time*nw)@np.linalg.lstsq(time*nw,novel*nw,rcond=None)[0]
        # Residual columns have weighted zero means before sqrt weighting;
        # use cross products to preserve the intended WLS metric exactly.
        gram=residual.T@residual
        corr=gram/np.sqrt(np.outer(np.diag(gram),np.diag(gram)))
        vif=np.diag(np.linalg.inv(corr))
        if regional:
            for block in np.unique(d.block):
                keep=d.block.to_numpy()!=block
                wx=np.sqrt(w[keep])[:,None]
                xx=X[keep]*wx
                if np.linalg.matrix_rank(xx)==X.shape[1]:
                    yy=centered[keep,raw.shape[1]:]*wx
                    fitted=np.linalg.lstsq(xx,yy,rcond=None)[0]
                    loo.append(fitted)
                    te=(~keep)&(new==1)
                    if te.any():
                        truth=centered[te,raw.shape[1]:]
                        time_fit=np.linalg.lstsq(xx[:,-3:],yy,rcond=None)[0]
                        err+=np.sum(w[te,None]*(truth-X[te]@fitted)**2,axis=0)
                        baseline_error+=np.sum(w[te,None]*(truth-X[te,-3:]@time_fit)**2,axis=0)
        for j,trait in enumerate(["angle","chroma","outline"]):
            f=sm.WLS(centered[:,raw.shape[1]+j],X,weights=w).fit(cov_type="cluster",cov_kwds={"groups":d.block},use_t=True)
            for k,axis in enumerate(axes):
                for term,indices in [("analogue",[k]),("novel",[k,k+len(axes)]),("difference",[k+len(axes)])]:
                    con=np.zeros(X.shape[1]);con[indices]=1
                    t=f.t_test(con);b=t.effect.item();lo,hi=t.conf_int().ravel()
                    rows.append(dict(species=sp,climate=model,regional=regional,trait=trait,axis=axis,term=term,beta=b,
                                     lo=lo,hi=hi,p=t.pvalue.item(),
                                     loo_min=min((con@b0[:,j] for b0 in loo),default=np.nan),
                                     loo_max=max((con@b0[:,j] for b0 in loo),default=np.nan),climate_sd=climate_sd[k],trait_sd=trait_sd[j],
                                     novel_VIF=vif[k],novel_within_region_cv_q2=1-err[j]/baseline_error[j] if regional else np.nan,
                                     native_like_cells=round(w[new==0].sum()),novel_cells=round(w[new==1].sum()),
                                     novel_blocks=d.loc[new==1,"block"].nunique()))
    return rows


def run(sp,native,intro,labeled,out):
    cells=[make_cells(native,True),make_cells(intro,True),make_cells(labeled)]
    data=[observations(g,c) for g,c in zip([native,intro,labeled],cells)]
    groups=[*data[:2],data[2][data[2].novel==0],data[2][data[2].novel==1]]
    cell_groups=[*cells[:2],cells[2][cells[2].novel==0],cells[2][cells[2].novel==1]]
    rows=[]
    for cols in [T,T[:2],*[T[k:k+1] for k in range(3)]]:
        mean=cells[0][cols].mean().to_numpy()
        sd=cells[0][cols].std(ddof=0).to_numpy()
        z=[(g[cols].to_numpy()-mean)/sd for g in groups]
        v=[(g[cols].to_numpy()-mean)/sd for g in cell_groups]
        for limits in ["original","all_values"]:
            pooled=np.vstack(v[:2] if limits=="original" else z[:2])
            bounds=np.quantile(pooled,[.005,.995] if limits=="original" else [0,1],axis=0)
            for h in [.35,.5,.75]:
                for representation,values,weights in [("records",z,[g.w for g in groups]),("medians",v,[None]*4)]:
                    densities=[density(m,w,bounds,h) for m,w in zip(values,weights)]
                    for k,group in enumerate(["all","analogue","novel"],1):
                        rows.append(dict(species=sp,traits="+".join(cols),dimensions=len(cols),limits=limits,h=h,
                                         representation=representation,group=group,D=np.minimum(densities[0],densities[k]).sum(),
                                         cells=len(cell_groups[k]),records=len(groups[k])))
    for name,records in [("record_density",rows),("novel_gradient",gradients(sp,data[2],out.parent/"temp"))]:
        path=out/(name+".csv")
        pd.DataFrame(records).to_csv(path,index=False,mode="a",header=not path.exists())


def finalize(out):
    p=out/"novel_gradient.csv"
    d=pd.read_csv(p)
    d["holm"]=d.groupby(["climate","regional","term"]).p.transform(lambda x:multipletests(x,method="holm")[1])
    d.to_csv(p,index=False)
