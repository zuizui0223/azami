"""Native-only predictions for the same introduced cell pairs as head overlap."""
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from analysis.v5.novel_head_overlap import make_cells, match_cells, T
from analysis.v5.rain_transfer import resample_blocks
from analysis.v5.seasonal_exposure import normals

MODELS = {"annual":["bio12","bio01","rsds_mean"],
          "head_month":["month_prec","month_tavg","month_srad"],
          "head_month_hot":["month_prec","month_tmax","month_srad"]}


def cells(frame, native, folder):
    d = frame.rename(columns=lambda c:c.removeprefix("chelsa_")).copy()
    d = d[d[T[0]].between(0,180)].dropna(subset=T)
    date = pd.to_datetime(d.observed_on,format="mixed",errors="coerce")
    monthly = normals(d,folder)
    month = date.dt.month.fillna(1).astype(int).to_numpy()-1
    for key,values in monthly.items():
        d["month_"+key] = np.where(date.notna(), values[np.arange(len(month)),month],np.nan)
    base = make_cells(d,native)
    d["i"] = np.minimum(179,np.floor(d.latitude+90)).astype(int)
    d["j"] = np.floor((d.longitude+180)%360).astype(int)
    axes = list(dict.fromkeys(a for model in MODELS.values() for a in model))
    return base.merge(d.groupby(["i","j"])[axes].median(),on=["i","j"],validate="one_to_one")


def coefficients(x,y):
    return np.linalg.lstsq(np.column_stack([np.ones(len(x)),x]),y,rcond=None)[0]


def run(sp,native,intro,out):
    n=cells(native,True,out.parent/"temp")
    c=cells(intro,False,out.parent/"temp")
    c=c[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c)
    a=c.iloc[np.concatenate([i for i,j in pairs])]
    b=c.iloc[np.concatenate([j for i,j in pairs])]
    rows=[]; shifts=[]; climate=[]; joint=[]
    for model,axes in MODELS.items():
        ref=n.dropna(subset=axes)
        ok=a[axes].notna().all(axis=1).to_numpy() & b[axes].notna().all(axis=1).to_numpy()
        left,right=a[ok],b[ok]
        x=ref[axes].to_numpy()
        mean=x.mean(0);sd=x.std(0)
        x=(x-mean)/sd
        y=ref[T].to_numpy()
        dx=(right[axes].to_numpy()-left[axes].to_numpy())/sd
        dy=right[T].to_numpy()-left[T].to_numpy()
        block=ref.block.to_numpy(); paired=left.block.to_numpy()
        fit=coefficients(x,y)
        prediction=dx@fit[1:]
        if model=="head_month":
            from analysis.v5.transfer_support import audit
            support, robustness=audit(sp,ref,left,right,axes,x,dx,dy,prediction,y[:,:2].std(0))
        cv=np.zeros_like(y);base=np.zeros_like(y)
        for tr,te in GroupKFold(5).split(x,groups=block):
            cv[te]=np.column_stack([np.ones(len(te)),x[te]])@coefficients(x[tr],y[tr])
            base[te]=y[tr].mean(0)
        def skill(obs,pred):
            return 1-((obs-pred)**2).sum(0)/(obs**2).sum(0)
        rng=np.random.default_rng(20261109); boot=[]
        sd_y=y[:,:2].std(0)
        def alignment(obs,pred):
            u,v=obs[:,:2]/sd_y,pred[:,:2]/sd_y
            length=np.linalg.norm(u,axis=1)*np.linalg.norm(v,axis=1)
            use=length>1e-12
            cos=(u[use]*v[use]).sum(1)/length[use]
            return [cos.mean(),np.mean(cos>0),1-np.sum((u-v)**2)/np.sum(u**2),len(cos)]
        directions=[]
        for _ in range(999):
            i=resample_blocks(block,rng); j=resample_blocks(paired,rng)
            p=dx[j]@coefficients(x[i],y[i])[1:]
            boot.append([dy[j].mean(0),p.mean(0),skill(dy[j],p)])
            directions.append(alignment(dy[j],p))
        result=alignment(dy,prediction)
        limits=np.quantile(directions,[.025,.975],axis=0)
        record=dict(species=sp,model=model,pairs=len(left),nonzero_pairs=result[3])
        for k,metric in enumerate(["mean_cosine","fraction_same_half_plane","joint_q2"]):
            record.update({metric:result[k],metric+"_lo":limits[0,k],metric+"_hi":limits[1,k]})
        joint.append(record)
        lo,hi=np.quantile(boot,[.025,.975],axis=0)
        for t,trait in enumerate(["angle","chroma","outline"]):
            rows.append(dict(species=sp,model=model,trait=trait,native=len(ref),pairs=len(left),blocks=len(set(paired)),
                native_cv_q2=1-((cv-y)**2).sum(0)[t]/((base-y)**2).sum(0)[t],
                observed=dy[:,t].mean(),obs_lo=lo[0,t],obs_hi=hi[0,t],predicted=prediction[:,t].mean(),
                pred_lo=lo[1,t],pred_hi=hi[1,t],q2_vs_no_shift=skill(dy,prediction)[t],q2_lo=lo[2,t],q2_hi=hi[2,t]))
            for k,key in enumerate(axes):
                shifts.append(dict(species=sp,model=model,trait=trait,axis=key,native_beta_per_sd=fit[1+k,t],
                                   predicted_contribution=dx[:,k].mean()*fit[1+k,t]))
        for k,key in enumerate(axes):
            climate.append(dict(species=sp,model=model,axis=key,native_mean=mean[k],native_sd=sd[k],
                                analogous_mean=left[key].mean(),novel_mean=right[key].mean(),change_sd=dx[:,k].mean()))
    from analysis.v5.transfer_support import single_endpoint_cohorts
    expanded=single_endpoint_cohorts(sp,native,intro,out.parent/"temp")
    for file,records in [("native_transfer_test",rows),("native_transfer_components",shifts),("paired_climate_directions",climate),("native_transfer_alignment",joint),("transfer_endpoint_cohorts",expanded),("transfer_support",support),("transfer_robustness",robustness)]:
        path=out/(file+".csv")
        pd.DataFrame(records).to_csv(path,index=False,header=not path.exists(),mode="a")
