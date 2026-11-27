"""Can the Part 1 water axes reconstruct Part 2's enrichment curves?

A fixed explanatory comparison, not a search over climate variables,
significance or thresholds. Species-specific native ranks are retained.
Condition on the pooled two-cell phenotype distribution; predict its
allocation to analogous/non-analogue cells from actual signed climate
differences. Refit the original six-term curve to predicted allocations.
No causal selection, fitness, inheritance or absolute-site prediction.
"""
from pathlib import Path
import json
from collections import Counter
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from analysis.v5.phenotype_curve_mechanism import (
    ROOT,SPECIES,DEFINITIONS,VARIABLES,RAIN,SEASON,NAMES,EXPECTED_BETA,
    authority_source,angular_cells,background_mesh,apply_original_five_CHELSA,
    load_three_axis,_background_species_audit,_pca_env_2d,run_species,
    X,K,fit,marginals,
)
OUT=ROOT/"phenotype_water_curve_reconstruction"
RIDGE=4.
CHROMA=NAMES.index("chroma")
ANGLE=NAMES.index("angle")
DIRECTION=(CHROMA,ANGLE)
CURVATURE=(CHROMA+3,ANGLE+3)
MODELS=("null","constant_B","shared_shift","species_shift","shared_shape","species_shape")
WATER_MODELS=MODELS[2:]
TERMS=("BIO15_x_chroma_rank","BIO12_x_angle_rank",
       "BIO15_x_chroma_rank_squared","BIO12_x_angle_rank_squared")

def inputs(d):
    # Both columns are signed NOVEL minus ANALOGOUS. The physically
    # comparable units are fixed to the earlier water analysis.
    x=np.column_stack(((d["right"][:,SEASON]-d["left"][:,SEASON])/20.,
                       (d["right"][:,RAIN]-d["left"][:,RAIN])/500.))
    a,b=d["hist_left"],d["hist_right"]
    np.testing.assert_allclose(a.sum(1),1.,atol=1e-10)
    np.testing.assert_allclose(b.sum(1),1.,atol=1e-10)
    return dict(a=a,b=b,x=x,g=d["blocks"])

def subset(d,mask):
    return {k:v[mask] for k,v in d.items()}

def valid(d):
    return len(d["a"])>=20 and len(np.unique(d["g"]))>=7

def design(x,shape):
    # Climate-level differences times rank features. Under a common
    # linear climate-conditioned exponential family, log p_N/p_A
    # has exactly this form. Tilt-only has no local quadratic term.
    x=np.asarray(x,float)
    if x.ndim!=2 or x.shape[1]!=2 or not np.isfinite(x).all():
        raise ValueError("invalid water differences")
    z=x[:,:,None]*X[:,DIRECTION].T[None,:,:]
    if shape:
        z=np.concatenate((z,x[:,:,None]*X[:,CURVATURE].T[None,:,:]),axis=1)
    return np.swapaxes(z,1,2)

def profile(a,b,eta):
    m=a+b
    if eta.shape!=a.shape or not np.isfinite(eta).all():
        raise ValueError("invalid allocation predictor")
    lower=-eta.max(1)-40.;upper=-eta.min(1)+40.
    for _ in range(48):
        mid=(lower+upper)/2
        p=expit(eta+mid[:,None])
        above=(m*p).sum(1)>1.
        upper=np.where(above,mid,upper)
        lower=np.where(above,lower,mid)
    v=eta+((lower+upper)/2)[:,None]
    p=expit(v)
    np.testing.assert_allclose((m*p).sum(1),1.,atol=1e-10,rtol=0)
    ls=(m*np.logaddexp(0.,v)-b*v).sum(1)
    return ls,p

def objective(beta,a,b,z,w,penalty):
    eta=z@beta
    ls,p=profile(a,b,eta)
    r=(a+b)*p-b
    grad=np.einsum("n,nk,nkj->j",w,r,z)+penalty*beta
    return float(w@ls+.5*penalty*(beta@beta)),grad

def fit_water(a,b,x,w,shape,penalty):
    z=design(x,shape)
    p=z.shape[-1]
    if not np.isfinite(w).all() or np.min(w)<=0:
        raise ValueError("invalid weights")
    # Identification diagnostic: after profiling pair intercepts, every
    # active column must have nonzero independent within-pair variation.
    m=a+b
    means=np.einsum("nk,nkj->nj",m,z)/2.
    centered=(z-means[:,None,:])*np.sqrt(m[:,:,None]*w[:,None,None])
    if np.linalg.matrix_rank(centered.reshape(-1,p))<p:
        raise ValueError("rank_deficient")
    ans=minimize(objective,np.zeros(p),args=(a,b,z,w,penalty),jac=True,
                 method="L-BFGS-B",options=dict(maxiter=300,ftol=1e-14,gtol=1e-9,maxls=40))
    val,grad=objective(ans.x,a,b,z,w,penalty)
    if (not np.isfinite(val) or not np.isfinite(ans.x).all() or
            np.max(np.abs(grad))>1e-5 or np.max(np.abs(ans.x))>12):
        raise ValueError("numerical_fit_failed")
    return ans.x

def train_models(train):
    # A shared fit uses the same objective and penalty as both species
    # fits under the restriction that their slopes are identical.
    # Scale species total training weight equally (N/2) so the larger
    # species cannot determine the common response alone.
    n=np.array([len(d["a"]) for d in train])
    w=[np.full(n[i],n.sum()/(2*n[i])) if n[i] else np.empty(0) for i in range(2)]
    out={name:[None,None] for name in MODELS}
    out["null"]=[np.zeros(6),np.zeros(6)]
    fail={name:["",""] for name in MODELS}
    for i,d in enumerate(train):
        if not valid(d):
            for model in MODELS[1:]:
                fail[model][i]="training_support"
            continue
        try:
            out["constant_B"][i]=fit(d["a"],d["b"],X,np.arange(6))[0]
        except (ValueError,np.linalg.LinAlgError):
            fail["constant_B"][i]="numerical_fit_failed"
        for model in ("species_shift","species_shape"):
            try:
                out[model][i]=fit_water(d["a"],d["b"],d["x"],w[i],
                                            model.endswith("shape"),RIDGE)
            except (ValueError,np.linalg.LinAlgError) as exc:
                fail[model][i]="rank_deficient" if str(exc)=="rank_deficient" else "numerical_fit_failed"
    for model in ("shared_shift","shared_shape"):
        if not all(valid(d) for d in train):
            fail[model]=["training_support","training_support"]
            continue
        try:
            a,b,x=(np.concatenate([d[k] for d in train],axis=0) for k in ("a","b","x"))
            beta=fit_water(a,b,x,np.concatenate(w),model.endswith("shape"),2*RIDGE)
            out[model]=[beta,beta]
        except (ValueError,np.linalg.LinAlgError) as exc:
            reason="rank_deficient" if str(exc)=="rank_deficient" else "numerical_fit_failed"
            fail[model]=[reason,reason]
    return out,fail

def allocate(d,model,beta):
    a,b=d["a"],d["b"]
    if model=="null":
        eta=np.zeros_like(a)
    elif model=="constant_B":
        eta=np.broadcast_to(X@beta,a.shape)
    else:
        eta=design(d["x"],model.endswith("shape"))@beta
    ls,p=profile(a,b,eta)
    m=a+b
    pa,pb=m*(1-p),m*p
    np.testing.assert_allclose(pa.sum(1),1.,atol=1e-10)
    np.testing.assert_allclose(pb.sum(1),1.,atol=1e-10)
    np.testing.assert_allclose(pa+pb,m,atol=1e-12,rtol=0)
    return ls,pa,pb

def reconstruct(a,b,pa,pb,observed):
    # Re-fit B on the predicted allocations, not on observed novel labels.
    # The observed B provides the comparison target only.
    bp=fit(pa,pb,X,np.arange(6))[0]
    if not np.isfinite(bp).all():
        raise ValueError("nonfinite reconstructed B")
    rows=[]
    for j,trait in enumerate(NAMES):
        for rank in range(-2,3):
            obs=np.exp(observed[j]*rank+observed[j+3]*rank**2)
            pred=np.exp(bp[j]*rank+bp[j+3]*rank**2)
            rows.append(dict(trait=trait,native_quintile=rank+3,
                    observed_linear=float(observed[j]),
                    observed_quadratic=float(observed[j+3]),
                    reconstructed_linear=float(bp[j]),
                    reconstructed_quadratic=float(bp[j+3]),
                    observed_relative_enrichment=float(obs),
                    reconstructed_relative_enrichment=float(pred)))
    return rows,bp

def summary(d,pred,sp,stage,meta,observed):
    skill=[];curves=[];profiles=[];coverage=[]
    for model in MODELS:
        ls,pa,pb,done,reasons=pred[model]
        complete=bool(np.all(done))
        row=dict(**meta,species=sp,stage=stage,model=model,
                 status="estimated" if complete else "incomplete",
                 valid_pairs=int(done.sum()),valid_regions=int(sum(np.all(done[d["g"]==g]) for g in np.unique(d["g"]))))
        coverage.append(dict(**row,**{k:int(reasons.get(k,0)) for k in ("training_support","rank_deficient","numerical_fit_failed","allocation_failed")}))
        if complete:
            byregion=lambda v:np.array([v[d["g"]==g].mean() for g in np.unique(d["g"])])
            row.update(mean_loss=float(ls.mean()),equal_region_loss=float(byregion(ls).mean()))
            for reference in ("null","constant_B","shared_shift","species_shift","shared_shape"):
                rl,_,_,rdone,_=pred[reference]
                if np.all(rdone):
                    row["skill_vs_"+reference+"_pct"]=float(100*(rl.mean()-ls.mean())/rl.mean())
                    row["equal_region_gain_vs_"+reference]=float((byregion(rl)-byregion(ls)).mean())
                    row["regions_better_than_"+reference]=int(((byregion(rl)-byregion(ls))>1e-12).sum())
            try:
                rr,beta=reconstruct(d["a"],d["b"],pa,pb,observed)
                for r in rr:
                    curves.append(dict(**meta,species=sp,stage=stage,model=model,status="estimated",**r))
                for j,trait in enumerate(NAMES):
                    row[trait+"_logcurve_rmse"]=float(np.sqrt(np.mean([
                        (np.log(r["observed_relative_enrichment"])-np.log(r["reconstructed_relative_enrichment"]))**2
                        for r in rr if r["trait"]==trait])))
            except (ValueError,np.linalg.LinAlgError):
                curves.append(dict(**meta,species=sp,stage=stage,model=model,status="reconstruction_failed"))
            aa,bb,ap,bp=(marginals(h).mean(0) for h in (d["a"],d["b"],pa,pb))
            for j,trait in enumerate(NAMES):
                for k in range(5):
                    profiles.append(dict(**meta,species=sp,stage=stage,model=model,trait=trait,native_quintile=k+1,
                        observed_analogous=float(aa[j,k]),observed_novel=float(bb[j,k]),
                        predicted_analogous=float(ap[j,k]),predicted_novel=float(bp[j,k]),
                        observed_difference=float(bb[j,k]-aa[j,k]),
                        predicted_difference=float(bp[j,k]-ap[j,k])))
        skill.append(row)
    return skill,curves,profiles,coverage

def empty(d):
    return {m:(np.full(len(d["a"]),np.nan),np.full_like(d["a"],np.nan),
               np.full_like(d["b"],np.nan),np.zeros(len(d["a"]),bool),Counter()) for m in MODELS}

def evaluate(d,indices,params,reasons,pred,i):
    if not np.any(indices):
        return
    test=subset(d,indices)
    for model in MODELS:
        ls,pa,pb,done,fail=pred[model]
        beta=params[model][i]
        if beta is None:
            fail[reasons[model][i] or "numerical_fit_failed"]+=1
            continue
        try:
            l,a,b=allocate(test,model,beta)
        except (ValueError,np.linalg.LinAlgError,AssertionError):
            fail["allocation_failed"]+=1
            continue
        ls[indices]=l;pa[indices]=a;pb[indices]=b;done[indices]=True

def run_definition(raw):
    data=[inputs(d) for d in raw]
    observed=[fit(d["a"],d["b"],X,np.arange(6))[0] for d in data]
    full_params,full_fail=train_models(data)
    parameters=[]
    for model in WATER_MODELS:
        for i,sp in enumerate(SPECIES):
            beta=full_params[model][i]
            for k,term in enumerate(TERMS[:4 if model.endswith("shape") else 2]):
                parameters.append(dict(species=sp,model=model,term=term,
                     status="estimated" if beta is not None else full_fail[model][i],
                     coefficient=float(beta[k]) if beta is not None else np.nan))
    ins=[empty(d) for d in data];held=[empty(d) for d in data]
    for i,d in enumerate(data):
        evaluate(d,np.ones(len(d["a"]),bool),full_params,full_fail,ins[i],i)
    # Hold out the same geographic region in both species, preventing
    # a shared model from learning a held-out region via the other species.
    groups=np.unique(np.concatenate([d["g"] for d in data]))
    for g in groups:
        train=[subset(d,d["g"]!=g) for d in data]
        params,reasons=train_models(train)
        for i,d in enumerate(data):
            evaluate(d,d["g"]==g,params,reasons,held[i],i)
    skill=[];curves=[];profiles=[];coverage=[]
    for i,sp in enumerate(SPECIES):
        d=data[i];meta=dict(pairs=len(d["a"]),regions=len(np.unique(d["g"])))
        for stage,p in (("in_sample",ins[i]),("heldout",held[i])):
            s,c,f,v=summary(d,p,sp,stage,meta,observed[i])
            skill.extend(s);curves.extend(c);profiles.extend(f);coverage.extend(v)
    return skill,curves,profiles,parameters,coverage

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=load_three_axis(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    skill=[];curves=[];profiles=[];pars=[];coverage=[]
    for definition in DEFINITIONS:
        data=[out[(sp,definition)] for sp in SPECIES]
        res=run_definition(data)
        for name,rows in zip(("skill","curves","profiles","parameters","coverage"),res):
            target={"skill":skill,"curves":curves,"profiles":profiles,"parameters":pars,"coverage":coverage}[name]
            target.extend(dict(definition=definition,**row) for row in rows)
    pd.DataFrame(skill).to_csv(OUT/"water_reconstruction_skill.csv",index=False)
    pd.DataFrame(curves).to_csv(OUT/"water_reconstructed_B.csv",index=False)
    pd.DataFrame(profiles).to_csv(OUT/"water_predicted_profiles.csv",index=False)
    pd.DataFrame(pars).to_csv(OUT/"water_response_parameters.csv",index=False)
    pd.DataFrame(coverage).to_csv(OUT/"water_reconstruction_coverage.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
