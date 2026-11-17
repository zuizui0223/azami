"""X-only climate feature selection; never reads floral traits, prior P values or RF rankings.

Greedy coverage of the standardized covariance of the INPUT X matrix.
The algorithm uses no biological outcome and is not Boruta, significance,
causal variable selection or an optimized response-prediction model.
"""
from __future__ import annotations
import numpy as np

SEED=20261107

def weighted_covariance(x,weights=None):
    x=np.asarray(x,float)
    if x.ndim!=2 or len(x)<12 or x.shape[1]<2:raise ValueError("X_ONLY_SHAPE_INVALID")
    finite=np.isfinite(x).all(axis=1)
    if finite.sum()<12:raise ValueError("X_ONLY_TOO_FEW_COMPLETE_ROWS")
    z=x[finite]
    if weights is None:w=np.ones(len(z))
    else:
        w=np.asarray(weights,float)[finite]
        if len(weights)!=len(x) or not np.isfinite(w).all() or (w<=0).any():
            raise ValueError("X_ONLY_WEIGHT_INVALID")
    w=w/w.sum()
    mean=np.sum(z*w[:,None],axis=0)
    sd=np.sqrt(np.sum((z-mean)**2*w[:,None],axis=0))
    if (sd<1e-12).any() or not np.isfinite(sd).all():
        raise ValueError("X_ONLY_CONSTANT_CLIMATE_COLUMN")
    standardized=(z-mean)/sd
    cov=(standardized*w[:,None]).T@standardized
    cov=(cov+cov.T)/2
    return cov,standardized,w,finite

def retained_covariance_fraction(cov,selected):
    if not selected:return 0.
    block=cov[np.ix_(selected,selected)]
    ridge=np.eye(len(selected))*1e-8
    explained=np.trace(cov[:,selected]@np.linalg.pinv(block+ridge,rcond=1e-10)@cov[selected,:])
    return float(np.clip(explained/np.trace(cov),0.,1.))

def greedy_rank(cov,names,*,max_k=12,max_abs_corr=.85):
    names=tuple(names);p=len(names)
    if cov.shape!=(p,p) or not np.isfinite(cov).all():raise ValueError("X_ONLY_COVAR_INVALID")
    if not 0<=max_abs_corr<=1:raise ValueError("INVALID_CORRELATION_CAP")
    selected=[];progress=[]
    while len(selected)<min(max_k,p):
        eligible=[i for i in range(p) if i not in selected
                  and (not selected or np.max(np.abs(cov[i,selected]))<max_abs_corr)]
        if not eligible:break
        scored=[(retained_covariance_fraction(cov,[*selected,i]),i) for i in eligible]
        best=sorted(scored,key=lambda t:(-t[0],names[t[1]]))[0]
        selected.append(best[1])
        progress.append({"variable":names[best[1]],"n_selected":len(selected),
                         "climate_reconstruction_fraction":float(best[0])})
    return selected,progress

def select_climate_variables(x,names,*,weights=None,target=.90,max_k=12,
                             max_abs_corr=.85,bootstraps=24,seed=SEED):
    names=tuple(names)
    if len(set(names))!=len(names) or len(names)!=np.asarray(x).shape[1]:
        raise ValueError("X_ONLY_NAMES_INVALID")
    cov,z,w,finite=weighted_covariance(x,weights)
    sequence,progress=greedy_rank(cov,names,max_k=max_k,max_abs_corr=max_abs_corr)
    if not sequence:raise ValueError("X_ONLY_NO_VARIABLES_AVAILABLE")
    selected_k=next((r["n_selected"] for r in progress
                     if r["climate_reconstruction_fraction"]>=target),len(progress))
    selected=sequence[:selected_k]
    # Bootstrap is strictly X-only: use row weights as sampling probabilities.
    frequencies=np.zeros(len(names),dtype=int)
    sizes=[]
    rng=np.random.default_rng(seed)
    for b in range(bootstraps):
        indices=rng.choice(len(z),size=min(len(z),900),replace=True,p=w)
        try:
            c,_,_,_=weighted_covariance(z[indices])
            seq,curves=greedy_rank(c,names,max_k=max_k,max_abs_corr=max_abs_corr)
            take=next((r["n_selected"] for r in curves
                       if r["climate_reconstruction_fraction"]>=target),len(curves))
        except ValueError:continue
        frequencies[seq[:take]]+=1;sizes.append(take)
    if bootstraps and not sizes:raise ValueError("X_ONLY_BOOTSTRAP_ALL_FAILED")
    return {"n_complete_X_rows":int(finite.sum()),
            "n_incomplete_X_rows":int((~finite).sum()),
            "n_candidate_predictors":len(names),
            "rule":"greedy_maximum_standardized_X_covariance_reconstruction",
            "target_reconstruction_fraction":float(target),
            "maximum_pairwise_abs_correlation_with_selected":float(max_abs_corr),
            "maximum_candidate_variable_count":int(max_k),
            "threshold_reached":bool(progress[selected_k-1]["climate_reconstruction_fraction"]>=target),
            "selected_variables":[names[i] for i in selected],
            "selected_count":len(selected),
            "captured_climate_variance_fraction":progress[selected_k-1]["climate_reconstruction_fraction"],
            "greedy_X_only_rank_to_capacity":progress,
            "bootstrap_X_only_selection_frequency":{
                name:float(frequencies[i]/len(sizes)) if sizes else None
                for i,name in enumerate(names)},
            "n_valid_X_only_bootstraps":len(sizes),
            "selection_never_uses_phenotype_or_previous_results":True,
            "not_RF_Boruta_VSURF_or_significance_test":True}
