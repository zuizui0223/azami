"""Matched multicoordinate scale estimates and candidate crossed resampling.

No empirical file entrypoint or ecological authorization. All response units,
cohort membership and nuisance columns must be supplied already qualified.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from .dependence_resampling import cohort_partition, crossed_draw
from .joint_partial_pooling import fit_joint_partial_pooling
from .model_design_diagnostics import _project_out

SCALES = ('within','among','among_minus_within')


def fit_matched_module(response,predictors,taxa,nuisance):
    y,x,g,b = np.asarray(response,float),np.asarray(predictors,float),np.asarray(taxa),np.asarray(nuisance,float)
    if (y.ndim!=2 or x.ndim!=2 or g.ndim!=1 or len(y)==0 or y.shape[1]==0 or x.shape[1]==0
            or len(x)!=len(y) or len(g)!=len(y) or pd.isna(g).any()
            or b.ndim not in (2,3) or len(b)!=len(y) or (b.ndim==3 and b.shape[1]!=y.shape[1])
            or not all(np.isfinite(v).all() for v in (y,x,b))):
        raise ValueError('Declare finite, aligned joint module support without row deletion')
    g = g.astype(str)
    if np.any(g==''):
        raise ValueError('Missing taxon')
    labels = tuple(sorted(set(g)))
    groups = [np.flatnonzero(g==label) for label in labels]
    means_x = np.array([x[i].mean(axis=0) for i in groups])
    means_y = np.array([y[i].mean(axis=0) for i in groups])
    within,among,diagnostics = [],[],[]
    for coordinate in range(y.shape[1]):
        nuisance_coordinate = b if b.ndim==2 else b[:,coordinate,:]
        means_b = np.array([nuisance_coordinate[i].mean(axis=0) for i in groups])
        among_nuisance = np.column_stack([np.ones(len(labels)),means_b])
        # One simultaneous projection uses the same nuisance span for X and Y.
        residual,nuisance_rank = _project_out(np.column_stack([means_x,means_y[:,coordinate]]),among_nuisance,np.ones(len(labels)))
        xx,yy = residual[:,:x.shape[1]],residual[:,-1]
        norms = np.linalg.norm(xx,axis=0)
        if (np.any(norms<=1e-12) or np.linalg.matrix_rank(xx/norms)!=x.shape[1]
                or len(labels)-nuisance_rank-x.shape[1]<=0):
            raise ValueError('Among exposure design is not identifiable; no predictor deletion')
        beta = np.linalg.lstsq(xx/norms,yy,rcond=None)[0]/norms
        fitted = fit_joint_partial_pooling(y[:,coordinate],x,g,nuisance_coordinate)
        if fitted.taxa != labels:
            raise ValueError('Within and among taxon support differ')
        within.append(fitted.hypermean)
        among.append(beta)
        diagnostics.append({'coordinate':coordinate,'within_retained_taxa':len(fitted.taxa),
                            'within_tau2':fitted.tau2.tolist(),'within_optimizer_attempts':fitted.optimizer_attempts,
                            'among_nuisance_rank':nuisance_rank,'among_residual_df':len(labels)-nuisance_rank-x.shape[1]})
    within,among = np.asarray(within),np.asarray(among)
    return {'coefficients':np.stack([within,among,among-within]),'taxa':labels,
            'observations':len(y),'diagnostics':diagnostics}


def quadratic_geometry(covariance):
    covariance = np.asarray(covariance,float)
    if (covariance.ndim!=2 or covariance.shape[0]!=covariance.shape[1] or not len(covariance)
            or not np.isfinite(covariance).all() or not np.allclose(covariance,covariance.T,rtol=1e-10,atol=1e-14)):
        raise ValueError('Finite symmetric covariance required')
    variance = np.diag(covariance)
    if np.any(variance<0):
        raise ValueError('Negative variance')
    scale = np.where(variance>0,np.sqrt(variance),1.)
    correlation = covariance/scale[:,None]/scale[None,:]
    values,vectors = np.linalg.eigh(correlation)
    tolerance = 1e-10*max(1.,float(abs(values).max()))
    if np.any(values < -tolerance):
        raise ValueError('Covariance is not positive semidefinite')
    keep = values>tolerance
    if not keep.any():
        raise ValueError('Rank-zero covariance is not a test')
    return scale,values[keep],vectors[:,keep]


def quadratic_values(effects, geometry):
    effects = np.atleast_2d(np.asarray(effects,float))
    scale,values,vectors = geometry
    if effects.shape[1]!=len(scale) or not np.isfinite(effects).all():
        raise ValueError('Malformed effects')
    standardized = effects/scale
    projection = standardized@vectors
    remainder = standardized-projection@vectors.T
    if np.any(np.linalg.norm(remainder,axis=1)>1e-8*(1+np.linalg.norm(standardized,axis=1))):
        raise ValueError('Effect outside covariance range; do not discard its unestimated component')
    return np.sum(projection**2/values,axis=1)


def summarize_draws(point,draws,blocks):
    point,draws = np.asarray(point,float),np.asarray(draws,float)
    if (point.ndim!=3 or point.shape[0]!=3 or draws.ndim!=4 or draws.shape[1:]!=point.shape
            or len(draws)<2 or not np.isfinite(point).all() or not np.isfinite(draws).all()):
        raise ValueError('Every planned draw must be estimable before summarizing inference')
    p = point.shape[-1]
    if sorted(i for block in blocks.values() for i in block)!=list(range(p)) or any(not block for block in blocks.values()):
        raise ValueError('Process blocks must partition every supplied predictor')
    vector = draws.reshape(len(draws),-1)
    covariance = np.atleast_2d(np.cov(vector,rowvar=False,ddof=1))
    low,high = np.quantile(draws,[.025,.975],axis=0)
    tests = []
    for scale,name in enumerate(SCALES):
        for block,indices in blocks.items():
            target = point[scale][:,indices].ravel()
            samples = draws[:,scale][:,:,indices].reshape(len(draws),-1)
            c = np.atleast_2d(np.cov(samples,rowvar=False,ddof=1))
            row = {'scale':name,'process':block,'coefficients_tested':len(target)}
            try:
                geometry = quadratic_geometry(c)
                observed = float(quadratic_values(target,geometry)[0])
                null = quadratic_values(samples-target,geometry)
                row.update(status='candidate_tail_estimated_not_calibrated',rank=len(geometry[1]),
                           statistic=observed,candidate_probability=(1+int(np.sum(null>=observed)))/(len(draws)+1))
            except ValueError as error:
                row.update(status='not_estimable',reason=str(error),candidate_probability=None)
            tests.append(row)
    return {'covariance':covariance,'basic_interval_low':2*point-high,'basic_interval_high':2*point-low,
            'candidate_tests':tests,'ecological_fitting_authorized':False}


def bootstrap_matched_module(response,predictors,nuisance,source,positions,*,seed,replicates,blocks,fit=fit_matched_module):
    """Paired source-factor bootstrap. Failed replicates are never replaced."""
    positions,_,_ = cohort_partition(source,positions)
    y,x,b = np.asarray(response,float),np.asarray(predictors,float),np.asarray(nuisance,float)
    if len(y)!=len(positions) or len(x)!=len(positions) or len(b)!=len(positions):
        raise ValueError('Source positions and module arrays differ')
    if not isinstance(replicates,int) or isinstance(replicates,bool) or replicates<2:
        raise ValueError('At least two planned replicates required')
    point = fit(y,x,source.taxa[positions],b)
    draws = np.full((replicates,*point['coefficients'].shape),np.nan)
    records = []
    for replicate in range(replicates):
        row = {'replicate':replicate,'seed':seed,'grid_degrees':source.grid_degrees}
        try:
            indices,copy_taxa = crossed_draw(source,positions,seed=seed,replicate=replicate)
            fitted = fit(y[indices],x[indices],copy_taxa,b[indices])
            if fitted['coefficients'].shape!=point['coefficients'].shape or not np.isfinite(fitted['coefficients']).all():
                raise ValueError('Invalid complete coefficient vector')
            draws[replicate] = fitted['coefficients']
            row.update(status='estimated',sampled_observations=len(indices),sampled_taxon_copies=len(fitted['taxa']),
                       diagnostics=fitted['diagnostics'])
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
            row.update(status='not_estimable',error_type=type(error).__name__,reason=str(error),
                       optimizer_attempts=getattr(error,'optimizer_attempts',[]))
        records.append(row)
    complete = all(r['status']=='estimated' for r in records)
    summary = summarize_draws(point['coefficients'],draws,blocks) if complete else None
    return {'point':point,'draws':draws,'records':records,'summary':summary,
            'planned_replicates':replicates,'estimable_replicates':sum(r['status']=='estimated' for r in records),
            'status':'candidate_bootstrap_complete_not_calibrated' if complete else 'bootstrap_incomplete_no_inference',
            'ecological_fitting_authorized':False}
