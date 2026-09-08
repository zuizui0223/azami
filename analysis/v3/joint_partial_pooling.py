"""Joint Gaussian random-slope numerics retaining locally aliased taxa.

Working homoskedastic likelihood only. These utilities neither choose a cohort
nor authorize ecological fitting. Spatial/component and multivariate inference
must use the separately validated joint resampling runner, not conditional SEs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import minimize, root


@dataclass
class PoolingFit:
    taxa: tuple[str, ...]
    hypermean: np.ndarray
    nuisance_coefficients: np.ndarray
    fixed_covariance: np.ndarray
    taxon_predictions: np.ndarray
    prediction_error_covariance: np.ndarray
    variance_ratios: np.ndarray
    tau2: np.ndarray
    residual_variance: float
    residual_df: int
    reml_criterion: float
    residuals: np.ndarray
    support: pd.DataFrame
    optimizer_attempts: list[dict]
    uncertainty_scope: str = 'conditional_working_model_not_spatially_calibrated'


class PoolingNotEstimable(RuntimeError):
    def __init__(self, message, attempts):
        super().__init__(message)
        self.optimizer_attempts = attempts


def _center(array):
    delta = array - array[0]
    return delta - delta.mean(axis=0)


class JointSlopeLikelihood:
    """Small taxon cross-products, never an N by all-taxon-slope matrix.

    Common-scale predictors are supplied unchanged. Redundant shared nuisance
    columns are represented by an orthogonal basis; exposures are never removed.
    Local rank deficiency is allowed. Global fixed-effect aliasing is rejected.
    """
    def __init__(self, response, predictors, taxa, nuisance):
        y, x, labels, b = np.asarray(response,float), np.asarray(predictors,float), np.asarray(taxa), np.asarray(nuisance,float)
        if (y.ndim!=1 or x.ndim!=2 or x.shape[1]==0 or b.ndim!=2 or labels.ndim!=1
                or len(y)==0 or len(x)!=len(y) or len(b)!=len(y) or len(labels)!=len(y)):
            raise ValueError('Joint hierarchy input dimensions differ')
        if not np.isfinite(y).all() or not np.isfinite(x).all() or not np.isfinite(b).all() or pd.isna(labels).any():
            raise ValueError('Declare finite aligned support; no silent row deletion')
        labels = labels.astype(str)
        if (labels=='').any():
            raise ValueError('Missing taxon label')
        self.taxa = tuple(sorted(set(labels)))
        self.p = x.shape[1]
        self.indices = [np.flatnonzero(labels==label) for label in self.taxa]
        yc, xc, bc = np.empty_like(y), np.empty_like(x), np.empty_like(b)
        for index in self.indices:
            yc[index], xc[index], bc[index] = _center(y[index]), _center(x[index]), _center(b[index])
        self.n_within = len(y)-len(self.taxa)
        norms = np.linalg.norm(bc,axis=0)
        keep = norms > 1e-12
        scaled = bc[:,keep]/norms[keep]
        u, singular, vt = np.linalg.svd(scaled,full_matrices=False)
        tol = (singular[0] if len(singular) else 0)*max(scaled.shape)*np.finfo(float).eps
        active = singular > tol
        self.nuisance_rank = int(active.sum())
        self.nuisance_map = np.zeros((b.shape[1],self.nuisance_rank))
        self.nuisance_map[keep] = vt[active].T / norms[keep,None] / singular[active] * np.sqrt(max(1,self.n_within))
        nuisance_basis = bc @ self.nuisance_map
        fixed = np.column_stack([xc,nuisance_basis])
        column_norms = np.linalg.norm(fixed,axis=0)
        if np.any(column_norms<=1e-12) or np.linalg.matrix_rank(fixed/column_norms) < fixed.shape[1]:
            raise ValueError('Common environmental slopes globally aliased with nuisance or other exposures')
        self.fixed_columns = fixed.shape[1]
        self.df = self.n_within-self.fixed_columns
        if self.df <= 0 or len(self.taxa)<2:
            raise ValueError('Insufficient joint residual degrees of freedom or taxa')
        self.x, self.y, self.fixed = xc, yc, fixed
        self.column_norms = column_norms
        self.gram = np.stack([xc[i].T@xc[i] for i in self.indices])
        self.xf = np.stack([xc[i].T@fixed[i] for i in self.indices])
        self.xy = np.stack([xc[i].T@yc[i] for i in self.indices])
        self.ff = sum(fixed[i].T@fixed[i] for i in self.indices)
        self.fy = sum(fixed[i].T@yc[i] for i in self.indices)
        self.yy = float(yc@yc)
        self.support = pd.DataFrame([{
            'taxon':label,'n_observations':len(i),'within_predictor_rank':int(np.linalg.matrix_rank(xc[i])),
            'predictor_count':self.p,'local_full_rank':bool(np.linalg.matrix_rank(xc[i])==self.p),
            'no_local_slope_information':bool(np.linalg.norm(xc[i])==0),
            'retained_in_joint_model':True,
        } for label,i in zip(self.taxa,self.indices)])

    def evaluate(self, variance_ratios):
        lam = np.asarray(variance_ratios,float)
        if lam.shape!=(self.p,) or not np.isfinite(lam).all() or np.any(lam<0):
            raise ValueError('Variance ratios must be finite and nonnegative')
        root = np.sqrt(lam)
        a = np.eye(self.p) + root[None,:,None]*self.gram*root[None,None,:]
        chol = np.linalg.cholesky(a)
        m = root[None,:,None]*np.linalg.solve(a,np.broadcast_to(np.diag(root),a.shape))
        info = self.ff - np.einsum('tpi,tpq,tqj->ij',self.xf,m,self.xf)
        info = (info+info.T)/2
        score = self.fy - np.einsum('tpi,tpq,tq->i',self.xf,m,self.xy)
        scaled = info / self.column_norms[:,None] / self.column_norms[None,:]
        info_chol = np.linalg.cholesky(scaled)
        inverse = np.linalg.solve(scaled,np.eye(self.fixed_columns)) / self.column_norms[:,None] / self.column_norms[None,:]
        inverse = (inverse+inverse.T)/2
        alpha = inverse@score
        q = self.yy - np.einsum('tp,tpq,tq->',self.xy,m,self.xy)
        rss = float(q-score@alpha)
        if rss <= max(1e-12, 1e-12*self.yy) or not np.isfinite(rss):
            raise ValueError('Residual scale is not numerically estimable')
        logdet_info = 2*np.log(np.diag(info_chol)).sum()
        logdet_info += 2*np.log(self.column_norms).sum()
        criterion = float(self.df*np.log(rss/self.df) + 2*np.log(np.diagonal(chol,axis1=1,axis2=2)).sum() + logdet_info)
        residual_score = self.xy-np.einsum('tpi,i->tp',self.xf,alpha)
        random_modes = np.einsum('tpq,tq->tp',m,residual_score)
        vxf = self.xf-self.gram@m@self.xf
        vxresid = residual_score-np.einsum('tpq,tq->tp',self.gram,random_modes)
        vxgram = self.gram-self.gram@m@self.gram
        trace = np.diagonal(vxgram,axis1=1,axis2=2).sum(axis=0)
        trace -= np.einsum('tpi,ij,tpj->p',vxf,inverse,vxf)
        gradient = trace-self.df/rss*(vxresid**2).sum(axis=0)
        return {'criterion':criterion,'gradient':gradient,'lambda':lam,'alpha':alpha,
                'fixed_inverse_information':inverse,'rss':rss,'sigma2':rss/self.df,
                'conditional_random_covariance_ratio':m,'random_modes':random_modes,'vxf':vxf}

    def result(self, evaluation, attempts=None):
        p, alpha, sigma2 = self.p, evaluation['alpha'], evaluation['sigma2']
        fixed_cov = sigma2*evaluation['fixed_inverse_information']
        selector = np.column_stack([np.eye(p),np.zeros((p,self.nuisance_rank))])
        mapping = selector[None,:,:] - evaluation['lambda'][None,:,None]*evaluation['vxf']
        stacked = mapping.reshape(-1,self.fixed_columns)
        prediction_cov = stacked@fixed_cov@stacked.T
        for t in range(len(self.taxa)):
            sl = slice(t*p,(t+1)*p)
            prediction_cov[sl,sl] += sigma2*evaluation['conditional_random_covariance_ratio'][t]
        prediction_cov = (prediction_cov+prediction_cov.T)/2
        residuals = self.y-self.fixed@alpha
        for t,index in enumerate(self.indices):
            residuals[index] -= self.x[index]@evaluation['random_modes'][t]
        return PoolingFit(self.taxa,alpha[:p],self.nuisance_map@alpha[p:],fixed_cov,
                          alpha[:p][None,:]+evaluation['random_modes'],prediction_cov,
                          evaluation['lambda'],sigma2*evaluation['lambda'],sigma2,self.df,
                          evaluation['criterion'],residuals,self.support.copy(),attempts or [])


def _projected_gradient(position, gradient):
    projected = gradient.copy()
    projected[(position<=1e-10)&(gradient>0)] = 0
    return projected


def _polish_stationarity(objective, position, ceiling):
    """Solve the score near a reported minimum, without loosening acceptance.

    Large-N objective subtraction can trigger relative-function convergence
    before the analytic score is small. Keep the same likelihood and bounds;
    only accept a stationary, positive-curvature, non-worse nearby solution.
    """
    before, gradient = objective(position)
    free = np.flatnonzero(position>1e-10)
    record = {'attempted':True,'accepted':False,'initial_projected_gradient_max':float(abs(_projected_gradient(position,gradient)).max())}
    if not len(free):
        return position, record
    def expand(values):
        candidate = position.copy(); candidate[free] = values
        if np.any(candidate<0) or np.any(candidate>=ceiling):
            raise ValueError('Score polishing left the original variance bounds')
        return candidate
    try:
        fitted = root(lambda values:objective(expand(values))[1][free],position[free],method='hybr',options={'xtol':1e-9})
        candidate = expand(fitted.x)
        after, score = objective(candidate)
        maximum = float(abs(_projected_gradient(candidate,score)).max())
        hessian = np.empty((len(free),len(free)))
        for column,index in enumerate(free):
            step = min(1e-5,candidate[index]/2,(ceiling-candidate[index])/2)
            if step<=1e-12:
                raise ValueError('Interior curvature cannot be checked at this boundary')
            delta = np.eye(len(position))[index]*step
            hessian[:,column] = (objective(candidate+delta)[1][free]-objective(candidate-delta)[1][free])/(2*step)
        curvature = float(np.linalg.eigvalsh((hessian+hessian.T)/2).min())
        accepted = bool(fitted.success and maximum<=1e-4 and curvature>0
                        and after<=before+max(1e-7,1e-10*abs(before)))
        record.update(success=bool(fitted.success),accepted=accepted,projected_gradient_max=maximum,
                      minimum_free_curvature=curvature,criterion_change=float(after-before),
                      maximum_log1p_parameter_change=float(abs(candidate-position).max()))
        return (candidate if accepted else position),record
    except (ValueError,np.linalg.LinAlgError) as error:
        record['error_type'] = type(error).__name__
        return position,record


def fit_joint_partial_pooling(response, predictors, taxa, nuisance):
    likelihood = JointSlopeLikelihood(response,predictors,taxa,nuisance)
    attempts, evaluations = [], []
    ceiling = float(np.log1p(1e6))
    def objective(value):
        out = likelihood.evaluate(np.expm1(value))
        return out['criterion'], out['gradient']*np.exp(value)
    for initial in (0., .1, 1., 10.):
        result = minimize(objective,np.full(likelihood.p,np.log1p(initial)),jac=True,method='L-BFGS-B',
                          bounds=[(0.,ceiling)]*likelihood.p,
                          options={'maxiter':500,'ftol':1e-12,'gtol':1e-6,'maxls':40})
        position = result.x
        polishing = {'attempted':False,'accepted':False}
        if result.success and np.max(abs(_projected_gradient(position,objective(position)[1])))>1e-4:
            position,polishing = _polish_stationarity(objective,position,ceiling)
        evaluation = likelihood.evaluate(np.expm1(position))
        gradient = evaluation['gradient']*np.exp(position)
        projected = _projected_gradient(position,gradient)
        stationary = float(np.max(np.abs(projected))) <= 1e-4
        upper = bool(np.any(position>=ceiling-1e-6))
        accepted = bool(result.success and stationary and not upper)
        attempts.append({'initial_variance_ratio':initial,'success':bool(result.success),'accepted':accepted,
                         'stationary':stationary,'upper_bound_hit':upper,'criterion':evaluation['criterion'],
                         'projected_gradient_max':float(np.max(np.abs(projected))),
                         'iterations':int(result.nit),'variance_ratios':evaluation['lambda'].tolist(),
                         'stationarity_polishing':polishing})
        evaluations.append(evaluation)
    accepted = [i for i,a in enumerate(attempts) if a['accepted']]
    if not accepted:
        raise PoolingNotEstimable('No converged interior-or-zero-boundary REML solution; do not report estimates',attempts)
    best = min(accepted,key=lambda i:attempts[i]['criterion'])
    if min(a['criterion'] for a in attempts) < attempts[best]['criterion']-1e-6:
        raise PoolingNotEstimable('Unresolved lower REML objective from failed attempt',attempts)
    return likelihood.result(evaluations[best],attempts)
