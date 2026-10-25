"""Fixed-covariance Gaussian GLS; paired hue rank-one multivariate hypotheses.

Conditional on covariance and Gaussian residual assumptions. No fitted lambda,
arbitrary tip-noise term, discrete transition process, or historical climate.
"""
import numpy as np
from scipy import linalg,stats

def fit(values,environment,covariance):
    y=np.asarray(values,float);e=np.asarray(environment,float);v=np.asarray(covariance,float)
    n=len(y)
    if y.ndim!=2 or y.shape[1] not in (1,2) or e.shape!=(n,4) or v.shape!=(n,n) or n<8:
        raise ValueError('Invalid dimensions')
    if not all(np.isfinite(x).all() for x in (y,e,v)) or not np.allclose(v,v.T):
        raise ValueError('Invalid finite covariance or observations')
    sd=e.std(0,ddof=1)
    if np.any(sd<=1e-12):raise ValueError('Constant environment')
    x=np.column_stack([np.ones(n),(e-e.mean(0))/sd])
    if np.linalg.matrix_rank(x)!=5:raise ValueError('Singular climate design')
    L=linalg.cholesky(v,lower=True)
    wx=linalg.solve_triangular(L,x,lower=True);wy=linalg.solve_triangular(L,y,lower=True)
    beta=linalg.lstsq(wx,wy)[0];res=wy-wx@beta;df=n-5;q=y.shape[1]
    sigma=res.T@res/df
    if np.linalg.matrix_rank(sigma)!=q:raise ValueError('Singular residual coordinate covariance')
    inverse=linalg.inv(wx.T@wx);critical=stats.t.ppf(.975,df)
    tests=[]
    for j in range(1,5):
        T2=float(beta[j]@linalg.solve(sigma,beta[j],assume_a='pos')/inverse[j,j])
        F=T2*(df-q+1)/(q*df)
        se=np.sqrt(np.diag(sigma)*inverse[j,j])
        tests.append(dict(p_value=float(stats.f.sf(F,q,df-q+1)),F=F,df1=q,df2=df-q+1,
            coefficient=beta[j].tolist(),conditional_se=se.tolist(),
            lower95=(beta[j]-critical*se).tolist(),upper95=(beta[j]+critical*se).tolist()))
    return dict(coefficients=beta.tolist(),tests=tests,residual_df=df,
        environment_means=e.mean(0).tolist(),environment_sample_sd=sd.tolist(),
        VIF=np.diag(linalg.inv(np.corrcoef(e,rowvar=False))).tolist(),
        whitened_design_condition_number=float(np.linalg.cond(wx)),
        residual_covariance=sigma.tolist(),coefficient_units='original trait coordinate per one within-panel environmental sample SD',
        joint_test='rank-one Hotelling-Lawley F for paired coordinates; scalar case equals squared t',
        conditional_assumptions='fixed covariance, independent whitened multivariate Gaussian errors; tree and photo sampling uncertainty not propagated')
