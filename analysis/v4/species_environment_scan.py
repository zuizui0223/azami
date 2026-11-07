"""Exploratory species-level atlas, not RF, confirmation, or source recovery.

Reads only already-authorized, physical-unit observations. No networking,
credentials, decryption, workflows, or individual-record export. Marginal HC3
results are diagnostic; supplied joint-covariance GLS is separately labelled.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy import linalg, stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests

FEATURES = {
    'angle': ('orientation_image_vertical_angle',),
    'lightness': ('corolla_lab_lightness',),
    'chroma': ('corolla_lab_chroma',),
    'hue': ('corolla_hue_sin', 'corolla_hue_cos'),
    'aspect': ('capitulum_outline_aspect_ratio',),
    'circularity': ('capitulum_outline_circularity',),
    'solidity': ('capitulum_outline_solidity',),
    'width_profile': ('capitulum_width_profile_cv',),
}
PREDICTORS = tuple(f'chelsa_bio{i:02}' for i in range(1,20)) + (
    'chelsa_rsds_mean','topo_elevation','chelsa_vpd_mean','chelsa_sfcwind_mean')
FAMILIES = {p: ('temperature' if int(p[-2:])<=11 else 'precipitation')
            for p in PREDICTORS[:19]}
FAMILIES.update(chelsa_rsds_mean='radiation',topo_elevation='elevation',
                chelsa_vpd_mean='atmospheric_dryness',chelsa_sfcwind_mean='wind')
PANELS=('marginal','marginal_covariance_matched','joint_covariance_GLS')
MIN_OBSERVATIONS,MIN_TAXA,EXPECTED_N=5,20,313350


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()


def validate_observations(d:pd.DataFrame)->None:
    if (not {'obs_id','accepted_key'} <= set(d) or d.obs_id.isna().any()
        or not d.obs_id.is_unique or d.accepted_key.isna().any()
        or d.obs_id.astype(str).str.strip().eq('').any()
        or d.accepted_key.astype(str).str.strip().eq('').any()):
        raise ValueError('SOURCE_IDENTITY_INVALID')


def aggregate_pair(d:pd.DataFrame,construct:str,predictor:str)->pd.DataFrame:
    """PRIVATE return: one taxon row; X/Y means use identical photo records.

    No minimum geographical-cell count, range extent, or 40-photo requirement.
    The five-observation floor is fixed, not chosen by the effect strength.
    """
    if construct not in FEATURES or predictor not in PREDICTORS:
        raise ValueError('UNREGISTERED_CONSTRUCT_OR_PREDICTOR')
    cols=FEATURES[construct]
    needed={'accepted_key',predictor,*cols,*(c+'__eligible' for c in cols)}
    output_cols=['taxon','x',*[f'y{i}' for i in range(len(cols))],'n_observations']
    if not needed <= set(d):
        out=pd.DataFrame(columns=output_cols)
        out.attrs['empty_reason']='source_columns_unavailable';return out
    v=d[[predictor,*cols]].apply(pd.to_numeric,errors='coerce')
    ok=np.isfinite(v.to_numpy(float)).all(axis=1)
    for c in cols:ok &= d[c+'__eligible'].eq(True).fillna(False).to_numpy(bool)
    if construct=='hue':
        norm=np.linalg.norm(v[list(cols)].to_numpy(float),axis=1)
        ok &= norm>=.20
        v.loc[ok,list(cols)]=v.loc[ok,list(cols)].to_numpy()/norm[ok,None]
    w=v.loc[ok].copy();w['taxon']=d.loc[ok,'accepted_key'].astype(str).to_numpy()
    g=w.groupby('taxon',sort=True,observed=True)
    out=g[[predictor,*cols]].mean().rename(columns={predictor:'x',**{c:f'y{i}' for i,c in enumerate(cols)}})
    out['n_observations']=g.size();out=out.loc[out.n_observations>=MIN_OBSERVATIONS].reset_index()
    if construct=='hue' and len(out):
        norm=np.linalg.norm(out[['y0','y1']].to_numpy(float),axis=1)
        out=out.loc[norm>=.20].copy()
        out[['y0','y1']]=out[['y0','y1']].to_numpy()/norm[norm>=.20,None]
    out=out[output_cols].reset_index(drop=True)
    out.attrs['empty_reason']='insufficient_paired_observations'
    return out


@dataclass
class JointCovariance:
    taxa:tuple[str,...]
    matrix:np.ndarray
    parameter_scope:str
    source_hashes:dict

    def subset(self,taxa:list[str])->np.ndarray:
        lookup={t:i for i,t in enumerate(self.taxa)}
        if not set(taxa)<=set(lookup):raise ValueError('COVARIANCE_TAXA_MISSING')
        ix=[lookup[t] for t in taxa]
        return self.matrix[np.ix_(ix,ix)]


def covariance_from_payload(d:dict)->JointCovariance:
    """Verify a supplied TOTAL residual covariance, not just an arbitrary tree.

    Its construction/fitting must be documented upstream. These checks do not
    prove its biological validity or propagate covariance-estimation uncertainty.
    """
    try:
        taxa=tuple(d['taxa']);V=np.asarray(d['matrix'],float)
        hs=d['source_hashes'];scope=d['parameter_scope']
        valid=(d['kind']=='phylogeny_and_space_total_residual'
            and len(taxa)>=2 and len(set(taxa))==len(taxa)
            and all(isinstance(x,str) and x.strip() for x in taxa)
            and V.shape==(len(taxa),len(taxa)) and np.isfinite(V).all()
            and np.allclose(V,V.T,atol=1e-10,rtol=1e-10)
            and scope in ('externally_fixed','null_model_estimated')
            and all(re.fullmatch('[a-f0-9]{64}',str(hs.get(k,''))) for k in ('phylogeny','space')))
        if not valid:raise ValueError()
        np.linalg.cholesky(V)
    except (KeyError,TypeError,ValueError,np.linalg.LinAlgError):
        raise ValueError('COVARIANCE_CONTRACT_INVALID') from None
    return JointCovariance(taxa,V,scope,hs)


def fit_pair(p:pd.DataFrame,construct:str,V:np.ndarray|None=None)->dict:
    n=len(p);q=len(FEATURES[construct])
    support=None if not n else [int(p.n_observations.min()),
        float(p.n_observations.median()),int(p.n_observations.max())]
    base={'n_taxa':int(n),'p_value':None,'joint_dimensions':q,'beta_raw':None,
          'n_photo_contributions':int(p.n_observations.sum()),
          'observations_per_taxon_min_median_max':support}
    if n<MIN_TAXA:return {**base,'status':p.attrs.get('empty_reason','insufficient_taxa') if n==0 else 'insufficient_taxa'}
    x=p.x.to_numpy(float);Y=p[[f'y{i}' for i in range(q)]].to_numpy(float)
    sx=float(x.std(ddof=1));sy=float(np.sqrt(np.var(Y,axis=0,ddof=1).sum()))
    if sx<1e-12:return {**base,'status':'constant_predictor'}
    if sy<1e-12:return {**base,'status':'constant_response'}
    X=np.column_stack((np.ones(n),(x-x.mean())/sx));df=n-2
    if V is None:
        A=X;B=Y
    else:
        V=np.asarray(V,float)
        if V.shape!=(n,n) or not np.isfinite(V).all() or not np.allclose(V,V.T):
            raise ValueError('INVALID_FIT_COVARIANCE')
        L=np.linalg.cholesky(V)
        A=linalg.solve_triangular(L,X,lower=True)
        B=linalg.solve_triangular(L,Y,lower=True)
    inv=np.linalg.inv(A.T@A);coef=inv@A.T@B;res=B-A@coef;b=coef[1]
    if V is None:
        h=np.einsum('ij,ij->i',A@inv,A)
        if np.max(h)>=1-1e-10:return {**base,'status':'singular_leverage'}
        influence=(A@inv)[:,1,None]*res/(1-h[:,None])
        cb=influence.T@influence
        method='HC3_t_scalar_or_asymptotic_joint_Wald'
    else:
        cb=(res.T@res/df)*inv[1,1]
        method='GLS_t_or_multivariate_F_conditional_on_supplied_total_covariance'
    if np.linalg.eigvalsh(cb).min()<=1e-15:
        return {**base,'status':'degenerate_residual_covariance'}
    if q==1:
        se=float(np.sqrt(cb[0,0]));pval=float(2*stats.t.sf(abs(b[0]/se),df))
        ci=(b[0]+np.array([-1,1])*stats.t.ppf(.975,df)*se)/sx
        details={'beta_raw':float(b[0]/sx),'beta_raw_ci95':ci.tolist(),
                 'beta_standardized':float(b[0]/sy),'slope_se_raw':se/sx}
    else:
        w=float(b@np.linalg.solve(cb,b))
        if V is None:pval=float(stats.chi2.sf(w,q))
        else:pval=float(stats.f.sf(w*(df-q+1)/(q*df),q,df-q+1))
        details={'beta_raw':None,'beta_raw_ci95':None,
            'joint_slope_raw':(b/sx).tolist(),
            'joint_slope_per_environment_sd':b.tolist(),
            'joint_standardized_magnitude':float(np.linalg.norm(b)/sy),
            'hue_has_no_single_signed_linear_coefficient':True}
    return {**base,**details,'status':'estimated','p_value':pval,
            'method':method,'environment_sd_among_taxa':sx,
            'taxon_equal_weighted':True,'causal_or_confirmed_discovery':False}


def adjust_family(rows:list[dict])->None:
    vals=np.array([r.get('p_value') if r.get('p_value') is not None else 1. for r in rows],float)
    if not len(vals):return
    if not np.isfinite(vals).all() or ((vals<0)|(vals>1)).any():raise ValueError('INVALID_P_VALUES')
    bh=multipletests(vals,method='fdr_bh')[1];by=multipletests(vals,method='fdr_by')[1]
    for r,a,b in zip(rows,bh,by):
        r.update(q_bh=None if r.get('p_value') is None else float(a),
                 q_by=None if r.get('p_value') is None else float(b),
                 multiplicity_family_size=len(rows),
                 unestimated_tests_retained_in_adjustment_as_one=True)


def environment_clusters(d:pd.DataFrame)->dict:
    """Outcome-blind complete-link |Spearman r| >= .8 environmental clusters.

    Means use environmental rows, not trait-specific QC. Exploratory only:
    unavailable pairs never justify merging; broad families remain separate labels.
    """
    means={}
    for p in PREDICTORS:
        if p in d:
            values=pd.to_numeric(d[p],errors='coerce').replace([np.inf,-np.inf],np.nan)
            means[p]=values.groupby(d.accepted_key.astype(str)).mean()
    e=pd.DataFrame(means);k=len(PREDICTORS);dist=np.ones((k,k));np.fill_diagonal(dist,0)
    diagnostics=[]
    for i,a in enumerate(PREDICTORS):
        for j in range(i+1,k):
            b=PREDICTORS[j];n=0;r=None
            if a in e and b in e:
                common=e[[a,b]].dropna();n=len(common)
                if n>=MIN_TAXA and (common.nunique()>1).all():
                    r=float(stats.spearmanr(common[a],common[b]).statistic)
                    dist[i,j]=dist[j,i]=np.clip(1-abs(r),0,1)
            diagnostics.append({'left':a,'right':b,'n_taxa':n,'spearman_r':r})
    Z=linkage(squareform(dist,checks=True),method='complete')
    ids=fcluster(Z,t=.200000000001,criterion='distance')
    groups={p:f'C{int(v):02}' for p,v in zip(PREDICTORS,ids)}
    return {'method':'complete_link_1_minus_abs_spearman_cut_0.2',
            'minimum_paired_taxa':MIN_TAXA,'trait_values_or_eligibility_accessed':False,
            'cluster_by_predictor':groups,'pairwise_environment_diagnostics':diagnostics,
            'predictor_columns_present':[p for p in PREDICTORS if p in d],
            'not_a_confirmatory_feature_selection':True}


def scan(d:pd.DataFrame,*,covariance:JointCovariance|None=None)->dict:
    validate_observations(d)
    d=d.sort_values('obs_id',key=lambda s:s.astype(str),kind='mergesort').reset_index(drop=True)
    cl=environment_clusters(d);rows=[]
    for trait in FEATURES:
        for predictor in PREDICTORS:
            pair=aggregate_pair(d,trait,predictor)
            common=(pair[pair.taxon.isin(covariance.taxa)].copy() if covariance is not None else pair.iloc[:0].copy())
            for panel in PANELS:
                if panel=='marginal':result=fit_pair(pair,trait)
                elif covariance is None:result={'status':'joint_covariance_not_supplied','n_taxa':0,'p_value':None,'beta_raw':None}
                elif panel=='marginal_covariance_matched':result=fit_pair(common,trait)
                else:result=fit_pair(common,trait,covariance.subset(common.taxon.tolist())) if len(common) else {'status':'no_covariance_overlap','n_taxa':0,'p_value':None,'beta_raw':None}
                rows.append({'construct':trait,'predictor':predictor,'panel':panel,
                    'environment_family':FAMILIES[predictor],
                    'environment_cluster':cl['cluster_by_predictor'][predictor],
                    'eligible_taxa_before_covariance':len(pair),**result})
    for panel in PANELS:adjust_family([r for r in rows if r['panel']==panel])
    return {'schema':'azami_species_environment_stage1_v1',
        'analysis_status':'exploratory_species_mean_scan','n_input_taxa':int(d.accepted_key.nunique()),
        'constructs':list(FEATURES),'predictors':list(PREDICTORS),
        'n_planned_trait_exposure_pairs':184,'panels':list(PANELS),
        'minimum_paired_observations_per_taxon':MIN_OBSERVATIONS,
        'minimum_taxa_per_fit':MIN_TAXA,'minimum_geographic_cells_per_taxon':None,
        'independent_confirmation':False,
        'marginal_tests_are_not_phylogenetically_or_spatially_corrected':True,
        'joint_correction_status':'not_estimated' if covariance is None else 'conditional_on_supplied_upstream_covariance',
        'covariance_parameter_uncertainty_propagated':False,
        'measurement_error_of_species_means_propagated':False,
        'privacy':'no original observations, locations, observer or taxon identifiers exported',
        'environment_clustering':cl,'results':rows}


def run(table:Path,receipt_path:Path,output:Path,*,expected_n:int=EXPECTED_N,
        covariance_path:Path|None=None)->dict:
    rec=json.loads(Path(receipt_path).read_text());digest=sha256(table)
    if (rec.get('schema')!='azami_authorized_joined_s1_v1'
        or rec.get('table_sha256')!=digest or rec.get('source_row_count')!=expected_n
        or rec.get('environment_values_in_physical_units') is not True
        or not re.fullmatch('[a-f0-9]{64}',str(rec.get('source_lineage_sha256','')))):
        raise ValueError('SOURCE_RECEIPT_OR_CHECKSUM_NOT_VERIFIED')
    d=pd.read_csv(table,dtype={'obs_id':str,'accepted_key':str},low_memory=False)
    if len(d)!=expected_n:raise ValueError('SOURCE_ROW_COUNT_MISMATCH')
    cov=None if covariance_path is None else covariance_from_payload(json.loads(Path(covariance_path).read_text()))
    out=scan(d,covariance=cov)
    out.update(source_row_count=expected_n,input_table_sha256=digest,
        source_lineage_sha256=rec['source_lineage_sha256'],
        provenance_scope='provided upstream receipt; original encrypted source and raster bytes not independently reauthenticated here',
        covariance_file_sha256=None if covariance_path is None else sha256(covariance_path))
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('table','receipt','output'):p.add_argument('--'+arg,type=Path,required=True)
    p.add_argument('--covariance',type=Path)
    a=p.parse_args();o=run(a.table,a.receipt,a.output,covariance_path=a.covariance)
    print(json.dumps({'status':o['analysis_status'],'planned_pairs':184,'output_rows':len(o['results'])}))

if __name__=='__main__':main()
