"""V4 S1 taxon-level multivariate RDA variation partition: climate and photograph sampling geography.

Exploratory two-block partial RDA, not an independent causal / adaptive inference.
PHYLOGENY IS NOT INCLUDED: it requires a separately authenticated V4-matched tree.
All original observations stay within the authorised original-source environment;
only anonymized aggregate statistics, dimensions, PC source weights are exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import lstsq

from analysis.v4.species_environment_spls2 import (
    aggregate_species_joint, CLIMATE, TRAITS, YCOLS, OUTCOMES, MIN_TAXA
)
from analysis.v4.species_environment_existing_inputs import load_existing_inputs

GEO_COLS=('geo_x','geo_y','geo_z','geo_xy','geo_xz','geo_yz','geo_x2_minus_y2','geo_3z2_minus1')
PRIMARY_PCS=6
SENSITIVITY_PCS=(4,8,22)
SEED=20261107


def spherical_geo(latitude:np.ndarray, longitude:np.ndarray)->np.ndarray:
    lat=np.asarray(latitude,float);lon=np.asarray(longitude,float)
    if lat.shape!=lon.shape or lat.ndim!=1 or not (np.isfinite(lat).all() and np.isfinite(lon).all()) or (np.abs(lat)>90).any() or (np.abs(lon)>180).any():
        raise ValueError('ORIGINAL_S1_COORDINATES_MISSING_OR_INVALID')
    la=np.deg2rad(lat);lo=np.deg2rad(lon)
    x=np.cos(la)*np.cos(lo);y=np.cos(la)*np.sin(lo);z=np.sin(la)
    return np.column_stack([x,y,z,x*y,x*z,y*z,x*x-y*y,3*z*z-1])


def aggregate_joint_geo(d:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    """Only the exact same photo-eligible cohort as original joint sPLS2."""
    required={'obs_id','accepted_key','latitude','longitude',*YCOLS,
              *(x+'__eligible' for x in YCOLS)}
    if not required <= set(d):raise ValueError('ORIGINAL_S1_COORDINATES_OR_TRAITS_MISSING')
    if (d.obs_id.isna().any() or not d.obs_id.is_unique or d.accepted_key.isna().any()):
        raise ValueError('ORIGINAL_S1_IDENTITIES_INVALID')
    d=d.sort_values('obs_id',key=lambda s:s.astype(str),kind='mergesort').reset_index(drop=True)
    means,audit=aggregate_species_joint(d)
    Y=d[list(YCOLS)].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    eligible=np.isfinite(Y).all(axis=1)
    for c in YCOLS:eligible &= d[c+'__eligible'].eq(True).fillna(False).to_numpy(bool)
    hue=OUTCOMES['hue'];eligible &= np.linalg.norm(Y[:,list(hue)],axis=1)>=.20
    if int(eligible.sum())!=audit['all_eight_joint_photo_eligible']:
        raise ValueError('JOINT_SPLS2_PHOTO_COHORT_CHANGED')
    points=d.loc[eligible,['accepted_key','latitude','longitude']].copy()
    z=points[['latitude','longitude']].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    vec=spherical_geo(z[:,0],z[:,1])
    geo=pd.DataFrame(vec,columns=GEO_COLS)
    geo['_taxon_private']=points.accepted_key.astype(str).to_numpy()
    aggregate=geo.groupby('_taxon_private',sort=True,observed=True).mean()
    geo_counts=geo.groupby('_taxon_private',sort=True,observed=True).size()
    out=means.merge(aggregate,left_on='_taxon_private',right_index=True,how='left',validate='one_to_one',sort=False)
    if not np.isfinite(out[list(GEO_COLS)].to_numpy(float)).all():
        raise ValueError('ORIGINAL_S1_COORDINATE_AGGREGATION_FAILED')
    check=geo_counts.reindex(out._taxon_private).to_numpy(int)
    if not np.array_equal(check,out.n_photos.to_numpy(int)):
        raise ValueError('GEOGRAPHIC_AND_PHENOTYPIC_SUPPORT_NOT_IDENTICAL')
    audit['geography']='eight_spherical_basis_terms_averaged_across_joint_QC_photos'
    audit['photo_geo_and_trait_matching_verified']=True
    return out,audit


def make_response_matrix(frame:pd.DataFrame, cols:list[str])->np.ndarray:
    """Each scalar trait variance 1; circular hue sin+cos *joint* variance 1."""
    cols=list(cols)
    y=frame[cols].to_numpy(float)
    if y.ndim!=2 or not np.isfinite(y).all() or y.shape[0]<MIN_TAXA:
        raise ValueError('JOINT_RDA_OUTCOMES_NOT_SUPPORTED')
    mu=y.mean(axis=0);sd=y.std(axis=0,ddof=1)
    if (sd<1e-12).any():raise ValueError('JOINT_RDA_CONSTANT_PHENOTYPIC_COMPONENT')
    for h in (('corolla_hue_sin','corolla_hue_cos'),):
        if set(h)<=set(cols):
            idx=[cols.index(j) for j in h]
            sd[idx]=np.sqrt(np.sum(sd[idx]**2))
    return (y-mu)/sd


def _xscaled(matrix:np.ndarray)->np.ndarray:
    x=np.asarray(matrix,float)
    if x.ndim!=2 or x.shape[0]<MIN_TAXA or not np.isfinite(x).all():
        raise ValueError('PREDICTOR_MATRIX_INVALID')
    sd=x.std(axis=0,ddof=1)
    if (sd<1e-13).any():raise ValueError('PREDICTOR_COLUMN_CONSTANT')
    return (x-x.mean(axis=0))/sd


def climate_pcs(frame:pd.DataFrame,k:int)->tuple[np.ndarray,dict]:
    """PCA uses X only; no choice of PC count from Y or model significance."""
    if not 1<=k<=len(CLIMATE):raise ValueError('INVALID_CLIMATE_PC_COUNT')
    X=_xscaled(frame[list(CLIMATE)].to_numpy(float))
    u,s,vt=np.linalg.svd(X,full_matrices=False)
    variance=s*s/(len(X)-1);fraction=variance/variance.sum()
    components=vt[:k].copy()
    for j in range(k):
        pivot=np.argmax(np.abs(components[j]));
        if components[j,pivot]<0:components[j]*=-1
    scores=X@components.T
    out={'n_registered_climate_axes':len(CLIMATE),
         'variance_ratio':[float(v) for v in fraction[:k]],
         'cumulative_variance_explained':float(fraction[:k].sum()),
         'learned_from_traits':False,
         'principal_component_loadings':[
             {'pc':int(j+1),'largest_absolute_X_loadings':[
                 {'variable':CLIMATE[int(i)],'loading':float(components[j,i])}
                 for i in np.argsort(-np.abs(components[j]),kind='stable')[:6]]}
             for j in range(k)]}
    return scores,out


def _fit(Y:np.ndarray,X:np.ndarray)->dict:
    n=len(Y);A=np.column_stack((np.ones(n),X))
    coef,_,rank,_=np.linalg.lstsq(A,Y,rcond=1e-12)
    fit=A@coef;resid=Y-fit;ss=float(np.square(resid).sum())
    total=float(np.square(Y-Y.mean(axis=0)).sum())
    if total<1e-12:raise ValueError('NO_RESPONSE_VARIATION')
    p=int(rank-1)
    if p>=n-2:raise ValueError('INSUFFICIENT_PARTIAL_RDA_RESIDUAL_DF')
    r2=1-ss/total
    adj=1-(ss/total)*(n-1)/(n-p-1)
    return {'ss':ss,'total_ss':total,'r2':float(r2),'adjR2':float(adj),
            'model_rank':int(rank),'predictor_rank':p,'fit':fit,'residual':resid}


def _permutation_p(Y:np.ndarray, full:np.ndarray, reduced:np.ndarray,
                   *,n_perm:int,seed:int, sector:np.ndarray|None)->float:
    base=_fit(Y,reduced);bfull=_fit(Y,full)
    difference=bfull['predictor_rank']-base['predictor_rank']
    if difference<=0:raise ValueError('ZERO_UNIQUE_CLIMATE_DESIGN_RANK')
    observed=max(0.,base['ss']-bfull['ss'])/difference/(bfull['ss']/(len(Y)-bfull['model_rank']))
    rng=np.random.default_rng(seed);ge=0
    sector=None if sector is None else np.asarray(sector,int)
    if sector is not None and (sector.shape!=(len(Y),) or len(np.unique(sector))<2):
        raise ValueError('INVALID_PERMUTATION_SECTOR_SUPPORT')
    for _ in range(n_perm):
        order=np.arange(len(Y))
        if sector is None:order=rng.permutation(len(Y))
        else:
            for g in np.unique(sector):
                ix=np.flatnonzero(sector==g);order[ix]=rng.permutation(ix)
        trial=base['fit']+base['residual'][order]
        null=_fit(trial,reduced);whole=_fit(trial,full)
        stat=max(0.,null['ss']-whole['ss'])/difference/(whole['ss']/(len(Y)-whole['model_rank']))
        ge+=stat>=observed-1e-12
    return float((ge+1)/(n_perm+1))


def varpart(Y:np.ndarray,C:np.ndarray,G:np.ndarray,*,n_permutations:int=499,
            seed:int=SEED,sector:np.ndarray|None=None)->dict:
    Y=np.asarray(Y,float); C=np.asarray(C,float);G=np.asarray(G,float)
    if (Y.ndim!=2 or C.ndim!=2 or G.ndim!=2 or len(Y)!=len(C) or len(Y)!=len(G)
            or not np.isfinite(Y).all() or not np.isfinite(C).all() or not np.isfinite(G).all()
            or n_permutations<9):raise ValueError('PARTIAL_RDA_INPUT_INVALID')
    m=_fit(Y,np.column_stack([C,G]));c=_fit(Y,C);g=_fit(Y,G)
    # Three-block (climate + space + history) is NOT estimated without a matched tree.
    shared=c['adjR2']+g['adjR2']-m['adjR2']
    base={'n_taxa':int(len(Y)),
          'n_climate_dimensions':int(C.shape[1]),'n_geo_dimensions':int(G.shape[1]),
          'rank_full_model':m['predictor_rank'],
          'rawR2_climate_only':c['r2'],'rawR2_geography_only':g['r2'],
          'rawR2_full':m['r2'],
          'adjR2_climate_only':c['adjR2'],'adjR2_geography_only':g['adjR2'],
          'adjR2_joint_climate_geography':m['adjR2'],
          'adjR2_unique_climate':float(m['adjR2']-g['adjR2']),
          'adjR2_unique_geography':float(m['adjR2']-c['adjR2']),
          'adjR2_shared':float(shared),'adjR2_unexplained':float(1-m['adjR2']),
          'negative_adjusted_fractions_retained':True,
          'permutation_count':int(n_permutations),
          'exploratory_p_climate_beyond_geo':_permutation_p(Y,np.column_stack([C,G]),G,n_perm=n_permutations,seed=seed,sector=None),
          'exploratory_p_climate_beyond_geo_sector_restricted':None if sector is None else _permutation_p(
              Y,np.column_stack([C,G]),G,n_perm=n_permutations,seed=seed+153,sector=sector),
          'p_values_are_posthoc_not_corrected_for_spatial_residual_dependence':True}
    return base


def analyze(frame:pd.DataFrame,*,n_permutations:int=499,seed:int=SEED)->dict:
    taxa,audit=aggregate_joint_geo(frame)
    n_before=len(taxa)
    complete=np.isfinite(taxa[list(CLIMATE)].to_numpy(float)).all(axis=1)
    taxa=taxa.loc[complete].copy()
    if len(taxa)<MIN_TAXA:raise ValueError('CLIMATE22_JOINT_PHENOTYPE_COHORT_TOO_SMALL')
    y=make_response_matrix(taxa,list(YCOLS));sector=taxa.sector.to_numpy(int)
    if (sector<0).any() or set(sector)!={0,1,2} or min(np.bincount(sector,minlength=3))<8:
        raise ValueError('JOINT_S1_THREE_SECTOR_SUPPORT_UNAVAILABLE')
    G=_xscaled(taxa[list(GEO_COLS)].to_numpy(float))
    families={'all_eight':list(range(len(YCOLS))),
              'color':[YCOLS.index(p) for k in ('lightness','chroma','hue') for p in
                       (('corolla_hue_sin','corolla_hue_cos') if k=='hue' else (('corolla_lab_lightness',) if k=='lightness' else ('corolla_lab_chroma',)))],
              'angle':[YCOLS.index('orientation_image_vertical_angle')],
              'outline':[YCOLS.index(p) for p in
                        ('capitulum_outline_aspect_ratio','capitulum_outline_circularity',
                         'capitulum_outline_solidity','capitulum_width_profile_cv')]}
    x4,diag4=climate_pcs(taxa,4);x6,diag6=climate_pcs(taxa,6);x8,diag8=climate_pcs(taxa,8)
    X22=_xscaled(taxa[list(CLIMATE)].to_numpy(float))
    designs=[('main_6_climate_pcs',x6),('sensitivity_4_climate_pcs',x4),
             ('sensitivity_8_climate_pcs',x8),('sensitivity_22_original_climate_axes',X22)]
    out={}
    for f,indices in families.items():
        out[f]={}
        for j,(label,C) in enumerate(designs):
            out[f][label]=varpart(y[:,indices],C,G,n_permutations=n_permutations,
                                 seed=seed+j*1000+len(indices),sector=sector)
    meta={k:v for k,v in audit.items() if k not in ('predictor_support',)}
    return {'schema':'azami_v4_s1_species_partial_rda_climate_vs_geography_v1',
        'analysis_status':'exploratory_partial_RDA_source_authentication_required',
        'n_original_observations':int(len(frame)),
        'n_joint_photo_supported_before_all22_climate_complete_case':int(n_before),
        'n_joint_supported_taxa':int(len(taxa)),
        'photo_source_support':meta,
        'model':'equal_taxon_multivariate_RDA_partition_with_adjusted_R2',
        'primary_climate_dimension':'6_X_only_PC_axes',
        'sensitivity_climate_dimensions':[4,8,22],
        'same_exact_taxa_in_all_comparisons':True,
        'sector_counts':{str(i):int((sector==i).sum()) for i in (0,1,2)},
        'climate_source_variables':list(CLIMATE),
        'climate_PCA_diagnostics':{'4':diag4,'6':diag6,'8':diag8},
        'geography_terms':'fixed_eight_spherical_coordinate_functions_photo_mean',
        'historical_phylogeny_status':'not_estimated_no_authenticated_V4_matched_tree',
        'no_cross_scale_legacy_tree_substitution':True,
        'independent_replication':False,
        'causal_environment_effects_identified':False,
        'multivariate_permutation_assumptions_are_verified':False,
        'permutations':'exploratory_Freedman_Lane_unrestricted_and_3_sector_restricted_null',
        'not_an_independent_phylogeny_or_species_plant_measurement':True,
        'individual_observation_taxon_camera_location_not_exported':True,
        'families':out}


def run_existing_partial_rda(native:Path,decoded_receipt:Path,seven:Path,nine:Path,
                             nine_manifest:Path,source_receipt:Path,output:Path,*,
                             source_receipt_sha256:str,n_permutations:int=499,
                             expected_n:int=313350)->dict:
    inputs=(native,decoded_receipt,seven,nine,nine_manifest,source_receipt)
    if Path(output).resolve() in {Path(p).resolve() for p in inputs}:
        raise ValueError('OUTPUT_OVERWRITES_SOURCE')
    d,receipt=load_existing_inputs(native,decoded_receipt,seven,nine,source_receipt,
                                   nine_manifest,source_receipt_sha256=source_receipt_sha256,
                                   expected_n=expected_n)
    loc=pd.read_csv(native,usecols=['obs_id','latitude','longitude'],dtype={'obs_id':str},low_memory=False)
    if len(loc)!=len(d) or loc.obs_id.isna().any() or not loc.obs_id.is_unique:
        raise ValueError('ORIGINAL_S1_GEOGRAPHY_IDENTITIES_INVALID')
    ids=d.obs_id.copy()
    d=d.merge(loc,on='obs_id',validate='one_to_one',how='left',sort=False)
    if not d.obs_id.equals(ids):raise ValueError('ORIGINAL_S1_GEOGRAPHY_JOIN_ORDER_CHANGED')
    out=analyze(d,n_permutations=n_permutations)
    out.update(analysis_status='V4_SOURCE_VERIFIED_EXPLORATORY_PARTIAL_RDA',
               original_S1_sha256=receipt['input_digests']['native'],
               source_receipt_sha256=source_receipt_sha256,
               previous_sPLS2_and_GEB_results_unchanged=True)
    forbidden={'obs_id','accepted_key','_taxon_private','taxon_name','user_id',
               'photo_id','photographer_id','latitude','longitude','individual_prediction'}
    def audit(item):
        if isinstance(item,dict):
            for k,v in item.items():
                if k in forbidden:raise ValueError('PRIVATE_RECORD_KEY_IN_RDA_OUTPUT_'+k)
                audit(v)
        elif isinstance(item,list):
            for v in item:audit(v)
    audit(out)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ('native','decoded-receipt','seven','nine','nine-manifest','source-receipt','output'):
        parser.add_argument('--'+arg,type=Path,required=True)
    parser.add_argument('--source-receipt-sha256',required=True)
    parser.add_argument('--permutations',type=int,default=499)
    args=vars(parser.parse_args())
    args['n_permutations']=args.pop('permutations')
    x=run_existing_partial_rda(**args)
    print(json.dumps({'status':x['analysis_status'],
          'n_source_taxa':x['n_joint_supported_taxa'],
          'climate_unique_adjR2':x['families']['all_eight']['main_6_climate_pcs']['adjR2_unique_climate']}))

if __name__=='__main__':main()
