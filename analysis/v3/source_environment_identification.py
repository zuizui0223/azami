"""Trait-blind calendar/spatial identification audit on the native source.

This is not the realized image-measurement cohort and lacks imaging covariates.
It does not change the already selected process variables or execute trait fits.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .environment_model import definition
from .hierarchical_ecology import SPATIAL_TERMS, spherical_basis
from .model_design_diagnostics import diagnose
from .protected_artifacts import new_json, require
from .recover_native_source_authority import private_directory
from .workflow import digest, text_digest


def run(matrix: Path,out: Path):
    spec=definition()
    require(digest(matrix)==spec['source_qc_matrix_sha256'],'Wrong source-QC numerical matrix')
    out=private_directory(out)
    require(not out.exists(),'Preserve prior source identification audit')
    out.mkdir(parents=True)
    columns=['obs_id','accepted_key','native_range_status','analysis_latitude','analysis_longitude',
             'observed_year','sin_doy','cos_doy','south_indicator','south_sin','south_cos',*spec['variables']]
    source=pd.read_csv(matrix,usecols=columns,dtype={'obs_id':str,'accepted_key':str})
    require(source['obs_id'].is_unique and source['native_range_status'].eq('native').all(),'Source membership differs')
    calendar=['sin_doy','cos_doy','south_indicator','south_sin','south_cos']
    required=[*spec['variables'],'analysis_latitude','analysis_longitude','observed_year',*calendar]
    finite=np.isfinite(source[required].to_numpy(float)).all(axis=1)
    data=source.loc[finite].copy().reset_index(drop=True)
    data['year_decades_since_2000']=(data['observed_year']-2000.0)/10.0
    data[list(SPATIAL_TERMS)]=spherical_basis(data['analysis_latitude'],data['analysis_longitude'])
    data['audit_weight']=1.0/data.groupby('accepted_key')['accepted_key'].transform('size')
    nuisance=[*calendar,'year_decades_since_2000',*SPATIAL_TERMS]
    result=diagnose(data,spec['variables'],nuisance,weight_column='audit_weight',
                    cohort_id='native_source_complete_retained_exposures_calendar_space_no_imaging')
    membership=out/'source_membership_private.csv'
    data[['obs_id','accepted_key']].to_csv(membership,index=False,mode='x',lineterminator='\n')
    report={'status':'NATIVE_SOURCE_EXPOSURE_IDENTIFICATION_AUDITED_NO_TRAITS',
            'source_matrix_sha256':digest(matrix),'source_rows':len(source),
            'excluded_missing_required_exposure_or_calendar':int((~finite).sum()),
            'membership_sha256':digest(membership),
            'implementation_sha256_text_lf':text_digest(Path(__file__)),
            'diagnostics':result,'source_variables_reselected':False,
            'trait_values_read':0,'ecological_models_executed':0,
            'limits':['This source-only audit is not the final module-specific design or an ecological result.',
                      'Calendar and spatial terms are included, but endpoint-matched imaging covariates and realized model weights are absent.',
                      'Per-taxon local nuisance-adjusted rank is conservative, not the rank of the joint shared-nuisance model.',
                      'Pre-1980 dates are retained; climatology is not observation-year weather. No missing value is imputed.']}
    new_json(out/'public_report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='diagnostics'}),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    run(args.matrix,args.out)
