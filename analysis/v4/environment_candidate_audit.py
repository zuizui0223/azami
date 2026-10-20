"""Outcome-blind environmental coverage and identifiability audit."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.observation_models import ENV, MODULES

CANDIDATES = ENV + ['chelsa_bio04', 'chelsa_bio15']
ROLES = {
    'chelsa_bio01': 'mean thermal conditions',
    'chelsa_bio12': 'water supply; not direct floral wetting',
    'chelsa_rsds_mean': 'radiation exposure; not flower-level irradiance',
    'chelsa_vpd_mean': 'atmospheric evaporative demand; coupled to temperature',
    'chelsa_bio04': 'temperature seasonality',
    'chelsa_bio15': 'precipitation seasonality',
}
DEFERRED = {
    'chelsa_sfcwind_mean': 'distinct mechanical exposure; requires separate hypothesis',
    'chelsa_gsp': 'overlapping precipitation summary; alternative exposure definition',
    'chelsa_npp': 'ecosystem response to climate; not an automatic climate confounder',
}


def matrix_summary(frame, columns):
    a=frame[columns].to_numpy(float)
    n=len(a)
    result={'n':n, 'columns':columns}
    if n<=len(columns) or not np.isfinite(a).all():
        return {**result, 'status':'insufficient_or_nonfinite'}
    sd=a.std(axis=0)
    if np.any(sd<=1e-12):
        return {**result, 'status':'constant', 'sd':sd.tolist()}
    z=(a-a.mean(axis=0))/sd
    correlation=z.T@z/n
    eigen=np.linalg.eigvalsh(correlation)
    result['correlation']=correlation.tolist()
    if eigen[0]<=1e-12:
        return {**result, 'status':'singular'}
    vif=np.diag(np.linalg.inv(correlation))
    return {**result, 'status':'estimated', 'vif':vif.tolist(),
            'condition_index':float(np.sqrt(eigen[-1]/eigen[0]))}


def audit(frame):
    # Species exposure means use the full native population, before trait masks.
    d=frame.copy()
    for v in CANDIDATES:
        d[v]=d[v].where(np.isfinite(d[v]))
        d[v+'__between']=d.groupby('accepted_key')[v].transform('mean')
        d[v+'__within']=d[v]-d[v+'__between']
    four=np.isfinite(d[ENV]).all(axis=1)
    six=np.isfinite(d[CANDIDATES]).all(axis=1)
    populations={'native_population':pd.Series(True,index=d.index)}
    for outcomes in MODULES.values():
        for name in outcomes:
            # Only eligibility/missingness, never response magnitude or association.
            populations[name]=d[name+'__eligible'].eq(True)&np.isfinite(d[name])
    groups={}
    for name,eligible in populations.items():
        rows=d.loc[eligible & six]
        groups[name]={
            'eligible_n':int(eligible.sum()),
            'four_environment_n':int((eligible & four).sum()),
            'six_environment_n':len(rows),
            'lost_from_four':int((eligible & four & ~six).sum()),
            'species_n':int(rows.accepted_key.nunique()),
            'designs':{level:matrix_summary(rows,columns) for level,columns in [
                ('raw',CANDIDATES),
                ('between_observation_weighted',[v+'__between' for v in CANDIDATES]),
                ('within',[v+'__within' for v in CANDIDATES]),
                ('full_twelve_terms',[v+s for v in CANDIDATES
                                      for s in ('__between','__within')])]},
        }
    return {'status':'candidate_audit_not_model_selection',
            'candidate_roles':ROLES, 'deferred_roles':DEFERRED,
            'finite_native_counts':{v:int(d[v].notna().sum()) for v in CANDIDATES},
            'populations':groups,
            'interpretation':'Unweighted design diagnostics, not causal adjustment '
            'certification. No automatic VIF deletion or phenotype-based selection.'}


def main():
    p=argparse.ArgumentParser()
    for name in ('table','receipt','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    receipt=json.loads(a.receipt.read_text())
    digest=hashlib.sha256(a.table.read_bytes()).hexdigest()
    if digest!=receipt['prepared_table_sha256']:
        raise ValueError('Prepared input checksum mismatch')
    result=audit(pd.read_csv(a.table,low_memory=False))
    result['input_sha256']=digest
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
