"""Predeclared six-dimension extension, retaining matched reference comparisons."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from analysis.v4.environment_candidate_audit import CANDIDATES
from analysis.v4.observation_models import ENV, MODULES, fit_module, spherical_basis


def execute(table, receipt, audit_path, output):
    digest=hashlib.sha256(table.read_bytes()).hexdigest()
    receipt=json.loads(receipt.read_text())
    audit=json.loads(audit_path.read_text())
    if digest!=receipt['prepared_table_sha256'] or digest!=audit['input_sha256']:
        raise ValueError('Input and candidate audit must match')
    frame=pd.read_csv(table,low_memory=False)
    if len(frame)!=receipt['native_observations'] or not frame.obs_id.is_unique:
        raise ValueError('Native observation population changed')
    scales={}
    for v in CANDIDATES:
        frame[v]=frame[v].where(np.isfinite(frame[v]))
        scales[v]={'mean':float(frame[v].mean()),'sd':float(frame[v].std(ddof=0))}
        if not np.isfinite(scales[v]['sd']) or scales[v]['sd']<=0:
            raise ValueError('Candidate lacks finite variation: '+v)
        frame[v+'__between']=frame.groupby('accepted_key')[v].transform('mean')
        frame[v+'__within']=frame[v]-frame[v+'__between']
    ranking=frame.accepted_key.value_counts().index
    geo=spherical_basis(frame.latitude.to_numpy(float),frame.longitude.to_numpy(float))
    if not np.isfinite(geo).all() or np.any(geo.std(axis=0)<=1e-12):
        raise ValueError('Invalid broad-space support')
    adjustments=[f'seasonal_space_{j}' for j in range(8)]
    frame[adjustments]=(geo-geo.mean(axis=0))/geo.std(axis=0)
    frame=frame.loc[np.isfinite(frame[CANDIDATES]).all(axis=1)]
    models={}
    for module,outcomes in MODULES.items():
        for name,predictors,drop,space in [
            ('matched_four',ENV,0,[]),('six',CANDIDATES,0,[]),
            ('six_drop2',CANDIDATES,2,[]),('six_drop5',CANDIDATES,5,[]),
            ('six_space',CANDIDATES,0,adjustments)]:
            print(f'Fitting {module}/{name}',flush=True)
            d=frame.loc[~frame.accepted_key.isin(ranking[:drop])]
            fit=fit_module(d,outcomes,scales,predictors,space)
            fit['excluded_species_count']=drop
            models[module+'/'+name]=fit
        for name in ('matched_four','six','six_space'):
            fit=models[module+'/'+name]
            if fit['status']=='estimated':
                for endpoint,c in fit['components'].items():
                    if c['n']!=audit['populations'][endpoint]['six_environment_n']:
                        raise ValueError('Matched-row count differs from candidate audit')
    # Fixed family sizes include failed tests as p=1; do not reward missing fits.
    for name,predictors in [('matched_four',ENV),('six',CANDIDATES),
                            ('six_drop2',CANDIDATES),('six_drop5',CANDIDATES),
                            ('six_space',CANDIDATES)]:
        tests=[models[module+'/'+name].get('tests',{}).get(v+s,{})
               for module in MODULES for v in predictors for s in ('__between','__within')]
        for test,q in zip(tests,multipletests([t.get('p',1.) for t in tests],method='holm')[1]):
            if test.get('status')=='estimated':
                test['holm_p_within_specification']=float(q)
    result={'input_sha256':digest,'candidate_audit_sha256':hashlib.sha256(audit_path.read_bytes()).hexdigest(),
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scales':scales,'models':models,
            'scope':'Retrospective six-environment extension; separate Holm families per '
                    'specification. No causal identification or significance-based selection.'}
    output.mkdir(parents=True,exist_ok=True)
    (output/'seasonal_models.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# 季節性を含む6環境モデル','',
           '既存4変数の結果を見た後の改訂解析。各形質の観測は平均せず、'
           '種間・種内に環境を分解した。変数は有意性によって選択していない。',
           '各仕様内でHolm補正（4変数24検定、6変数36検定）。仕様間は独立な再現ではない。','',
           '|形質|環境・尺度|同じ観測の4変数|6変数|上位2種除外|上位5種除外|6変数＋広域空間|',
           '|---|---|---:|---:|---:|---:|---:|']
    for module in MODULES:
        for v in CANDIDATES:
            for s in ('__between','__within'):
                cells=[]
                for name in ('matched_four','six','six_drop2','six_drop5','six_space'):
                    p=models[module+'/'+name].get('tests',{}).get(v+s,{}).get('holm_p_within_specification')
                    cells.append(f'{p:.4g}' if p is not None else '対象外／評価不能')
                lines.append(f'|{module}|{v+s}|'+ '|'.join(cells)+'|')
    lines+=['','係数・点別95%区間・観測数はseasonal_models.json。'
            '空間項は低次の球面傾向であり、残差空間相関や系統の調整を保証しない。']
    (output/'seasonal_conclusions_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('table','receipt','audit','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    execute(a.table,a.receipt,a.audit,a.output)
