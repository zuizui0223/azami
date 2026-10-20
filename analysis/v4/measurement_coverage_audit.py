"""Trace inherited phyllary coverage; publish aggregate counts, never source rows."""
import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import zipfile

import pandas as pd
from reproducibility.run_current_analysis import INPUTS, restore_zenodo


def summarize_heads(heads, native_ids):
    if heads.obs_id.isna().any() or heads.annotation_unit_id.duplicated().any():
        raise ValueError('Invalid source head identity')
    out={}
    for metric in ('involucre_projection_p95','involucre_spread_fraction'):
        status=heads[metric+'_status']
        if status.isna().any() or status.astype(str).str.strip().eq('').any():
            raise ValueError('Unrecorded measurement disposition')
        groups={}
        for name,mask in [('all_measured',pd.Series(True,index=heads.index)),
                          ('current_native',heads.obs_id.isin(native_ids))]:
            d=heads.loc[mask]
            s=d[metric+'_status']
            usable=s.eq('usable')
            reasons=Counter(reason for v in s.loc[~usable]
                            for reason in set(str(v).split(';')) if reason)
            groups[name]={
                'heads':len(d),'observations':int(d.obs_id.nunique()),
                'usable_heads':int(usable.sum()),
                'usable_observations':int(d.loc[usable,'obs_id'].nunique()),
                'unusable_heads':int((~usable).sum()),
                'exclusive_status_combinations':s.value_counts().to_dict(),
                'overlapping_failure_reasons':dict(reasons)}
        out[metric]=groups
    return out


def execute(table, work, output):
    work=Path(work)
    archives,inputs=work/'archives',work/'inputs'
    restore_zenodo(archives,inputs)
    artifact,expected,_=INPUTS['continuous']
    raw=(archives/f'artifact-{artifact}-continuous.zip').read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise ValueError('Changed inherited measurement archive')
    native=pd.read_csv(table,usecols=['obs_id'],dtype=str)
    if native.obs_id.isna().any() or not native.obs_id.is_unique:
        raise ValueError('Invalid native population')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        cohort=json.loads(z.read('extended_cohort_audit.json'))
        rebuild=json.loads(z.read('rebuilt/continuous_trait_rebuild_report.json'))
        source=json.loads(z.read('extended_measurements/extended_continuous_provenance.json'))
        heads=pd.read_csv(io.BytesIO(z.read('extended_measurements/extended_continuous_head_level.csv')),
                          dtype={'obs_id':str},low_memory=False)
        source_commit=z.read('git_sha.txt').decode().strip()
    if len(heads)!=source['n_heads'] or len(heads)!=rebuild['n_rebuilt_heads']:
        raise ValueError('Published measurement stage counts disagree')
    if rebuild['n_rebuilt_heads']+rebuild['n_failed_heads']!=cohort['n_selected_heads']:
        raise ValueError('Reconstruction stage is not conserved')
    summary=summarize_heads(heads,set(native.obs_id))
    for metric,groups in summary.items():
        if groups['all_measured']['usable_heads']!=source['endpoint_usable_heads'][metric]:
            raise ValueError('Endpoint dispositions disagree with published provenance')
    result={'source_archive_sha256':expected,'source_commit':source_commit,
            'native_population_observations':len(native),
            'historical_selection':cohort,'historical_crop_reconstruction':rebuild,
            'endpoint_dispositions':summary,
            'scope':'historical measurement attrition, not new measurements',
            'full_native_measurement_completed':False,
            'claim_boundary':'An absent archived measurement is not biological absence '
                             'or evidence that an original photograph is unmeasurable.',
            'failure_reason_counts_overlap':True}
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    Path(output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('table','work','output'):
        p.add_argument('--'+name,required=True,type=Path)
    args=p.parse_args()
    execute(args.table,args.work,args.output)
