"""Aggregate taxon membership using the frozen S1 selection functions.
Run with PYTHONPATH pointing to the d2c45b0 checkout, never source row export.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v4.s1_s2_extreme_lego import ALL_ENV, EXTRA, FEATURES, complete_lego_cohort
from analysis.v4.s1_s2_vif5_si10_alltraits import selected_source
from analysis.v4.s1_s2_extreme_environment import retrieve, open_raster, sample_grid


def main():
    p=argparse.ArgumentParser()
    for name in ('native','out','temporary'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);a.temporary.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.native,dtype={'obs_id':str,'accepted_key':str},low_memory=False)
    assert len(d)==313350 and d.obs_id.is_unique
    receipts=[]
    for variable in EXTRA:
        archive,meta=retrieve(variable,a.temporary)
        raster=open_raster(archive,variable,a.temporary)
        d[variable],diagnostic=sample_grid(raster,variable,d)
        receipts.append({**meta,**diagnostic})
        if raster!=archive:raster.unlink()
        archive.unlink()
        print(json.dumps({'sampled':variable,'finite':diagnostic['n_finite_extracted']}),flush=True)
    cohorts={}
    for name,cols in FEATURES.items():
        part,cover=selected_source(d,cols)
        assert part is not None, (name,cover)
        cohorts[name]=part
    common=complete_lego_cohort(d)
    common=common.loc[np.isfinite(common[ALL_ENV].to_numpy(float)).all(axis=1)].copy()
    counts=common.groupby('accepted_key').size()
    common=common.loc[common.accepted_key.isin(counts[counts>=5].index)]
    assert len(cohorts['angle'])==168907
    assert len(common)==156683 and common.accepted_key.nunique()==174
    cohorts['common8']=common
    rows=[]; aggregates=[]; names=[c for c in ('taxon_name','source_taxon','species','scientific_name') if c in d]
    for name,part in cohorts.items():
        for key,g in part.groupby('accepted_key'):
            row={'cohort':name,'accepted_key':key,'n_observations':len(g)}
            for col in names:row[col]=';'.join(sorted(set(g[col].dropna().astype(str))))
            rows.append(row)
            trait_columns=list(dict.fromkeys(c for cs in FEATURES.values() for c in cs)) if name=='common8' else list(FEATURES[name])
            for coordinate in [*trait_columns,*ALL_ENV]:
                values=g[coordinate].to_numpy(float)
                assert len(values)>=5 and np.isfinite(values).all(), (name,key,coordinate)
                aggregates.append({'cohort':name,'accepted_key':key,'coordinate':coordinate,
                    'n_observations':len(values),'mean':float(values.mean()),
                    'sample_sd':float(values.std(ddof=1))})
    pd.DataFrame(rows).to_csv(a.out/'S1_FROZEN_PRIMARY_TAXON_COHORTS.csv',index=False)
    aggregate_path=a.out/'S1_FROZEN_PRIMARY_TAXON_MOMENTS.csv'
    pd.DataFrame(aggregates).to_csv(aggregate_path,index=False)
    receipt={'frozen_code':'d2c45b00018a0c04e7c69a29857207dc1ac06d0c','native_file_sha256':hashlib.sha256(a.native.read_bytes()).hexdigest(),'cohorts':{k:{'n':len(v),'taxa':v.accepted_key.nunique()} for k,v in cohorts.items()},'raster_sources':receipts,'private_observations_exported':False,'taxon_name_columns':names}
    receipt['aggregate_moments']={'sha256':hashlib.sha256(aggregate_path.read_bytes()).hexdigest(),
        'rows':len(aggregates),'method':'Arithmetic means and sample SD on each exact frozen cohort; hue retained jointly as sine/cosine coordinates, not an arithmetic mean of degrees.',
        'no_observation_ids_locations_or_users':True,'environment_columns':ALL_ENV,
        'claim_limit':'Aggregate contemporary photographic traits and niches; no historical climate or independent botanical direction states.'}
    (a.out/'S1_FROZEN_PRIMARY_TAXON_COHORTS_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')

if __name__=='__main__':main()
