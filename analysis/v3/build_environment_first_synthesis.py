"""Reorganize pinned PR93 results by environmental gradient without new tests."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

MODULES={'presentation_angle':'presentation','floral_lightness':'colour','floral_chroma':'colour','floral_hue':'colour','head_elongation':'head_form','head_compactness':'head_form','involucre_form':'involucre','projection_prominence':'involucre','projection_pattern':'involucre'}

def build(source,out):
    out.mkdir(parents=True,exist_ok=True);tables=[]
    for scale,file in [('among_taxon','biological_axes_among_min5.csv'),('within_taxon','biological_axes_within.csv')]:
        f=pd.read_csv(source/file);assert len(f)==90 and not f.duplicated(['construct_id','predictor']).any()
        f=f[f.construct_id.isin(MODULES)].copy();assert len(f)==81
        f['scale']=scale;f['module']=f.construct_id.map(MODULES)
        f['direction']=np.where(f.kind.eq('joint'),'joint_no_single_direction',np.where(f.beta_std>0,'positive','negative'))
        tables.append(f)
    f=pd.concat(tables,ignore_index=True)
    sampling=pd.read_csv(source/'sampling.csv');spatial=pd.read_csv(source/'spatial.csv');historical=pd.read_csv(source/'historical.csv')
    f=f.merge(sampling[['scale','construct_id','predictor','sampling_composition_stability_class']],on=['scale','construct_id','predictor'],how='left',validate='one_to_one')
    f=f.merge(spatial[['scale','construct_id','predictor','broad_spatial_sensitivity_pass']],on=['scale','construct_id','predictor'],how='left',validate='one_to_one')
    historical['scale']='among_taxon'
    f=f.merge(historical[['scale','construct_id','predictor','historical_placement_sensitivity_pass']],on=['scale','construct_id','predictor'],how='left',validate='one_to_one')
    def status(r):
        if not r.fdr_0_05:return 'not_fdr_supported'
        sampling_ok=r.sampling_composition_stability_class=='stable_all_declared_scenarios'
        if pd.isna(r.sampling_composition_stability_class) or pd.isna(r.broad_spatial_sensitivity_pass):return 'sensitivity_incomplete'
        if not sampling_ok or not r.broad_spatial_sensitivity_pass:return 'sensitivity_not_passed'
        if r.scale=='within_taxon':return 'sampling_and_space_pass_phylogeny_not_applicable'
        if pd.isna(r.historical_placement_sensitivity_pass):return 'historical_not_evaluated'
        return 'full_declared_chain_pass' if r.historical_placement_sensitivity_pass else 'historical_not_passed'
    f['robustness_status']=f.apply(status,axis=1)
    f=f.sort_values(['predictor','module','construct_id','scale']);f.to_csv(out/'trait_environment_map.csv',index=False)
    rows=[]
    for (scale,predictor),one in f.groupby(['scale','predictor']):
        supported=one[one.fdr_0_05];passed=one[one.robustness_status.isin(['full_declared_chain_pass','sampling_and_space_pass_phylogeny_not_applicable'])]
        rows.append(dict(scale=scale,predictor=predictor,tested_constructs=len(one),fdr_constructs=len(supported),supported_modules='|'.join(sorted(supported.module.unique())),robust_constructs=len(passed),robust_names='|'.join(passed.construct_id)))
    pd.DataFrame(rows).to_csv(out/'gradient_overview.csv',index=False)
    report={'upstream_pr':93,'upstream_commit':'4deae0815880baf5db3e2795bd79369344e8bace','rows':len(f),
      'question':'Which continuous capitulum traits covary along environmental gradients, and which relationships withstand declared sensitivity analyses?',
      'scope':'nine non-surface biological constructs; inherited q values from ten-construct 90-test families at each scale; no recalculated FDR or new tests',
      'whole_capitulum':'compare constituent trait responses along the same gradients; shared significance alone is not syndrome evidence; use existing common-cohort integration separately',
      'boundary':'marginal associations among correlated predictors; no independent causal effects, adaptation, pigment amount or fully resolved phylogenetic correction',
      'status_counts':f.robustness_status.value_counts().to_dict(),
      'robust_rows':f[f.robustness_status.eq('full_declared_chain_pass')][['construct_id','predictor','scale','beta_std','q_bh']].to_dict('records'),
      'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.glob('*.csv'))}}
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(build(a.source,a.out),indent=2))
