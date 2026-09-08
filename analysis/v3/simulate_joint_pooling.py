"""Run the complete pinned synthetic hierarchy check, never empirical fitting."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np

from .hierarchical_ecology import spherical_basis
from .joint_partial_pooling import fit_joint_partial_pooling
from .protected_artifacts import new_json, require
from .workflow import ROOT, canonical_digest, digest, text_digest

CONTRACT = ROOT/'analysis/v3/joint_pooling_simulation_contract.json'


def generate(scenario, seed, taxa=32):
    rng = np.random.default_rng(seed)
    p = 9 if scenario=='iid_nine_predictors' else 2
    counts = np.r_[1,3,np.clip(np.rint(rng.lognormal(3.6,.7,taxa-2)),12,180).astype(int)]
    g = np.repeat(np.arange(taxa),counts)
    n = len(g)
    cell = rng.integers(0,24,n)
    lat_centers, lon_centers = rng.uniform(-65,65,24), rng.uniform(-170,170,24)
    lat = lat_centers[cell]+rng.normal(0,.3,n)
    lon = lon_centers[cell]+rng.normal(0,.3,n)
    factor = rng.normal(size=n)
    x = np.sqrt(.65)*factor[:,None]+np.sqrt(.35)*rng.normal(size=(n,p))
    x += .3*np.sin(np.radians(lat))[:,None]
    x[g==1,1] = x[g==1,0]
    nuisance = np.column_stack([.4*x[:,0]+rng.normal(size=n),rng.normal(size=n),spherical_basis(lat,lon)])
    mu = np.r_[0.,np.full(p-1,.3)]
    random = rng.normal(0,.45,size=(taxa,p))
    error = rng.normal(size=n)
    if scenario=='spatially_correlated_errors_two_predictors':
        error = np.sqrt(.35)*error+np.sqrt(.65)*rng.normal(size=24)[cell]
    y = rng.normal(size=taxa)[g]+np.sum(x*(mu+random[g]),axis=1)+nuisance@np.linspace(.1,.4,nuisance.shape[1])+error
    return y,x,g,nuisance


def wilson(successes,total):
    if not total:
        return [None,None]
    z = 1.959963984540054
    p = successes/total
    center = (p+z*z/(2*total))/(1+z*z/total)
    half = z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    return [float(center-half),float(center+half)]


def run(out):
    spec = json.loads(CONTRACT.read_text(encoding='utf-8'))
    out = out.resolve()
    require(out.is_relative_to(ROOT/'local_data') and not out.exists(),'Use a new ignored synthetic execution directory')
    out.mkdir(parents=True)
    execution = {'specification':spec,'specification_canonical_sha256':canonical_digest(spec),
                 'pooling_contract_canonical_sha256':canonical_digest(json.loads((ROOT/'analysis/v3/joint_partial_pooling_contract.json').read_text())),
                 'implementation_sha256_text_lf':text_digest(ROOT/'analysis/v3/joint_partial_pooling.py'),
                 'simulation_sha256_text_lf':text_digest(Path(__file__))}
    new_json(out/'execution_contract.json',execution)
    records = []
    with (out/'synthetic_replicates.jsonl').open('x',encoding='utf-8',newline='\n') as saved:
        for s,scenario in enumerate(spec['scenarios']):
            for replicate in range(spec['replicates_per_scenario']):
                seed = spec['seed']+10000*s+replicate
                data = generate(scenario,seed,spec['taxa'])
                row = {'scenario':scenario,'replicate':replicate,'seed':seed,'observations':len(data[0]),'planned_taxa':spec['taxa']}
                try:
                    fit = fit_joint_partial_pooling(*data)
                    se = float(np.sqrt(fit.fixed_covariance[0,0]))
                    require(np.isfinite(se) and se>0,'Undefined conditional SE')
                    row.update(status='estimated',retained_taxa=len(fit.taxa),hypermean=float(fit.hypermean[0]),
                               conditional_se=se,tau=float(np.sqrt(fit.tau2[0])),
                               rejects_zero=bool(abs(fit.hypermean[0]/se)>1.959963984540054),
                               optimizer_attempts=fit.optimizer_attempts)
                except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
                    row.update(status='not_estimable',error_type=type(error).__name__,error_message=str(error),
                               optimizer_attempts=getattr(error,'optimizer_attempts',[]))
                records.append(row)
                saved.write(json.dumps(row,allow_nan=False)+'\n'); saved.flush()
                if (replicate+1)%20==0:
                    print(json.dumps({'scenario':scenario,'replicates_recorded':replicate+1}),flush=True)
    summaries = []
    for scenario in spec['scenarios']:
        rows = [r for r in records if r['scenario']==scenario]
        success = [r for r in rows if r['status']=='estimated']
        rejections = sum(r['rejects_zero'] for r in success)
        summaries.append({'scenario':scenario,'planned_replicates':len(rows),'estimable_replicates':len(success),
                          'unestimable_replicates':len(rows)-len(success),'zero_rejections':rejections,
                          'rejection_rate_among_estimable':rejections/len(success) if success else None,
                          'rejection_rate_wilson95':wilson(rejections,len(success)),
                          'coverage_among_estimable':1-rejections/len(success) if success else None,
                          'hypermean_bias':float(np.mean([r['hypermean'] for r in success])) if success else None,
                          'mean_tau':float(np.mean([r['tau'] for r in success])) if success else None,
                          'all_taxa_retained_in_every_estimable_replicate':all(r['retained_taxa']==spec['taxa'] for r in success)})
    report = {'status':'JOINT_POOLING_SYNTHETIC_REFERENCE_CHECK_EXECUTED_NOT_SPATIAL_INFERENCE_ADMISSION',
              'execution_contract':execution,'replicate_records':len(records),
              'replicate_file_sha256':digest(out/'synthetic_replicates.jsonl'),'summaries':summaries,
              'ecological_models_executed':0,'empirical_trait_environment_values_read':0,
              'ecological_fitting_authorized':False,'limits':spec['scope']}
    new_json(out/'public_report.json',report)
    print(json.dumps(report),flush=True)
    return report


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    run(parser.parse_args().out)
