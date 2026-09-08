"""Execute a pinned slice of synthetic paired within/among calibration only."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .dependence_resampling import source_partition
from .hierarchical_ecology import spherical_basis
from .module_ecology import bootstrap_matched_module,SCALES
from .protected_artifacts import new_json,require
from .workflow import ROOT,canonical_digest,digest,text_digest

SPEC = ROOT/'analysis/v3/crossed_bootstrap_simulation_contract.json'


def generate(scenario,seed):
    rng = np.random.default_rng(seed)
    taxa,p,d = 64,9,2
    counts = np.r_[1,3,np.clip(np.rint(rng.lognormal(3.7,.4,taxa-2)),18,100).astype(int)]
    g = np.repeat(np.arange(taxa),counts)
    n = len(g)
    region = rng.integers(0,16,n)
    centers_lat = np.repeat([27.5,37.5,47.5,57.5],4)
    centers_lon = np.tile([-112.5,-52.5,27.5,117.5],4)
    side = rng.choice([-1.,1.],n)
    lat = centers_lat[region]+side+rng.uniform(-.1,.1,n)
    lon = centers_lon[region]+side+rng.uniform(-.1,.1,n)
    common = rng.normal(size=(n,1))
    x = np.sqrt(.65)*common+np.sqrt(.35)*rng.normal(size=(n,p))+rng.normal(size=(taxa,p))[g]
    b = np.column_stack([.2*x[:,0]+rng.normal(size=n),rng.normal(size=n),spherical_basis(lat,lon)])
    components = np.arange(n).astype(str)
    error = rng.normal(size=(n,d))@np.array([[1.,.6],[0.,.8]])
    if scenario!='iid_null':
        error = .5*error + rng.normal(size=(16,d))[region]@np.array([[.8,.5],[0.,.65]])
    # Duplicate measurement/source units before defining realized taxon means.
    for taxon in range(2,taxa,8):
        i,j = np.flatnonzero(g==taxon)[:2]
        components[j] = components[i]
        x[j],b[j],error[j],lat[j],lon[j] = x[i],b[i],error[i],lat[i],lon[i]
    beta_w = np.zeros((d,p)); beta_w[:,5:] = [[.25,.25,.25,.25],[.35,.35,.35,.35]]
    beta_a = beta_w.copy()
    if scenario=='scale_difference':
        beta_a[:,4] = [.45,-.3]
    deviation = (np.sqrt(.4)*rng.normal(size=(taxa,1,1))+np.sqrt(.6)*rng.normal(size=(taxa,d,p)))
    tau = np.array([.15,.9]) if scenario=='heterogeneous_slopes' else np.array([.45,.45])
    deviation *= tau[None,:,None]
    intercept = rng.normal(size=(taxa,d))@np.array([[.6,.3],[0.,.6]])
    gamma = np.linspace(.05,.25,b.shape[1])[:,None]*np.array([[1.,1.3]])
    y = np.empty((n,d))
    for taxon in range(taxa):
        i = np.flatnonzero(g==taxon)
        mean_x = x[i].mean(axis=0)
        y[i] = (x[i]-mean_x)@(beta_w+deviation[taxon]).T + mean_x@beta_a.T + b[i]@gamma + intercept[taxon]+error[i]
    truth = np.stack([beta_w,beta_a,beta_a-beta_w])
    return y,x,g,b,components,lat,lon,truth


def json_ready(value):
    if isinstance(value,np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value,(np.floating,float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value,np.integer):
        return int(value)
    if isinstance(value,dict):
        return {key:json_ready(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):
        return [json_ready(item) for item in value]
    return value


def run(out,scenario,start,end):
    spec = json.loads(SPEC.read_text())
    require(scenario in spec['scenarios'] and 0<=start<end<=spec['replicates_per_scenario'],'Slice outside specified simulation')
    require(not out.exists(),'Preserve earlier synthetic run')
    out.mkdir(parents=True)
    code = ['analysis/v3/dependence_resampling.py','analysis/v3/module_ecology.py',
            'analysis/v3/joint_partial_pooling.py','analysis/v3/hierarchical_ecology.py',
            'analysis/v3/model_design_diagnostics.py','analysis/v3/simulate_crossed_bootstrap.py']
    execution = {'specification':spec,'specification_canonical_sha256':canonical_digest(spec),
                 'implementation_sha256_text_lf':{path:text_digest(ROOT/path) for path in code},
                 'resampling_contract_canonical_sha256':canonical_digest(json.loads((ROOT/'analysis/v3/dependence_resampling_contract.json').read_text())),
                 'scenario':scenario,'start':start,'end':end,'ecological_fitting_authorized':False}
    new_json(out/'execution_contract.json',execution)
    compact = []
    with (out/'synthetic_nested_records.jsonl').open('x',encoding='utf-8',newline='\n') as records:
        for replicate in range(start,end):
            seed = spec['seed']+10000*spec['scenarios'].index(scenario)+replicate
            y,x,g,b,components,lat,lon,truth = generate(scenario,seed)
            for degrees in spec['grid_degrees']:
                started = time.perf_counter()
                row = {'scenario':scenario,'outer_replicate':replicate,'seed':seed,'grid_degrees':degrees,
                       'observations':len(y),'generating_coefficients':truth.tolist()}
                try:
                    source = source_partition(np.arange(len(y)),g,components,lat,lon,grid_degrees=degrees)
                    result = bootstrap_matched_module(y,x,b,source,np.arange(len(y)),seed=seed,
                              replicates=spec['bootstrap_replicates'],blocks=spec['process_indices'])
                    row.update(status=result['status'],planned_bootstrap_replicates=result['planned_replicates'],
                               estimable_bootstrap_replicates=result['estimable_replicates'],
                               source_taxon_blocks=len(np.unique(source.taxon_blocks)),source_spatial_blocks=len(np.unique(source.spatial_blocks)),
                               point_coefficients=result['point']['coefficients'].tolist())
                    if result['summary'] is not None:
                        summary = result['summary']
                        row['interval_coverage'] = ((summary['basic_interval_low']<=truth)&(truth<=summary['basic_interval_high'])).tolist()
                        row['candidate_tests'] = []
                        for test in summary['candidate_tests']:
                            test = test.copy()
                            block = spec['process_indices'][test['process']]
                            target = truth[SCALES.index(test['scale'])][:,block]
                            test['generating_null'] = bool(np.all(target==0))
                            test['rejects_at_point05'] = (test['candidate_probability']<.05 if test['candidate_probability'] is not None else None)
                            row['candidate_tests'].append(test)
                    records.write(json.dumps(json_ready({'identity':row,'result':result}),allow_nan=False)+'\n')
                except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
                    row.update(status='point_or_design_not_estimable',error_type=type(error).__name__,reason=str(error),
                               planned_bootstrap_replicates=spec['bootstrap_replicates'],estimable_bootstrap_replicates=0)
                    records.write(json.dumps(json_ready({'identity':row,'optimizer_attempts':getattr(error,'optimizer_attempts',[])}),allow_nan=False)+'\n')
                records.flush()
                row['elapsed_seconds'] = time.perf_counter()-started
                compact.append(row)
                new_json(out/f'summary-{replicate:03d}-{degrees}.json',row)
                print(json.dumps({key:row[key] for key in ('scenario','outer_replicate','grid_degrees','status','estimable_bootstrap_replicates','elapsed_seconds')}),flush=True)
    public = {'status':'PRELIMINARY_CROSSED_MODULE_SIMULATION_SLICE_EXECUTED_NO_ECOLOGY',
              'execution_contract':execution,'results':compact,'nested_records_sha256':digest(out/'synthetic_nested_records.jsonl'),
              'ecological_models_executed':0,'empirical_trait_environment_values_read':0,'ecological_fitting_authorized':False}
    new_json(out/'public_report.json',public)
    print(json.dumps({'status':public['status'],'planned_cases':(end-start)*2,'recorded_cases':len(compact)}),flush=True)
    return public


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--scenario',required=True)
    parser.add_argument('--start',type=int,required=True)
    parser.add_argument('--end',type=int,required=True)
    args = parser.parse_args()
    run(args.out,args.scenario,args.start,args.end)
