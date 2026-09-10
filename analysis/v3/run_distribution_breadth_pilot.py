"""Exploratory full-range density breadth and spatial-held-out prediction pilot.

Reuses frozen v2 values. No new pixels, native filter, or PR93 integration tests.
KDE volumes describe shared PCA projections, not true biological trait-space volume.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.spatial.distance import cdist, pdist
from scipy.stats import spearmanr

ENV = ['chelsa_bio12', 'chelsa_rsds_mean', 'chelsa_vpd_mean', 'chelsa_sfcwind_mean']
GROUPS = {
    'orientation': ['orientation_image_vertical_angle'],
    'colour': ['corolla_lab_lightness', 'corolla_lab_chroma', 'corolla_hue_sin', 'corolla_hue_cos'],
    'outline': ['capitulum_outline_aspect_ratio', 'capitulum_outline_circularity', 'capitulum_outline_solidity', 'capitulum_width_profile_cv'],
    'involucre': ['involucre_length_width_ratio', 'involucre_apical_taper_ratio', 'involucre_basal_taper_ratio', 'bract_projection_roughness', 'bract_projection_p95', 'bract_projection_maximum', 'bract_spread_fraction', 'bract_projection_peak_density', 'bract_projection_asymmetry'],
}
GROUPS['whole'] = sum(GROUPS.values(), [])


def rng_for(*parts):
    return np.random.default_rng(int.from_bytes(hashlib.sha256('|'.join(map(str, parts)).encode()).digest()[:8], 'little'))


def weighted_scale(x, weights):
    if not np.isfinite(x).all() or not np.isfinite(weights).all():
        raise ValueError('Non-finite coordinate or weight before scaling')
    mu = np.average(x, axis=0, weights=weights)
    sd = np.sqrt(np.average((x-mu)**2, axis=0, weights=weights))
    if np.any(sd <= 0):
        raise ValueError('Constant coordinate')
    return mu, sd


def log_kde(points, sample, h):
    d = sample.shape[1]
    return logsumexp(-cdist(points, sample, 'sqeuclidean')/(2*h*h), axis=1) - np.log(len(sample)) - d*np.log(h*np.sqrt(2*np.pi))


def kde_volume(sample, h, rng, draws):
    # ponytail: shared isotropic kernels; estimate 90%-probability volume by Monte Carlo.
    # Kernels are unbounded in a shared projection. Do not interpret as anatomical support.
    d = sample.shape[1]
    a = sample[rng.integers(len(sample), size=draws)] + rng.normal(size=(draws, d))*h
    cutoff = np.quantile(log_kde(a, sample, h), .10)
    b = sample[rng.integers(len(sample), size=draws)] + rng.normal(size=(draws, d))*h
    density = log_kde(b, sample, h)
    return float(np.mean(np.where(density >= cutoff, np.exp(-density), 0)))


def spatial_basis(lat, lon):
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    x, y, z = np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)
    return np.column_stack([np.ones(len(x)), x, y, z, x*y, x*z, y*z, x*x, y*y])


def geographic_rms(lat, lon):
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    xyz = np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
    angular = 2*np.arcsin(np.clip(pdist(xyz)/2, 0, 1))
    return float(np.sqrt(np.mean((6371.0088*angular)**2)))


def predictive_gain(frame, y, e, group, taxon, seed):
    if len(frame) < 100:
        return None
    cells = list(zip(np.floor(frame.latitude/5).astype(int), np.floor(frame.longitude/5).astype(int)))
    unique = sorted(set(cells))
    if len(unique) < 10:
        return None
    rng = rng_for(seed, group, taxon, 'folds')
    order = rng.permutation(len(unique))
    mapping = {unique[index]: int(i % 5) for i, index in enumerate(order)}
    folds = np.array([mapping[cell] for cell in cells])
    s = spatial_basis(frame.latitude.to_numpy(), frame.longitude.to_numpy())
    full = np.column_stack([s, e])
    baseline_error = augmented_error = 0.
    ranks = []
    for k in range(5):
        train, test = folds != k, folds == k
        if train.sum() <= full.shape[1] + 5 or not test.any():
            return None
        b0, _, _, _ = np.linalg.lstsq(s[train], y[train], rcond=1e-8)
        b1, _, rank, _ = np.linalg.lstsq(full[train], y[train], rcond=1e-8)
        ranks.append(int(rank))
        baseline_error += float(np.sum((y[test]-s[test]@b0)**2))
        augmented_error += float(np.sum((y[test]-full[test]@b1)**2))
    return dict(group=group, taxon_name=taxon, observations=len(frame), geographic_cells=len(unique),
                min_training_rank=min(ranks), design_columns=full.shape[1],
                heldout_gain=1-augmented_error/baseline_error if baseline_error > 0 else None,
                spatial_sse=baseline_error, environment_plus_spatial_sse=augmented_error)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--traits', type=Path, required=True)
    p.add_argument('--environment', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--replicates', type=int, default=30)
    p.add_argument('--draws', type=int, default=512)
    p.add_argument('--seed', type=int, default=20260910)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    contract = dict(status='exploratory_pilot_not_confirmatory', sample_size=50, replicates=args.replicates,
                    seed=args.seed, draws=args.draws, groups=GROUPS, environment=ENV,
                    native_filter=False, cohort='frozen v2 spatially thinned full range',
                    shared_pca='equal taxon weighted; 90% target, maximum 3 axes, no whitening; report retained fraction',
                    hue='sin/cos share one pooled scale preserving chord geometry',
                    bandwidth='50**(-1/(dimension+4)); same within module for all taxa and resamples',
                    volume='90% probability KDE projection volume; unbounded kernels, not anatomical support',
                    breadth='environment RMS pairwise distance in four global standardized exposures; geography RMS great-circle distance',
                    uncertainty='taxon bootstrap for breadth correlation; image rarefaction is not a spatial independence correction',
                    breadth_gate='at least 15 eligible taxa; fewer descriptive only',
                    prediction='>=100 observations and >=10 occupied 5-degree cells; five held-out cell folds; compare spatial trend against same trend plus all four environments',
                    limitations=['sampled range not complete species range', 'no phylogenetic correction of breadth association', 'KDE bandwidth sensitivity pending', 'no causal plasticity or adaptation claim', 'predictive gain is not a significance test'])
    contract['sha256'] = {key: hashlib.sha256(path.read_bytes()).hexdigest() for key,path in [('traits',args.traits),('environment',args.environment)]}
    if contract['sha256']['environment'] != 'e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7':
        raise ValueError('Not the frozen full-range environment')
    (args.out/'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    t = pd.read_csv(args.traits, usecols=['obs_id','taxon_name','endpoint_id','value','measurement_available'], dtype={'obs_id':str})
    e = pd.read_csv(args.environment, dtype={'obs_id':str})
    assert len(e)==46276 and e.taxon_name.nunique()==259 and not e.obs_id.duplicated().any()
    t = t[t.measurement_available.astype(str).str.lower().isin(['true','1']) & np.isfinite(t.value)]
    assert t.endpoint_id.nunique()==22 and not t.duplicated(['obs_id','endpoint_id']).any()
    assert set(t.obs_id) <= set(e.obs_id)
    wide = t.pivot(index=['obs_id','taxon_name'], columns='endpoint_id', values='value').reset_index()
    e_complete = e.replace([np.inf,-np.inf],np.nan).dropna(subset=[*ENV,'latitude','longitude'])
    ew = 1/e_complete.groupby('taxon_name').taxon_name.transform('size').to_numpy(float)
    em, es = weighted_scale(e_complete[ENV].to_numpy(), ew)
    inventory, bases, rare, predictions = [], [], [], []
    for group, members in GROUPS.items():
        f = wide[['obs_id','taxon_name',*members]].dropna().merge(e_complete, on=['obs_id','taxon_name'], validate='one_to_one').reset_index(drop=True)
        counts = f.groupby('taxon_name').size()
        inventory.extend(dict(group=group,taxon_name=taxon,observations=int(n),eligible_50=bool(n>=50)) for taxon,n in counts.items())
        w = 1/f.groupby('taxon_name').taxon_name.transform('size').to_numpy(float)
        raw = f[members].to_numpy(float)
        mu, sd = weighted_scale(raw,w)
        if 'corolla_hue_sin' in members:
            si, ci = members.index('corolla_hue_sin'), members.index('corolla_hue_cos')
            sd[si] = sd[ci] = np.sqrt((sd[si]**2+sd[ci]**2)/2)
        z = (raw-mu)/sd
        eigen, vectors = np.linalg.eigh(z.T@(z*(w/w.sum())[:,None]))
        order = np.argsort(eigen)[::-1]; eigen, vectors = eigen[order], vectors[:,order]
        fraction = eigen/eigen.sum(); d = min(3,int(np.searchsorted(np.cumsum(fraction),.90)+1))
        scores = z@vectors[:,:d]
        env = (f[ENV].to_numpy()-em)/es
        h = 50**(-1/(d+4))
        bases.append(dict(group=group,observations=len(f),taxa=len(counts),eligible_taxa_50=int((counts>=50).sum()),
                          input_dimensions=len(members),dimensions=d,retained_variance=float(fraction[:d].sum()),
                          members=members,centre=mu.tolist(),scale=sd.tolist(),loadings=vectors[:,:d].tolist(),bandwidth=h))
        print(json.dumps({k:v for k,v in bases[-1].items() if k not in ['loadings','centre','scale','members']}),flush=True)
        for taxon, idx0 in f.groupby('taxon_name').groups.items():
            idx = np.asarray(list(idx0))
            if len(idx)<50: continue
            r = rng_for(args.seed,group,taxon,'rarefaction')
            for rep in range(args.replicates):
                chosen = r.choice(idx,50,replace=False)
                volume = kde_volume(scores[chosen],h,r,args.draws)
                rare.append(dict(group=group,taxon_name=taxon,replicate=rep,kde_volume=volume,
                                 environment_breadth=float(np.sqrt(np.mean(pdist(env[chosen])**2))),
                                 geographic_breadth_km=geographic_rms(f.latitude.to_numpy()[chosen],f.longitude.to_numpy()[chosen])))
            row = predictive_gain(f.iloc[idx],scores[idx],env[idx],group,taxon,args.seed)
            if row: predictions.append(row)
        pd.DataFrame(rare).to_csv(args.out/'rarefaction.csv',index=False)
    pd.DataFrame(inventory).to_csv(args.out/'support_inventory.csv',index=False)
    (args.out/'shared_bases.json').write_text(json.dumps(bases,indent=2)+'\n')
    pd.DataFrame(predictions).to_csv(args.out/'within_taxon_prediction.csv',index=False)
    rr = pd.DataFrame(rare)
    med = rr.groupby(['group','taxon_name'])[['kde_volume','environment_breadth','geographic_breadth_km']].median().reset_index()
    med.to_csv(args.out/'taxon_breadths.csv',index=False)
    result=[]
    for group,f in med.groupby('group'):
        for predictor in ['environment_breadth','geographic_breadth_km']:
            row=dict(group=group,predictor=predictor,taxa=len(f),status='descriptive_only_fewer_than_15_taxa')
            if len(f)>=15:
                x,y=f[predictor].to_numpy(),f.kde_volume.to_numpy()
                rho=float(spearmanr(x,y).statistic)
                r=rng_for(args.seed,group,predictor,'taxon_bootstrap'); boot=[]
                for _ in range(1000):
                    i=r.integers(len(f),size=len(f)); boot.append(spearmanr(x[i],y[i]).statistic)
                row.update(status='exploratory_no_spatial_or_phylogenetic_correction',rho=rho,
                           bootstrap_low=float(np.nanquantile(boot,.025)),bootstrap_high=float(np.nanquantile(boot,.975)))
            result.append(row)
    pd.DataFrame(result).to_csv(args.out/'breadth_associations.csv',index=False)
    summary=dict(status='pilot_complete',cohort_observations=len(e),cohort_taxa=e.taxon_name.nunique(),environment_complete_observations=len(e_complete),
                 breadth_results=result,predictions_evaluated=len(predictions),
                 prediction_positive_by_group={g:dict(taxa=len(f),positive=int((f.heldout_gain>0).sum()),median_gain=float(f.heldout_gain.median())) for g,f in pd.DataFrame(predictions).groupby('group')})
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
