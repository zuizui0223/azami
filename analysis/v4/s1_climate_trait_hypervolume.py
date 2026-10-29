"""S1 exploratory ecology: climate occupancy breadth versus floral morphospace breadth.

Predeclared post-freeze third question, independent of S1/S2 regression families.
Uses source-verified private native S1 photographs only. Three functional-candidate
photographic axes: camera-relative capitulum angle, CIELAB chroma, 2D outline
aspect. Three realized climate occupancy axes: BIO12, BIO01 and mean RSDS.

Equal geographic-cell medians rather than photo abundance define each taxon's
points. A rarefied regularised Gaussian-ellipsoid log-volume is a *proxy* for
occupied hypervolume, NOT a KDE Hutchinsonian niche, fitness, flexibility,
plasticity or function measured directly. Evolutionary convergence/radiation
requires source-resolved independent phylogeny, not tested here.

No source observation, taxon, uploader or geographic identifiers leave Actions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
import statsmodels.api as sm

CLIMATE = ("chelsa_bio12", "chelsa_bio01", "chelsa_rsds_mean")
PHENOTYPE = ("orientation_image_vertical_angle", "corolla_lab_chroma",
             "capitulum_outline_aspect_ratio")
AXES = CLIMATE + PHENOTYPE
SOURCE_N = 313350
SEED = 20261029
MIN_PHOTOS = 40
MIN_CELLS = 8
RAREFIED_CELLS = 8
N_DRAWS = 96
N_PERMUTATIONS = 399
RIDGE = 0.025
MIN_TAXA = 25
MIN_CROSS_REGION_PAIRS = 50
REGION_NAMES = ("western", "central", "eastern")


def _required_cols():
    return ["obs_id", "accepted_key", "latitude", "longitude", *AXES,
            *(p + "__eligible" for p in PHENOTYPE)]


def load_source(table_path: Path, receipt_path: Path):
    receipt=json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    source=Path(table_path)
    sha=hashlib.sha256(source.read_bytes()).hexdigest()
    if (receipt.get("prepared_table_sha256")!=sha or
        receipt.get("native_observations")!=SOURCE_N or
        receipt.get("environment_decoding_status")!="physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("SOURCE_ARCHIVE_CHECKSUM_OR_DECODE_NOT_VALIDATED")
    d=pd.read_csv(source, usecols=_required_cols(),
                  dtype={"obs_id":str},low_memory=False)
    if len(d)!=SOURCE_N or d.obs_id.isna().any() or d.obs_id.duplicated().any():
        raise ValueError("ORIGINAL_NATIVE_OBSERVATION_GRAIN_DRIFT")
    mask=np.isfinite(d[list(AXES)+["latitude","longitude"]].to_numpy(float)).all(axis=1)
    for name in PHENOTYPE:
        mask &= d[name+"__eligible"].fillna(False).eq(True).to_numpy(bool)
    mask &= d[PHENOTYPE[0]].between(0,180).to_numpy(bool)
    d=d.loc[mask,["accepted_key","latitude","longitude",*AXES]].copy()
    if (len(d)<50000 or d.accepted_key.nunique()<75 or
        not d.latitude.between(-90,90).all() or
        not d.longitude.between(-180,180).all()):
        raise ValueError("INSUFFICIENT_MATCHED_HEAD_CLIMATE_SOURCE")
    return d,{"decoded_sha256":sha,"source_observations":SOURCE_N,
              "quality_eligible_matched_observations":len(d),
              "matched_source_taxa":int(d.accepted_key.nunique())}


def geo_aggregate(d: pd.DataFrame, degrees: int):
    if degrees not in (1,5):
        raise ValueError("CELL_SIZE_NOT_PREDECLARED")
    work=d.copy()
    work["_latcell"]=np.floor((work.latitude+90)/degrees).astype(int)
    work["_loncell"]=np.floor((work.longitude+180)/degrees).astype(int)
    grouped=work.groupby(["accepted_key","_latcell","_loncell"],sort=True)
    cells=grouped.agg(
        **{name:(name,"median") for name in AXES},
        n_photos=(PHENOTYPE[0],"size"),
        latitude=("latitude","median"),
        longitude=("longitude","median"),
    ).reset_index()
    counts=cells.groupby("accepted_key",sort=True).agg(
        n_cells=(PHENOTYPE[0],"size"),
        n_photos=("n_photos","sum"))
    qualified=counts.index[(counts.n_cells>=MIN_CELLS)&
                            (counts.n_photos>=MIN_PHOTOS)]
    cells=cells.loc[cells.accepted_key.isin(qualified)].copy()
    return cells,{
        "degrees":degrees,
        "source_photos":len(work),
        "source_taxa":int(work.accepted_key.nunique()),
        "source_taxa_with_cell_and_photo_support":len(qualified),
        "eligible_taxon_cell_pairs":len(cells),
        "n_cells_gate":MIN_CELLS,
        "n_photos_gate":MIN_PHOTOS,
    }


def standardise_cells(cells: pd.DataFrame):
    """Across-taxonomy pooled z scales with equal taxon total weight."""
    x=cells[list(AXES)].to_numpy(float)
    counts=cells.groupby("accepted_key").size()
    w=cells.accepted_key.map(1/counts).to_numpy(float)
    w/=w.sum()
    mean=np.average(x,axis=0,weights=w)
    sd=np.sqrt(np.average((x-mean)**2,axis=0,weights=w))
    if not np.isfinite(sd).all() or np.any(sd<=1e-10):
        raise ValueError("ZERO_OR_NAN_CLIMATE_OR_TRAIT_GLOBAL_SCALING")
    result=cells.copy()
    for idx,name in enumerate(AXES):
        result["_z"+str(idx)]=(x[:,idx]-mean[idx])/sd[idx]
    return result,{"global_mean":dict(zip(AXES,[float(v) for v in mean])),
                   "global_sd":dict(zip(AXES,[float(v) for v in sd])),
                   "equal_taxon_total_weight":True}


def ellipsoid_log_volume(samples: np.ndarray) -> float:
    """Regularised ellipsoidal log-volume; constant ellipsoid factor omitted."""
    x=np.asarray(samples,float)
    if x.ndim!=2 or x.shape[1]!=3 or len(x)<RAREFIED_CELLS:
        raise ValueError("THREE_AXIS_HYPERVOLUME_SUPPORT_REQUIRED")
    cov=np.cov(x,rowvar=False,ddof=1)
    sign,logdet=np.linalg.slogdet(cov+RIDGE*np.eye(3))
    if sign<=0 or not np.isfinite(logdet):
        raise ValueError("REGULARISED_VOLUME_NOT_POSITIVE")
    return float(logdet/2)


def center_spread_km(frame: pd.DataFrame) -> float:
    lat,lon=np.deg2rad(frame.latitude.to_numpy(float)),np.deg2rad(frame.longitude.to_numpy(float))
    xyz=np.column_stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    unit=xyz.mean(axis=0)
    unit/=max(np.linalg.norm(unit),1e-12)
    angle=np.arccos(np.clip(xyz@unit,-1,1))
    return float(6371*np.sqrt(np.mean(angle**2)))


def draw_taxon_scores(cells: pd.DataFrame, *, draws=N_DRAWS):
    """Equal cell opportunity, same sampled cells for climate and phenotype."""
    rng=np.random.default_rng(SEED)
    rows=[]
    for taxon,part in cells.groupby("accepted_key",sort=True):
        if len(part)<RAREFIED_CELLS:
            continue
        x=part[["_z0","_z1","_z2"]].to_numpy(float)
        y=part[["_z3","_z4","_z5"]].to_numpy(float)
        e=[];p=[]
        for k in range(draws):
            indices=rng.choice(len(part),RAREFIED_CELLS,replace=False)
            e.append(ellipsoid_log_volume(x[indices]))
            p.append(ellipsoid_log_volume(y[indices]))
        dominant=int(np.argmax([
            int((part.longitude<-30).sum()),
            int(((part.longitude>=-30)&(part.longitude<60)).sum()),
            int((part.longitude>=60).sum()),
        ]))
        rows.append({
            "_taxon":str(taxon),
            "climate_log_volume":float(np.mean(e)),
            "trait_log_volume":float(np.mean(p)),
            "climate_log_volume_sd":float(np.std(e)),
            "trait_log_volume_sd":float(np.std(p)),
            "n_cells":len(part),
            "n_photos":int(part.n_photos.sum()),
            "geo_spread_km":center_spread_km(part),
            "sector":dominant,
            **{f"clim_centroid_{i}":float(v) for i,v in enumerate(x.mean(axis=0))},
            **{f"trait_centroid_{i}":float(v) for i,v in enumerate(y.mean(axis=0))},
        })
    return pd.DataFrame(rows)


def _design(scores):
    """Keep climate breadth first; exclude constant/collinear nuisance covariates.

    This is not predictor selection for a discovery claim. The three declared
    sampling/geography controls are attempted in order and skipped only when
    algebraically non-identifiable on the actual taxon eligibility cohort.
    """
    climate=scores.climate_log_volume.to_numpy(float)
    base=np.column_stack([np.ones(len(scores)), climate])
    if not np.isfinite(base).all() or np.linalg.matrix_rank(base)!=2:
        raise ValueError("CLIMATE_BREADTH_FOCAL_TERM_NOT_IDENTIFIABLE")
    controls=[
        ("log geographic spread",np.log1p(scores.geo_spread_km.to_numpy(float))),
        ("log observed photos",np.log1p(scores.n_photos.to_numpy(float))),
        ("log occupied geographic cells",np.log1p(scores.n_cells.to_numpy(float))),
    ]
    group=scores.sector.to_numpy(int)
    controls.extend((f"sector_{j}",(group==j).astype(float)) for j in (0,1))
    kept=[];dropped=[]
    x=base
    for name,col in controls:
        if not np.isfinite(col).all():
            raise ValueError("NONFINITE_COVARIATE")
        candidate=np.column_stack([x,col])
        if (np.linalg.matrix_rank(candidate)==x.shape[1]+1 and
            np.linalg.cond(candidate)<1e9):
            x=candidate;kept.append(name)
        else:
            dropped.append(name)
    if len(scores)<=x.shape[1]+5:
        raise ValueError("HYPERVOLUME_REGRESSION_NOT_IDENTIFIABLE")
    return x,kept,dropped


def breadth_regression(scores: pd.DataFrame, permutations=N_PERMUTATIONS):
    if len(scores)<MIN_TAXA:
        return {"status":"insufficient_taxonomic_support","n_taxa":len(scores),"required":MIN_TAXA}
    x,active,dropped=_design(scores)
    null=x[:,[0,*range(2,x.shape[1])]]
    y=scores.trait_log_volume.to_numpy(float)
    fit=sm.OLS(y,x).fit(cov_type="HC3")
    beta=float(fit.params[1]);ci=[float(q) for q in fit.conf_int(alpha=.05)[1]]
    # Freedman-Lane residual permutations within the response-blind geo sectors.
    nullfit=sm.OLS(y,null).fit()
    prediction=nullfit.fittedvalues
    residual=nullfit.resid
    rng=np.random.default_rng(SEED+1)
    sector=scores.sector.to_numpy(int)
    permutations_beta=[]
    for k in range(permutations):
        u=residual.copy()
        for category in np.unique(sector):
            idx=np.where(sector==category)[0]
            u[idx]=rng.permutation(residual[idx])
        permutations_beta.append(float(np.linalg.lstsq(x,prediction+u,rcond=None)[0][1]))
    p=(1+sum(abs(q)>=abs(beta) for q in permutations_beta))/(permutations+1)
    return {"status":"estimated","n_taxa":len(scores),
            "trait_logvolume_slope_per_climate_logvolume":beta,
            "HC3_95CI":ci,
            "HC3_p_pointwise":float(fit.pvalues[1]),
            "sector_conditioned_residual_permutation_p_exploratory":p,
            "n_permutations":permutations,
            "adjusted_R2":float(fit.rsquared_adj),
            "controls_retained":active,
            "controls_omitted_algebraically":dropped,
            "not_causal_plasticity_or_adaptation":True}


def cross_region_convergence(scores: pd.DataFrame, permutations=N_PERMUTATIONS):
    if len(scores)<MIN_TAXA:
        return {"status":"insufficient_taxonomic_support","n_taxa":len(scores)}
    climate=scores[[f"clim_centroid_{i}" for i in range(3)]].to_numpy(float)
    phenotype=scores[[f"trait_centroid_{i}" for i in range(3)]].to_numpy(float)
    group=scores.sector.to_numpy(int)
    pairs=np.triu(group[:,None]!=group[None,:],1)
    n=int(pairs.sum())
    if n<MIN_CROSS_REGION_PAIRS or len(set(group))<2:
        return {"status":"insufficient_cross_region_pairs","pairs":n}
    x=squareform(pdist(climate))[pairs]
    y=squareform(pdist(phenotype))[pairs]
    observed=float(spearmanr(x,y).statistic)
    if not np.isfinite(observed):
        return {"status":"degenerate_centroid_distances","pairs":n}
    rng=np.random.default_rng(SEED+2)
    permuted=[]
    for k in range(permutations):
        z=phenotype.copy()
        for g in np.unique(group):
            ix=np.where(group==g)[0]
            z[ix]=phenotype[rng.permutation(ix)]
        dist=squareform(pdist(z))[pairs]
        permuted.append(float(spearmanr(x,dist).statistic))
    p=(1+sum(abs(v)>=abs(observed) for v in permuted))/(permutations+1)
    return {"status":"estimated",
            "cross_sector_taxon_pairs":n,
            "climate_distance_vs_trait_distance_spearman":observed,
            "within_sector_taxon_label_permutation_p_exploratory":p,
            "n_permutations":permutations,
            "cross_sector_centroids_not_shared_climate_ancestry_control":True,
            "phylogenetic_evolutionary_convergence_NOT_tested":True}


def run(source, receipt, out):
    d,provenance=load_source(source,receipt)
    results={}
    for resolution in (1,5):
        cells,coverage=geo_aggregate(d,resolution)
        if coverage["source_taxa_with_cell_and_photo_support"]<MIN_TAXA:
            results[f"{resolution}degree"]={"status":"insufficient_taxa",**coverage}
            continue
        cells,scale=standardise_cells(cells)
        scores=draw_taxon_scores(cells)
        if scores._taxon.duplicated().any():
            raise ValueError("DUPLICATE_TAXON_HYPERVOLUME")
        breadth=breadth_regression(scores)
        convergence=cross_region_convergence(scores)
        results[f"{resolution}degree"]={
            "status":"source_matched_exploratory",
            "source_population":coverage,
            "scaling":scale,
            "n_taxa_eligible":len(scores),
            "sector_taxa_counts":{name:int((scores.sector==j).sum())
                                  for j,name in enumerate(REGION_NAMES)},
            "climate_versus_multitrait_breadth":breadth,
            "cross_region_centroid_similarity":convergence,
            "sampling_design":{
                "equal_occupied_geographical_cells":True,
                "each_taxon_rarefied_to_cells":RAREFIED_CELLS,
                "repeated_without_replacement_draws":N_DRAWS,
                "three_climate_axes":list(CLIMATE),
                "three_photographic_morphology_axes":list(PHENOTYPE),
                "covariance_ellipsoid_regulariser":RIDGE,
                "not_KDE_true_hutchinsonian_hypervolume":True,
                "observer_lighting_and_photo_angle_not_fully_removed":True,
                "climate_occupancy_not_fundamental_niche":True,
                "phylogenetic_branch_distances_unavailable":True,
            },
        }
    result={"schema":"azami_S1_climate_trait_morphospace_breadth_20261029_v1",
            "original_decoded_source":provenance,
            "predeclared_primary":"1degree",
            "secondary_resolution":"5degree",
            "analysis":results,
            "inference":{
                "exploratory_post_freeze_third_ecological_question":True,
                "separate_from_original_S1_80_and_S2_128_Holm_families":True,
                "single_model_family_2_biological_contrasts_not_independent_evolutionary_replication":True,
                "no_evolutionary_radiation_or_convergence_without_source_matched_phylogeny":True,
                "no_taxon_ids_source_coords_or_photo_identifiers_exported":True,
            }}
    dst=Path(out)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    # console and artifact contain taxon aggregate counts only, never scored taxa.
    print(json.dumps({"status":"S1_HYPERVOLUME_ECOLOGICAL_ANALYSIS",
                      "source":provenance,
                      "analyses":{k:{"taxa":v.get("n_taxa_eligible",0),
                                    "status":v["status"],
                                    "breadth":v.get("climate_versus_multitrait_breadth"),
                                    "similarity":v.get("cross_region_centroid_similarity")}
                                  for k,v in results.items()}},allow_nan=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.source,args.receipt,args.out)


if __name__=="__main__":
    main()
