"""Unified Gates B/C/D — shared native observations, explicit source and scope.

All three stages consume exactly the same pre-declared accepted-taxon × coarse
cell × original observer candidate frame from Gate A. Existing PR144 five-
exposure models are reused for conditional environmental adjustment; this file
does not re-select climate variables or refit 52 tree placements. KDE/TPD and W2
are *actual* numerical densities/distances, never Gaussian-ellipsoid proxies.
"""
from __future__ import annotations

import hashlib
import itertools
import json

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.special import logsumexp

from analysis.v4.native_integration import FEATURES, PAIRS, analyse
from analysis.v4.native_gates_ae import CORE9

MODULES = {
    "presentation": ("orientation_image_vertical_angle",),
    "colour": ("corolla_lab_lightness", "corolla_lab_chroma",
               "corolla_hue_sin", "corolla_hue_cos"),
    "outline": ("capitulum_outline_aspect_ratio", "capitulum_outline_circularity",
                "capitulum_outline_solidity", "capitulum_width_profile_cv"),
}
FIRST_SEED = 20261021


def rng_for(*parts):
    digest = hashlib.sha256("|".join(map(str, parts)).encode("utf-8")).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


def fixed_sample(candidate: pd.DataFrame, n: int = 50, replicate: int = 0):
    if "user_id" not in candidate or "_cell" not in candidate:
        raise ValueError("Gate A observer/cell eligibility required")
    samples = []
    for taxon, f in candidate.groupby("accepted_key", sort=True):
        if len(f) < n:
            raise ValueError("Attempted rarefaction of under-sampled taxon")
        choice = rng_for(FIRST_SEED, taxon, replicate, "A-E-common").choice(
            len(f), size=n, replace=False)
        samples.append(f.iloc[choice])
    if not samples:
        raise ValueError("No shared native taxa passed Gate A")
    return pd.concat(samples, ignore_index=True)


def adjusted_coordinates(sample: pd.DataFrame, models: dict, columns=CORE9):
    """Subtract PR144 fixed effects without re-centering environmental means."""
    selection = models["selection"]
    predictors, scales = selection["selected"], selection["scales"]
    if predictors != ["chelsa_bio01","chelsa_bio04","chelsa_bio12",
                      "chelsa_bio15","chelsa_rsds_mean"]:
        raise ValueError("Active five-variable PR144 model changed")
    components = {k: v for m in models["models"].values()
                  for k,v in m["components"].items()}
    adjusted = sample.loc[:, list(columns)].to_numpy(float).copy()
    for index, name in enumerate(columns):
        v = components[name]
        if v["status"] != "estimated":
            raise ValueError("Environmental coefficient is missing for " + name)
        for predictor in predictors:
            for suffix in ("__between", "__within"):
                term=predictor+suffix
                expected = scales[predictor]["mean"] if suffix == "__between" else 0
                adjusted[:, index] -= (
                    (sample[term].to_numpy(float)-expected)
                    /scales[predictor]["sd"]
                    *v["effects"][term]["raw_units_per_environment_sd"]
                )
    if not np.isfinite(adjusted).all():
        raise ValueError("Non-finite adjusted native phenotype")
    return adjusted


def gate_b(candidate: pd.DataFrame, models: dict, integration: dict,
           *, replicates: int = 200) -> dict:
    """New matched cohort is shared with Gates C and D, not 157654 raw rows."""
    sample = fixed_sample(candidate)
    columns = tuple(CORE9)
    denom = np.array([integration["coordinate_scales"][name] for name in columns], float)
    if np.any(~np.isfinite(denom)) or np.any(denom<=0):
        raise ValueError("Original nine-coordinate scales changed")
    raw = sample.loc[:, list(columns)].to_numpy(float)/denom
    adjusted = adjusted_coordinates(sample,models)/denom
    result = analyse(raw,adjusted,sample.accepted_key.to_numpy(),list(columns),
                     PAIRS,np.random.default_rng(FIRST_SEED),replicates=replicates,equal_n=True)
    if result["status"]!="estimated":
        raise ValueError("Fewer than 20 qualifying accepted taxa: no gate B generalisation")
    cross = []
    module_by_feature = {
        k: next(group for group,names in MODULES.items() if set(fields) <= set(names))
        for k,fields in FEATURES.items()
    }
    for pair in result["pairs"]:
        if module_by_feature[pair["left"]]==module_by_feature[pair["right"]]:
            continue
        x=pair["estimates"]["adjusted/between_minus_within_rv"]
        cross.append({
            "pair":pair["left"]+"/"+pair["right"],
            "adjusted_among_minus_within_rv":float(x["estimate"]),
            "pointwise_interval":[float(x["interval"]["lower"]),
                                  float(x["interval"]["upper"])],
        })
    return {
        "status":"completed_exploratory_matched_native_gate_b",
        "cohort":{"observations":len(sample),"taxa":int(sample.accepted_key.nunique()),
                  "one_of_50_per_taxon_in_all_three_phenotype_modules":True},
        "raw_vs_five_environment_adjusted":True,
        "species_cluster_bootstrap":replicates,
        "B_C_share_first_rarefaction_D_shares_same_candidate_frame":True,
        "cross_module_pairs":cross,
        "full_28_pair_descriptive":result["pairs"],
        "equal_n_diagnostic":result["equal_n_diagnostic"],
        "not_genetic_or_functional_covariance":True,
    }


def _basis(values: np.ndarray, dimension: int):
    values = np.asarray(values,float)
    centre=values.mean(axis=0)
    scale=values.std(axis=0,ddof=0)
    if np.any(~np.isfinite(scale)) or np.any(scale<=1e-12):
        raise ValueError("Degenerate photographic phenotype projection")
    z=(values-centre)/scale
    if dimension == len(scale):
        return centre,scale,np.eye(dimension),1.0
    _, eig, vh=np.linalg.svd(z,full_matrices=False)
    variance=eig**2
    retained=float(variance[:dimension].sum()/variance.sum())
    return centre,scale,vh[:dimension].T,retained


def _kde_log(points: np.ndarray, sample: np.ndarray, h: float):
    d=sample.shape[1]
    return (logsumexp(-cdist(points,sample,"sqeuclidean")/(2*h*h),axis=1)
            -np.log(len(sample))-d*np.log(h*np.sqrt(2*np.pi)))


def _kde_draw(x: np.ndarray,h:float,rng:np.random.Generator,draws:int):
    return x[rng.integers(len(x),size=draws)] + rng.normal(
        size=(draws,x.shape[1]))*h


def _kde_probability_volume(x,h,rng,draws=256):
    a=_kde_draw(x,h,rng,draws)
    threshold=float(np.quantile(_kde_log(a,x,h),.10))
    b=_kde_draw(x,h,rng,draws)
    densities=_kde_log(b,x,h)
    return float(np.mean(np.where(densities>=threshold,np.exp(-densities),0)))


def _density_overlap(x,y,h,draws=128,seed=0):
    rng=np.random.default_rng(seed)
    dx=_kde_draw(x,h,rng,draws);dy=_kde_draw(y,h,rng,draws)
    px=_kde_log(dx,x,h);qx=_kde_log(dx,y,h)
    py=_kde_log(dy,y,h);qy=_kde_log(dy,x,h)
    value=np.r_[
        2*np.exp(np.minimum(px,qx)-np.logaddexp(px,qx)),
        2*np.exp(np.minimum(py,qy)-np.logaddexp(py,qy))
    ]
    return float(np.mean(value))


def wasserstein_parts(x:np.ndarray,y:np.ndarray)->tuple[float,float,float]:
    if x.shape!=y.shape or x.ndim!=2 or len(x)<2 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("W2 needs equal-size finite samples")
    d=cdist(x,y,"sqeuclidean")
    i,j=linear_sum_assignment(d)
    square=float(d[i,j].mean())
    centre=float(np.square(x.mean(axis=0)-y.mean(axis=0)).sum())
    if square<centre-1e-9:
        raise ValueError("Wasserstein moment identity invalid")
    return np.sqrt(square),np.sqrt(centre),np.sqrt(max(0.,square-centre))


def _group_scores(candidate: pd.DataFrame,first:pd.DataFrame,fields,dimensions):
    # Shared PCA fit on equal numbers per taxon, never on environment response.
    fields=list(fields)
    v=first[fields].to_numpy(float)
    centre=v.mean(axis=0)
    scale=v.std(axis=0,ddof=0)
    if "corolla_hue_sin" in fields:
        ix,iy=[fields.index(q) for q in ("corolla_hue_sin","corolla_hue_cos")]
        shared=np.sqrt((scale[ix]**2+scale[iy]**2)/2)
        scale[ix]=scale[iy]=shared  # one circular colour unit, no double trait
    if not np.isfinite(scale).all() or np.any(scale<=1e-12):
        raise ValueError("Degenerate shared colour/outline coordinate scale")
    z=(v-centre)/scale
    if dimensions>=len(fields):
        loadings=np.eye(len(fields)); retained=1.
    else:
        _,sv,vh=np.linalg.svd(z,full_matrices=False)
        loadings=vh[:dimensions].T
        retained=float((sv[:dimensions]**2).sum()/(sv**2).sum())
    raw=candidate[fields].to_numpy(float)
    return ((raw-centre)/scale)@loadings,retained

def gate_c(candidate:pd.DataFrame,*,replicates:int=12,n:int=50,
           kde_draws:int=128,max_taxa:int=40)->dict:
    """True cross-taxon KDE overlap, 90% density-content volume and exact W2.
    
    This is a screening run on independent human×cell strata, not a full
    ecological/functional hypervolume or multiple-test corrected result.
    """
    if replicates<2 or n<10 or kde_draws<32:
        raise ValueError("Insufficient numerical KDE screening settings")
    taxa=sorted(candidate.accepted_key.unique())
    if len(taxa)<5 or len(taxa)>max_taxa:
        return {"status":"insufficient_or_excessive_taxon_support",
                "available_taxa":len(taxa),"required_minimum":5,"maximum":max_taxa}
    first=fixed_sample(candidate,n=n,replicate=0)
    groups={}
    for name,fields in MODULES.items():
        d=1 if name=="presentation" else 3
        scores,retained=_group_scores(candidate,first,fields,d)
        data = {taxon: scores[np.flatnonzero(candidate.accepted_key.to_numpy()==taxon)]
                for taxon in taxa}
        # Taxon-specific finite-sample baseline requires two disjoint sets of 50.
        null={}
        for taxon,values in data.items():
            if len(values)<2*n:
                continue
            total=[];shape=[]
            for rep in range(min(20,replicates)):
                r=rng_for("same",name,taxon,rep)
                ids=r.choice(len(values),size=2*n,replace=False)
                w,_,centred=wasserstein_parts(values[ids[:n]],values[ids[n:]])
                total.append(w);shape.append(centred)
            null[taxon]={"full":float(np.quantile(total,.95)),
                         "centred":float(np.quantile(shape,.95))}
        all_overlap,all_w2,all_centres,all_centreless,all_volumes=[],[],[],[],[]
        test_all,test_centred,referenced_pairs=0,0,0
        n_pairs=len(taxa)*(len(taxa)-1)//2
        h=n**(-1/(d+4))
        for rep in range(replicates):
            s={};vol=[]
            for taxon in taxa:
                r=rng_for("sample",name,taxon,rep)
                x=data[taxon][r.choice(len(data[taxon]),size=n,replace=False)]
                s[taxon]=x
                vol.append(_kde_probability_volume(x,h,rng_for("vol",name,taxon,rep),
                                                    draws=kde_draws))
            all_volumes.append(float(np.median(vol)))
            for i,a in enumerate(taxa):
                for b in taxa[i+1:]:
                    x,y=s[a],s[b]
                    ov=_density_overlap(x,y,h,draws=kde_draws,
                           seed=int(rng_for("density",name,a,b,rep).integers(2**31)))
                    w,centre,shape=wasserstein_parts(x,y)
                    all_overlap.append(ov);all_w2.append(w)
                    all_centres.append(centre);all_centreless.append(shape)
                    if rep==0 and a in null and b in null:
                        referenced_pairs+=1
                        test_all+=int(w>max(null[a]["full"],null[b]["full"]))
                        test_centred+=int(shape>max(null[a]["centred"],null[b]["centred"]))
        if not all_overlap:
            raise ValueError("No KDE comparisons for sampled taxa")
        groups[name]={
            "dimensions":d,
            "retained_PCA_variance":retained,
            "taxa":len(taxa),"taxon_pairs":n_pairs,
            "replicates":replicates,
            "median_density_overlap":float(np.median(all_overlap)),
            "median_empirical_wasserstein_w2":float(np.median(all_w2)),
            "median_centroid_distance":float(np.median(all_centres)),
            "median_centre_removed_w2":float(np.median(all_centreless)),
            "median_90pct_KDE_probability_volume":float(np.median(all_volumes)),
            "kde_shared_bandwidth":h,
            "same_taxon_baseline_taxa":len(null),
            "taxon_pairs_both_with_same_taxon_null":referenced_pairs,
            "reference_exceedances_raw_W2":test_all,
            "reference_exceedances_centre_removed_W2":test_centred,
            "baseline_warning":"0 reference pairs means finite-sample distinctiveness not evaluable",
        }
    return {
        "status":"completed_native_density_screen_not_functional_hypervolume",
        "shared_taxon_sample_every_module":True,
        "observer_cell_balanced_source":True,
        "per_taxon_rarefaction":n,"replicates":replicates,
        "modules":groups,
        "scope_limits":[
            "Uncalibrated image-relative phenotypes, not anatomical physical dimensions.",
            "Cell and observer strata reduce repeated-source bias, not proven individual dependence.",
            "PCA projections share first rarefaction's measurement-only coordinate system.",
            "Hue sine and cosine are coordinates of one circular colour axis.",
            "KDE volumes across differing dimensions or retained fractions are not commensurable.",
            "Species pairs are correlated and rarefaction ranges are not P values.",
            "No significant evolutionary or functional equivalence inferred.",
        ]
    }


def gate_d_same_sample(candidate:pd.DataFrame,full_native:pd.DataFrame,models:dict,
                       *,n:int=50,minimum_side:int=10)->dict:
    """Within-taxon BIO12/RSDS morphospace shift in exactly C's taxon set.

    Gradients are descriptive climate-space contrasts; the primary conditioned
    climate fits, spatial models and phylogeny are separately kept in PR144.
    """
    first=fixed_sample(candidate,n=n)
    results={}
    for var in ("chelsa_bio12","chelsa_rsds_mean"):
        # Source-wide medians are computed WITHOUT trait eligibility selection.
        med=full_native.groupby("accepted_key")[var].median()
        per_module={}
        for name,fields in MODULES.items():
            d=1 if name=="presentation" else 3
            x,retained=_group_scores(candidate,first,fields,d)
            marked=candidate.copy()
            marked["_basis_index"]=np.arange(len(marked))
            contributions=[]
            for taxon,group in marked.groupby("accepted_key",sort=True):
                eligible=group[var].notna()
                v=group.loc[eligible]
                low=v.loc[v[var]<=med.loc[taxon]]
                high=v.loc[v[var]>med.loc[taxon]]
                if min(len(low),len(high))<minimum_side:
                    continue
                if min(low._cell.nunique(),high._cell.nunique())<3:
                    continue
                nn=min(len(low),len(high),n//2)
                a=rng_for("regime",var,name,taxon,0)
                li=a.choice(low["_basis_index"].to_numpy(int),size=nn,replace=False)
                hi=a.choice(high["_basis_index"].to_numpy(int),size=nn,replace=False)
                w,centre,shape=wasserstein_parts(x[li],x[hi])
                contributions.append((w,centre,shape))
            per_module[name]={
                "taxa_with_both_regimes_and_three_cells_per_side":len(contributions),
                "retained_measurement_space_variance":retained,
                "median_low_high_empirical_w2":float(np.median([v[0] for v in contributions]))
                   if contributions else None,
                "median_low_high_centre_distance":float(np.median([v[1] for v in contributions]))
                   if contributions else None,
                "median_low_high_shape_w2":float(np.median([v[2] for v in contributions]))
                   if contributions else None,
            }
        results[var]=per_module
    return {
        "status":"descriptive_native_environment_density_structuring_screen",
        "model_owner":"PR144 unchanged five-variable within/between adjusted fit",
        "environment_split":"Within each accepted taxon, source-native unconditional median (not phenotype-selected)",
        "same_candidate_observer_cell_frame_as_gates_B_C_D":True,
        "source_quality":"All native records contribute to environment median when finite",
        "predictors":results,
        "scope":"Not a controlled environmental effect, phylogenetic result, or spatially external prediction. Geographic climate sorting and observer selection remain possible.",
    }
