"""H3c: geography-separated climatic analogues and 3D photo-trait alignment.

Does a *given pair* of source taxa from different geographical sectors look more
similar at environmental-analogue cell pairs than when both taxa's observed
cells are paired without environmental matching?  The null exchanges complete
3D floral phenotype vectors across cells *within each source taxon*, sharing
the same permutation across all cross-taxon pairs.  This is ecological matching,
not phylogenetic convergence, selection, response plasticity or fitness.

Source: the exact decoded original confidential S1 photographs.  Export ONLY
aggregate statistics; no taxon or observer identifiers, source coordinates or
individual photographic observations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from analysis.v4.s1_climate_trait_hypervolume import (
    AXES, CLIMATE, PHENOTYPE, SOURCE_N,
    load_source, geo_aggregate, standardise_cells
)

SEED = 20261029
RAREFIED_CELLS = 8
N_DRAWS = 24
N_PERMUTATIONS = 199
CLIMATE_ANALOG_Z = 0.9
KDE_BANDWIDTH = 0.75  # Fixed z-scaled 3D floral-space density bandwidth
MIN_ANALOG_CELL_PAIRS = 3
MIN_DISTINCT_CELLS_PER_SIDE = 2
MIN_ELIGIBLE_TAXON_PAIRS = 40
MIN_SUCCESSFUL_DRAWS = 12
MIN_SECTOR_PAIR_TYPES = 2
RESOLUTIONS = (1, 5)


def geographic_sector(longitude):
    lon=np.asarray(longitude, float)
    if not np.isfinite(lon).all() or np.any(abs(lon)>180):
        raise ValueError("INVALID_SOURCE_GEOGRAPHY")
    return np.select((lon < -30, lon < 60), (0, 1), default=2).astype(int)


def eligible_taxa(cells: pd.DataFrame):
    """Make one region assignment per taxon; never publish source taxon IDs."""
    taxa=[]
    for taxon,part in cells.groupby("accepted_key",sort=True):
        if len(part)<RAREFIED_CELLS:
            continue
        sector=geographic_sector(part.longitude.to_numpy(float))
        counts=np.bincount(sector,minlength=3)
        dominant=int(np.argmax(counts))
        taxa.append({
            "taxon":str(taxon),  # private in runner; never exported
            "sector":dominant,
            "x":part[["_z0","_z1","_z2"]].to_numpy(float),
            "y":part[["_z3","_z4","_z5"]].to_numpy(float),
        })
    return taxa


def _sample_taxa(taxa, rng):
    x=[];y=[];sec=[]
    for t in taxa:
        pick=rng.choice(len(t["x"]),RAREFIED_CELLS,replace=False)
        x.append(t["x"][pick]);y.append(t["y"][pick]);sec.append(t["sector"])
    return np.array(x),np.array(y),np.array(sec,int)


def matched_pairs(x: np.ndarray, y: np.ndarray, sector: np.ndarray):
    """Return pairwise original morphology distances and climate-only masks.

    Pair cell counts and supported sector pairs are all source blind to y.
    """
    if (x.shape!=y.shape or x.ndim!=3 or x.shape[1:]!=(RAREFIED_CELLS,3)
        or sector.shape!=(x.shape[0],)):
        raise ValueError("INVALID_MULTIVARIATE_SOURCE")
    group_i=[];group_j=[];distances=[];masks=[];counts=[];types=[]
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            if sector[i]==sector[j]:
                continue
            d_env=np.linalg.norm(x[i,:,None,:]-x[j,None,:,:],axis=2)
            matched=(d_env<=CLIMATE_ANALOG_Z)
            if (int(matched.sum())<MIN_ANALOG_CELL_PAIRS
                or np.flatnonzero(matched.any(axis=1)).size<MIN_DISTINCT_CELLS_PER_SIDE
                or np.flatnonzero(matched.any(axis=0)).size<MIN_DISTINCT_CELLS_PER_SIDE):
                continue
            d_trait=np.linalg.norm(y[i,:,None,:]-y[j,None,:,:],axis=2)
            group_i.append(i);group_j.append(j)
            distances.append(d_trait)
            masks.append(matched)
            counts.append(int(matched.sum()))
            types.append(tuple(sorted((int(sector[i]),int(sector[j])))))
    if not group_i:
        return None,{"status":"no_analogous_pair_support","n_pairs":0,
                     "n_sector_pair_types":0}
    unique=len(set(types))
    meta={"status":"qualified" if
          (len(group_i)>=MIN_ELIGIBLE_TAXON_PAIRS and unique>=MIN_SECTOR_PAIR_TYPES)
          else "insufficient_analogue_pairs",
          "n_pairs":len(group_i),
          "n_sector_pair_types":unique,
          "n_analog_cell_pairs":sum(counts),
          "mean_analog_cell_pairs_per_taxon_pair":float(np.mean(counts)),
          "n_taxa_represented":len(set(group_i)|set(group_j))}
    return {
        "i":np.array(group_i,int),
        "j":np.array(group_j,int),
        "d":np.stack(distances),
        "mask":np.stack(masks),
        "analog_n":np.array(counts,float),
    }, meta


def observed_and_permutations(pair, n_taxa: int, rng, reps=N_PERMUTATIONS):
    """One vote per source taxon pair; permute whole phenotype vectors per taxon."""
    if pair is None:
        raise ValueError("SOURCE_ANALOG_PAIR_REQUIRED")
    dist=pair["d"];mask=pair["mask"]
    if dist.shape!=mask.shape or dist.shape[1:]!=(RAREFIED_CELLS,RAREFIED_CELLS):
        raise ValueError("PAIR_METRIC_SHAPE_DRIFT")
    i,j=pair["i"],pair["j"]
    n=np.arange(len(dist))[:,None,None]
    all_means=dist.mean(axis=(1,2))
    analog_means=(dist*mask).sum(axis=(1,2))/pair["analog_n"]
    observed=float(np.mean(all_means-analog_means))
    null=np.empty(reps,float)
    for k in range(reps):
        perms=np.argsort(rng.random((n_taxa,RAREFIED_CELLS)),axis=1)
        rearranged=dist[n,perms[i][:,:,None],perms[j][:,None,:]]
        shuffled=(rearranged*mask).sum(axis=(1,2))/pair["analog_n"]
        null[k]=float(np.mean(all_means-shuffled))
    if not np.isfinite(null).all():
        raise ValueError("NONFINITE_ANALOGUE_NULL")
    return observed,null


def kde_overlap_and_permutations(pair, y: np.ndarray, rng,
                                 reps=N_PERMUTATIONS,
                                 bandwidth=KDE_BANDWIDTH):
    """Cross-taxon overlap of *conditional Gaussian-KDE distributions*.

    The overlap is the normalised L2 inner product between two KDE densities:
      <f_i,f_j>/sqrt(<f_i,f_i><f_j,f_j>), mathematically in [0,1].
    Each KDE is based on the subset of a taxon's 8 occupied cells participating
    in at least one cross-taxon climatic analogue, not on an estimated bounded
    or multi-modal Hutchinsonian niche set.  Full 8-cell KDEs supply the paired
    comparator. The null permutes intact 3D phenotype vectors within taxon.
    """
    if y.ndim!=3 or y.shape[1:]!=(RAREFIED_CELLS,3):
        raise ValueError("ORIGINAL_3D_FLORAL_CELLS_REQUIRED")
    pi,pj=pair["i"],pair["j"]
    if len(pi)==0 or bandwidth<=0:
        raise ValueError("INVALID_KERNEL_OVERLAP_SOURCE")
    kernel_self=np.exp(
        -np.sum((y[:,:,None,:]-y[:,None,:,:])**2,axis=-1)/
        (4*bandwidth**2))
    kernel_cross=np.exp(-(pair["d"]**2)/(4*bandwidth**2))
    full=kernel_cross.mean(axis=(1,2))/np.sqrt(
         kernel_self[pi].mean(axis=(1,2))*kernel_self[pj].mean(axis=(1,2)))
    row=pair["mask"].any(axis=2).astype(float)
    col=pair["mask"].any(axis=1).astype(float)
    na,nb=row.sum(axis=1),col.sum(axis=1)
    if (na<MIN_DISTINCT_CELLS_PER_SIDE).any() or (nb<MIN_DISTINCT_CELLS_PER_SIDE).any():
        raise ValueError("KDE_CONDITIONAL_SUPPORT_MISSING")
    select_cross=row[:,:,None]*col[:,None,:]
    select_i=row[:,:,None]*row[:,None,:]
    select_j=col[:,:,None]*col[:,None,:]
    def cos(kcross,ki,kj):
        dot=(kcross*select_cross).sum(axis=(1,2))/(na*nb)
        ni=(ki*select_i).sum(axis=(1,2))/(na**2)
        nj=(kj*select_j).sum(axis=(1,2))/(nb**2)
        answer=dot/np.sqrt(ni*nj)
        if not np.isfinite(answer).all() or np.any(answer<-.000001) or np.any(answer>1.000001):
            raise ValueError("KDE_L2_OVERLAP_OUTSIDE_UNIT_INTERVAL")
        return np.clip(answer,0,1)
    conditional=cos(kernel_cross,kernel_self[pi],kernel_self[pj])
    if not np.isfinite(full).all() or np.any(full>1.000001):
        raise ValueError("KDE_FULL_OVERLAP_INVALID")
    obs=float((conditional-full).mean())
    target=np.arange(len(pi))[:,None,None]
    permuted=np.empty(reps,float)
    for k in range(reps):
        rearrange=np.argsort(rng.random((len(y),RAREFIED_CELLS)),axis=1)
        p_i=rearrange[pi]
        p_j=rearrange[pj]
        kc=kernel_cross[target,p_i[:,:,None],p_j[:,None,:]]
        ki=kernel_self[pi][target,p_i[:,:,None],p_i[:,None,:]]
        kj=kernel_self[pj][target,p_j[:,:,None],p_j[:,None,:]]
        permuted[k]=float((cos(kc,ki,kj)-full).mean())
    return {
        "conditional_mean":float(conditional.mean()),
        "unconditional_mean":float(full.mean()),
        "conditional_minus_unconditional_mean":obs,
    },permuted


def analogue_analysis(cells: pd.DataFrame, *, draws=N_DRAWS,
                      permutations=N_PERMUTATIONS, seed=SEED):
    """Draws are repeated subsamples of the *same* taxa, not independent tests."""
    taxa=eligible_taxa(cells)
    sector_counts={f"sector_{k}":sum(t["sector"]==k for t in taxa)
                   for k in range(3)}
    if len(taxa)<25 or min(np.unique([t["sector"] for t in taxa]).size,3)<2:
        return {"status":"insufficient_taxa_or_sectors",
                "n_taxa":len(taxa),"sector_taxa_counts":sector_counts}
    rng=np.random.default_rng(seed)
    measurements=[];nulls=[];supports=[];n_qualified=0
    overlaps=[];overlap_nulls=[]
    for rep in range(draws):
        x,y,sector=_sample_taxa(taxa,rng)
        pair,meta=matched_pairs(x,y,sector)
        supports.append(meta)
        if meta["status"]!="qualified":
            continue
        estimate,null=observed_and_permutations(pair,len(taxa),rng,permutations)
        overlap,null_overlap=kde_overlap_and_permutations(
            pair,y,rng,permutations)
        measurements.append(estimate)
        nulls.append(null)
        overlaps.append(overlap)
        overlap_nulls.append(null_overlap)
        n_qualified+=1
    if n_qualified<min(MIN_SUCCESSFUL_DRAWS,draws):
        return {
            "status":"insufficient_overlapping_climatic_analogues",
            "n_taxa":len(taxa),
            "sector_taxa_counts":sector_counts,
            "successful_rarefaction_draws":n_qualified,
            "requested_rarefaction_draws":draws,
            "mean_cross_sector_pairs":float(np.mean([m["n_pairs"] for m in supports])),
            "min_required_pairs_per_draw":MIN_ELIGIBLE_TAXON_PAIRS,
            "min_sector_pair_types":MIN_SECTOR_PAIR_TYPES,
            "no_effect_or_p_value_estimated":True,
        }
    estimates=np.array(measurements)
    nulls=np.array(nulls)
    observed=float(estimates.mean())
    permuted=nulls.mean(axis=0)
    expected=float(permuted.mean())
    p=float((1+np.count_nonzero(np.abs(permuted-expected)>=abs(observed-expected)))/
            (permutations+1))
    source_pairs=[q["n_pairs"] for q in supports if q["status"]=="qualified"]
    overall_overlap=np.mean([q["conditional_minus_unconditional_mean"]
                             for q in overlaps])
    null_overlap=np.mean(np.array(overlap_nulls),axis=0)
    null_overlap_mean=float(null_overlap.mean())
    p_overlap=float((1+np.count_nonzero(
        np.abs(null_overlap-null_overlap_mean)>=
        abs(overall_overlap-null_overlap_mean)))/(permutations+1))
    p_adjusted=multipletests([p,p_overlap],method="holm")[1]
    return {
        "status":"source_matched_exploratory_analogue_alignment",
        "n_taxa":len(taxa),
        "sector_taxa_counts":sector_counts,
        "successful_rarefaction_draws":n_qualified,
        "requested_rarefaction_draws":draws,
        "mean_eligible_cross_sector_taxon_pairs":float(np.mean(source_pairs)),
        "min_eligible_cross_sector_taxon_pairs":min(source_pairs),
        "max_eligible_cross_sector_taxon_pairs":max(source_pairs),
        "mean_eligible_analogue_cell_pairs":float(np.mean(
            [q["n_analog_cell_pairs"] for q in supports if q["status"]=="qualified"])),
        "primary_trait_distance_advantage_of_climate_matched_cells":observed,
        "null_mean":expected,
        "null_standard_deviation":float(permuted.std()),
        "two_sided_taxon_block_cell_permutation_p_exploratory":p,
        "holm2_distance_vs_overlap_p":float(p_adjusted[0]),
        "climate_analogue_conditional_KDE_L2_overlap":{
            "status":"estimated",
            "kernel_bandwidth_sd":KDE_BANDWIDTH,
            "mean_conditional_analogue_support":float(np.mean(
                [q["conditional_mean"] for q in overlaps])),
            "mean_unconditional_eight_cells":float(np.mean(
                [q["unconditional_mean"] for q in overlaps])),
            "conditional_minus_all_cells_overlap":float(overall_overlap),
            "null_mean":null_overlap_mean,
            "two_sided_taxon_block_permutation_p":p_overlap,
            "holm2_distance_vs_overlap_p":float(p_adjusted[1]),
            "bandwidth_sensitive_small_cell_KDE_not_true_hypervolume_intersection":True,
            "conditional_clouds_selected_by_climate_only":True,
        },
        "n_permutations":permutations,
        "min_climate_analog_distance":CLIMATE_ANALOG_Z,
        "sampling_sd_of_rarefaction_estimates":float(estimates.std()),
        "observed_positive_means_more_similar_under_climate_matching":True,
        "within_taxon_cell_shuffles_preserve_taxon_mean":True,
        "taxon_pairs_not_independent_replicates":True,
        "not_phylogenetic_evolutionary_convergence_or_adaptation":True,
    }


def run(source: Path, receipt: Path, out: Path):
    frame,provenance=load_source(source,receipt)
    result={}
    for degree in RESOLUTIONS:
        cells,coverage=geo_aggregate(frame,degree)
        if coverage["source_taxa_with_cell_and_photo_support"]<25:
            result[f"{degree}degree"]={"status":"insufficient_source_taxa",
                                      "source":coverage}
            continue
        z,scale=standardise_cells(cells)
        fit=analogue_analysis(z,seed=SEED+degree)
        result[f"{degree}degree"]={
            "status":fit["status"],"source":coverage,"analysis":fit,
            "global_axis_standardisation_equal_taxon_weight":scale["equal_taxon_total_weight"],
        }
    outdoc={
        "schema":"azami_H3c_photo_cell_climate_analogue_morphological_alignment_v1",
        "status":"SOURCE_DERIVED_EXPLORATORY_NOT_ORIGINAL_V4_HOLM",
        "original_source":provenance,
        "registered_design":{
            "climate_dimensions":list(CLIMATE),
            "visible_flower_dimensions":list(PHENOTYPE),
            "each_taxon_rarefied_cells":RAREFIED_CELLS,
            "n_rarefaction_draws":N_DRAWS,
            "climate_analogue_Euclidean_z_threshold":CLIMATE_ANALOG_Z,
            "n_permutations":N_PERMUTATIONS,
            "KDE_L2_overlap_fixed_bandwidth":KDE_BANDWIDTH,
            "two_predeclared_ecological_analogue_metrics_Holm2":True,
            "min_cross_region_taxon_pairs":MIN_ELIGIBLE_TAXON_PAIRS,
            "min_distinct_climatic_region_pairs":MIN_SECTOR_PAIR_TYPES,
            "same_taxon_same_picture_climate_trait_QC":True
        },
        "results":result,
        "limits":{
            "no_tree_or_independent_origins_test":True,
            "no_causal_environmental_filtering_test":True,
            "photographic_angle_not_gravity_calibrated":True,
            "only_one_species_taxon_cross_sector_pair_vote":True,
            "source_cell_and_observer_sampling_remain_potential_bias":True,
            "original_S1_80_S2_128_Holm_still_frozen":True,
            "no_identifiers_coordinates_or_images_published":True,
        }
    }
    dest=Path(out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(outdoc,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                    encoding="utf-8")
    print(json.dumps({
       "status":outdoc["status"],
       "source":provenance,
       "results":{name:d["analysis"] if "analysis" in d else d
                  for name,d in result.items()}
    },allow_nan=False),flush=True)
    return outdoc


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.source,args.receipt,args.output)


if __name__=="__main__":
    main()
