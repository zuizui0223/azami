"""Native Gate C robustness: 30 matched rarefactions, bandwidth, taxon-block
uncertainty and observer/geographic source-deletion sensitivities.

Operates on exactly the prior authenticated 313350-row native source and the
same prephenotype WCVP/observer/cell support. No image acquisition, new source
selection, environmental re-fitting, or publication of taxon/observer/obs IDs.
All uncertainty is exploratory: selected taxon/phylogeny and camera error persist.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analysis.v4.native_gates_ae import authority_mapping, gate_a
from analysis.v4.native_gates_bcd import (
    MODULES, _group_scores, _density_overlap, fixed_sample, rng_for,
    wasserstein_parts,
)
from analysis.v4.run_native_abcde import (
    EXPECTED_SHA, validate_inputs, FIXED_PREDICTORS,
)

SEED = 20261021
N_PER_TAXON = 50
RANDOM_LONGITUDE_SECTORS = (-180, -90, 0, 90, 180)


def source_scores(candidate: pd.DataFrame, taxa: list[str]) -> tuple[dict, dict]:
    """Fixed measurement-only PCA basis at same first-50/taxon design as A-E."""
    base=fixed_sample(candidate,n=N_PER_TAXON,replicate=0)
    label=candidate.accepted_key.astype(str).to_numpy()
    modules,retained={},{}
    for name,columns in MODULES.items():
        dim=1 if name=="presentation" else 3
        score,r=_group_scores(candidate,base,columns,dim)
        modules[name]={taxon:score[label==taxon] for taxon in taxa}
        retained[name]=float(r)
    return modules,retained


def _reference(values:np.ndarray,tag:str,taxon:str,
               reps:int,n:int=N_PER_TAXON)->dict|None:
    if len(values)<2*n:
        return None
    full,shape=[],[]
    for rep in range(reps):
        rng=rng_for(SEED,"reference",tag,taxon,rep)
        chosen=rng.choice(len(values),size=2*n,replace=False)
        w,_,s=wasserstein_parts(values[chosen[:n]],values[chosen[n:]])
        full.append(w);shape.append(s)
    return {"full":float(np.quantile(full,.95)),
            "shape":float(np.quantile(shape,.95))}


def pair_metrics(modules:dict, taxa:list[str],*,reps:int,
                 null_reps:int,n:int=N_PER_TAXON)->tuple[dict,dict]:
    """All modules use the SAME taxon-paired source image draws per replicate."""
    pair_index=list(itertools.combinations(range(len(taxa)),2))
    result,private={},{}
    for name,data in modules.items():
        null={t:_reference(data[t],name,t,null_reps,n=n) for t in taxa}
        shape=np.empty((len(pair_index),reps),float)
        full=np.empty((len(pair_index),reps),float)
        center=np.empty((len(pair_index),reps),float)
        for rep in range(reps):
            sampled={}
            for t in taxa:
                x=data[t]
                if len(x)<n:
                    raise ValueError("Frozen Gate A sampling no longer supports fixed n")
                rng=rng_for(SEED,"paired-source-rows",t,rep)
                sampled[t]=x[rng.choice(len(x),size=n,replace=False)]
            for ix,(i,j) in enumerate(pair_index):
                w,c,s=wasserstein_parts(sampled[taxa[i]],sampled[taxa[j]])
                full[ix,rep]=w;center[ix,rep]=c;shape[ix,rep]=s
        median_shape=np.median(shape,axis=1)
        median_full=np.median(full,axis=1)
        matched=np.array([
            null[taxa[i]] is not None and null[taxa[j]] is not None
            for i,j in pair_index],bool)
        shape_flag=np.zeros(len(pair_index),dtype=bool)
        full_flag=np.zeros(len(pair_index),dtype=bool)
        if matched.any():
            targets=np.flatnonzero(matched)
            for ix in targets:
                i,j=pair_index[ix]
                shape_flag[ix]=median_shape[ix]>max(null[taxa[i]]["shape"],null[taxa[j]]["shape"])
                full_flag[ix]=median_full[ix]>max(null[taxa[i]]["full"],null[taxa[j]]["full"])
        n_referenced=int(matched.sum())
        result[name]={
            "dimensions":next(iter(data.values())).shape[1],
            "total_pairs":len(pair_index),"taxa_with_disjoint_50_plus_50_null":sum(x is not None for x in null.values()),
            "pairs_with_both_taxa_null":n_referenced,
            "rarefaction_replicates":reps,"same_taxon_null_replicates":null_reps,
            "median_over_pairs_of_30_replicate_median_full_W2":float(np.median(median_full)),
            "median_over_pairs_of_30_replicate_median_centred_W2":float(np.median(median_shape)),
            "median_centroid_W2":float(np.median(np.median(center,axis=1))),
            "reference_exceedances_full_W2":int(full_flag.sum()),
            "reference_exceedances_centred_W2":int(shape_flag.sum()),
            "centred_exceedance_fraction":float(shape_flag.sum()/n_referenced) if n_referenced else None,
            "full_exceedance_fraction":float(full_flag.sum()/n_referenced) if n_referenced else None,
            "reference_uses_pair_median_not_first_rep_only":True,
        }
        private[name]={"reference":null,"pair_index":pair_index,"matched":matched,
                       "centred_flags":shape_flag,"full_flags":full_flag}
    return result,private


def taxon_paired_difference(private:dict,taxa:list[str],candidate:pd.DataFrame,
                            *,bootstrap:int=1000)->dict:
    """Classify pair before bootstrapping taxa, never treat shared-taxon pairs iid."""
    a,b=private["colour"],private["outline"]
    if a["pair_index"]!=b["pair_index"] or not np.array_equal(a["matched"],b["matched"]):
        raise ValueError("Colour and outline do not share same original taxon pairs")
    target=np.flatnonzero(a["matched"])
    if len(target)<10:
        return {"status":"too_few_same_taxon_null_pairs"}
    ta=np.array(taxa)
    pairs=np.array(a["pair_index"],dtype=int)[target]
    colour=a["centred_flags"][target].astype(float)
    outline=b["centred_flags"][target].astype(float)
    n=len(taxa)
    reference_taxon_indices=np.array(
        [i for i,t in enumerate(taxa) if a["reference"][t] is not None],dtype=int
    )
    nc=int(len(reference_taxon_indices))
    diff=float(np.mean(colour-outline))
    both=int(np.sum((colour==1)&(outline==1)))
    colour_only=int(np.sum((colour==1)&(outline==0)))
    outline_only=int(np.sum((colour==0)&(outline==1)))
    neither=int(np.sum((colour==0)&(outline==0)))
    if sum((both,colour_only,outline_only,neither))!=len(target):
        raise AssertionError("Paired contrast pair-count identity failed")
    # Re-sample original taxonomy units, not the many dependent taxon pairs.
    rng=rng_for(SEED,"paired-colour-vs-outline-taxon-cluster-bootstrap")
    draws=[]
    for rep in range(bootstrap):
        w=np.bincount(rng.choice(reference_taxon_indices,size=nc,replace=True),
                      minlength=n)
        ew=w[pairs[:,0]]*w[pairs[:,1]]
        den=ew.sum()
        if den>0:
            draws.append(float(ew@(colour-outline)/den))
    if len(draws)<.95*bootstrap:
        raise ValueError("Insufficient valid taxon block resamples")
    # Nonparametric leave-one-taxon sensitivity, same paired comparisons.
    loo=[]
    for k in range(n):
        keep=(pairs[:,0]!=k)&(pairs[:,1]!=k)
        if keep.any():
            loo.append(float(np.mean(colour[keep]-outline[keep])))
    # Already selected, high-observation taxa exclusion holds bases, samples,
    # denominators and original species pair comparisons fixed.
    abundance=candidate.accepted_key.astype(str).value_counts()
    ordered=sorted(taxa,key=lambda t:(-int(abundance[t]),t))
    omit={taxa.index(t) for t in ordered[:2]}
    no_top=np.array([i not in omit and j not in omit for i,j in pairs],bool)
    if no_top.sum()<10:
        raise ValueError("No valid pair support after deleting top taxa")
    series=np.asarray(draws,float)
    return {
        "status":"paired_colour_vs_outline_same_38_taxa_exploratory",
        "shared_species_pairs":len(target),
        "reference_supported_taxa":nc,
        "centred_W2_colour_above_null":int(colour.sum()),
        "centred_W2_outline_above_null":int(outline.sum()),
        "pair_discordance":{"both_above":both,"colour_only":colour_only,
                            "outline_only":outline_only,"neither_above":neither},
        "observed_difference_fraction_colour_minus_outline":diff,
        "taxon_cluster_bootstrap_valid":len(draws),
        "taxon_cluster_bootstrap_q025_q50_q975":[
            float(x) for x in np.quantile(series,[.025,.5,.975])
        ],
        "taxon_cluster_bootstrap_positive_fraction":float(np.mean(series>0)),
        "leave_one_taxon_out_difference_min_max":[float(min(loo)),float(max(loo))],
        "exclude_two_most_reported_taxa":{
            "remaining_pairs":int(no_top.sum()),
            "difference_fraction_colour_minus_outline":float(
                np.mean((colour-outline)[no_top])
            ),
            "selection":"original source reporting count, outcome-blind",
        },
        "uncertainty_scope":"Taxon-cluster resampling conditional on the 40 selected taxa, the same camera-derived shared PCA coordinates and n=50. It does not address phylogenetic or photographic correlation and is not a calibrated p-value.",
    }


def bandwidth_sensitivity(modules:dict,taxa:list[str],*,pairs_sample:int=100,
                          reps:int=10,draws:int=128)->dict:
    if reps<2 or pairs_sample<2 or draws<32:
        raise ValueError("Bandwidth numerical settings are too small")
    pairs=list(itertools.combinations(taxa,2))
    k=min(pairs_sample,len(pairs))
    indices=np.sort(rng_for(SEED,"outcome-blind-bandwidth-pair-selection").choice(
        len(pairs),size=k,replace=False))
    selected=[pairs[int(i)] for i in indices]
    result={}
    for name,data in modules.items():
        dimensions=next(iter(data.values())).shape[1]
        h=N_PER_TAXON**(-1/(dimensions+4))
        values=np.empty((k,reps,3),float)
        for rep in range(reps):
            involved=sorted({t for a,b in selected for t in (a,b)})
            samples={}
            for t in involved:
                r=rng_for(SEED,"paired-source-rows",t,rep)
                samples[t]=data[t][r.choice(len(data[t]),size=N_PER_TAXON,replace=False)]
            for ix,(a,b) in enumerate(selected):
                seed=int(rng_for(SEED,"matched-bandwidth-monte-carlo",name,a,b,rep
                       ).integers(0,2**31))
                for ih,mult in enumerate((.75,1.,1.25)):
                    values[ix,rep,ih]=_density_overlap(
                        samples[a],samples[b],h*mult,draws=draws,seed=seed)
        paired=np.median(values,axis=1)
        if not np.isfinite(paired).all() or (paired<0).any() or (paired>1).any():
            raise ValueError("Density estimator returned invalid mass overlap")
        rho={}
        for col,label in ((0,"0.75"),(2,"1.25")):
            stat=float(spearmanr(paired[:,1],paired[:,col]).statistic)
            rho[label]=stat if np.isfinite(stat) else None
        result[name]={
            "shared_sampled_taxon_pairs":k,
            "paired_rarefaction_replicates":reps,
            "monte_carlo_draws_per_distribution":draws,
            "KDE_bandwidth_multipliers":[.75,1.,1.25],
            "median_pairwise_KDE_overlap":[float(x) for x in np.median(paired,axis=0)],
            "rank_spearman_vs_base":rho,
            "base_bandwidth_h":float(h),
            "same_pairs_and_same_30_based_subsampling_seed_across_bandwidths":True,
        }
    return {
        "status":"completed_matched_bandwidth_screen",
        "same_taxon_pair_indices_for_all_three_modules":True,
        "modules":result,
        "scope":"100 outcome-blind pairs, not every pair; 10 paired repetitions and fixed normalized bandwidths. Spearman rank robustness is not absolute overlap calibration.",
    }


def deletion_sensitivity(candidate:pd.DataFrame,modules:dict,taxa:list[str],
                         *,n:int=N_PER_TAXON,reps:int=6,null_reps:int=10)->dict:
    """Rerun actual fixed-basis colour/outline W2 after observer/region deletion."""
    original=pd.DataFrame({"taxon":candidate.accepted_key.astype(str),
                           "observer":candidate.user_id.astype(str)})
    if not original.observer.notna().all():
        raise ValueError("Missing private original photographer strata")
    abundance=original.observer.value_counts()
    top=max(1,int(math.ceil(.01*len(abundance))))
    top_people=set(abundance.index[:top])
    lon=candidate.longitude.to_numpy(float)
    masks={"top_one_pct_observers_removed":~original.observer.isin(top_people).to_numpy()}
    for k,(low,high) in enumerate(zip(RANDOM_LONGITUDE_SECTORS[:-1],
                                      RANDOM_LONGITUDE_SECTORS[1:])):
        in_sector=(lon>=low)&(lon<high)
        masks["exclude_fixed_longitude_sector_"+str(k)]=~in_sector
    names=original.taxon.to_numpy()
    result={}
    for label,mask in masks.items():
        taxa_eligible=[t for t in taxa if int(np.sum(mask&(names==t)))>=2*n]
        count_obs=int(mask.sum())
        if len(taxa_eligible)<10:
            result[label]={"status":"not_comparable_under_predeclared_min_ten_taxa",
                           "remaining_candidate_rows":count_obs,
                           "taxa_with_100_records_after_exclusion":len(taxa_eligible)}
            continue
        arrays={}
        for module in ("colour","outline"):
            arrays[module]={t:modules[module][t][mask[names==t]] for t in taxa_eligible}
        summary,private=pair_metrics(arrays,taxa_eligible,reps=reps,null_reps=null_reps,n=n)
        a=summary["colour"];b=summary["outline"]
        if a["pairs_with_both_taxa_null"]!=b["pairs_with_both_taxa_null"]:
            raise ValueError("Cross-module deletion sensitivity has different taxon universe")
        result[label]={
            "status":"executed_fixed_projection_two_module_subcohort_screen",
            "remaining_candidate_rows":count_obs,
            "taxa_with_100_records_after_exclusion":len(taxa_eligible),
            "pairs_evaluated":a["pairs_with_both_taxa_null"],
            "colour_centred_exceedance_fraction":a["centred_exceedance_fraction"],
            "outline_centred_exceedance_fraction":b["centred_exceedance_fraction"],
            "colour_minus_outline_difference":float(
                a["centred_exceedance_fraction"]-b["centred_exceedance_fraction"]
            ),
            "source_subset_defines_taxa":"taxa dropped on insufficient remaining records; not strictly matched to full 38-taxon primary",
            "PCA_basis":"frozen from original full selected 40-taxon Gate A frame",
            "source_deletion_reps":reps,
            "within_species_null_reps":null_reps,
        }
    return {
        "status":"five_a_priori_source_deletion_analyses",
        "original_original_observer_percentile_removed":.01,
        "original_observers_removed_count":top,
        "fixed_longitude_sectors_degrees":list(RANDOM_LONGITUDE_SECTORS),
        "scenarios":result,
        "limitations":"Geographic quadrants are not independent clades or spatially buffered holdouts; scenario taxon subsets differ from main estimator when n<100.",
    }


def analyse(candidate:pd.DataFrame,*,reps:int=30,null_reps:int=30,
            bandwidth_pairs:int=100,bandwidth_reps:int=10,
            bootstrap:int=1000,scenario_reps:int=6)->dict:
    required={"accepted_key","user_id","latitude","longitude","_cell","obs_id"}
    if not required<=set(candidate):
        raise ValueError("Full Gate A candidate and original observer identity required")
    taxa=sorted(candidate.accepted_key.astype(str).unique().tolist())
    if len(taxa)!=40:
        raise ValueError("Primary registered native taxa changed; do not move the goalposts")
    if min(candidate.accepted_key.astype(str).value_counts())<N_PER_TAXON:
        raise ValueError("Previously selected 40 taxa no longer meet fixed 50 rarefaction")
    if reps<6 or null_reps<8 or bootstrap<40 or bandwidth_reps<3:
        raise ValueError("Insufficient declared robustness audit repetitions")
    modules,retained=source_scores(candidate,taxa)
    main,private=pair_metrics(modules,taxa,reps=reps,null_reps=null_reps)
    paired=taxon_paired_difference(private,taxa,candidate,bootstrap=bootstrap)
    bw=bandwidth_sensitivity(modules,taxa,pairs_sample=bandwidth_pairs,
                             reps=bandwidth_reps)
    deletions=deletion_sensitivity(candidate,modules,taxa,reps=scenario_reps)
    return {
        "status":"executed_same_40_taxa_native_colour_vs_outline_robustness",
        "fixed_sampling_plan":{"same_original_native_frame":True,
            "same_registered_40_taxa":True,
            "sample_per_taxon":N_PER_TAXON,
            "original_human_observer_x_spatial_cell_strata":True,
            "same_source_rows_for_all_module_rarefactions":True,
            "pairwise_rare_replicates":reps,
            "same_species_null_replicates":null_reps,
            "taxon_cluster_bootstrap":bootstrap,
            "selection_uses_outcome":False,
            "retained_shared_projection_variance":retained,
        },
        "modules":main,
        "colour_vs_outline":paired,
        "KDE_bandwidth_sensitivity":bw,
        "observer_region_deletions":deletions,
        "privacy":"Only aggregate counts, effect summaries, source SHA and fixed design constants; no taxon, observer, source observation or geographic point ID exported.",
        "scientific_limits":[
            "Taxon bootstrap is conditional on selection of 40 well-sampled taxa.",
            "All taxon pairs share species, and source photos may repeat individual plants.",
            "Differing module projection dimensions prevent raw orientation-versus-colour overlap test.",
            "Camera pose, illumination and colour/photo technical uncertainty are not calibrated.",
            "Same-species distance thresholds are descriptive finite-n references, not p values.",
            "The genetic ancestry and trait function link must be established by independent EAzami/aza3 data.",
        ],
    }


def execute(table:Path,receipt:Path,models:Path,integration:Path,spatial:Path,
            mapping:Path,observers:Path,baseline:Path,output:Path,
            *,reps=30,null_reps=30,bootstrap=1000,scenario_reps=6)->dict:
    validate_inputs(table,receipt,models,integration,spatial)
    data=pd.read_csv(table,dtype={"obs_id":str,"accepted_key":str},low_memory=False)
    obs=pd.read_csv(observers,dtype={"obs_id":str,"user_id":str},low_memory=False)
    concepts=authority_mapping(mapping)
    prior=json.loads(baseline.read_text(encoding="utf-8"))
    a,candidate=gate_a(data,concepts,obs,max_taxa=40)
    if a["status"]!="observer_cell_balanced_selection_constructed":
        raise ValueError("No original-observer-verified common native A gate")
    if prior["A"]["predeclared_selection"]["selected_taxa"]!=40:
        raise ValueError("Original integrated run had a different taxon-support gate")
    if (a["matched_nine_output_observations"],
        a["matched_rows_in_collided_accepted_keys"],
        a["predeclared_selection"]["selected_nonredundant_observations"])!=(
        prior["A"]["matched_nine_output_observations"],
        prior["A"]["matched_rows_in_collided_accepted_keys"],
        prior["A"]["predeclared_selection"]["selected_nonredundant_observations"]
    ):
        raise ValueError("Reconstructed Gate A cohort is not prior real-data sample")
    result=analyse(candidate,reps=reps,null_reps=null_reps,
                   bootstrap=bootstrap,scenario_reps=scenario_reps)
    result["input_provenance"]={
        "verified_native_csv_sha256":EXPECTED_SHA,
        "gate_a_prior_source_sha256":prior["native_source_sha256"],
        "source_integrated_baseline_SHA256":hashlib.sha256(baseline.read_bytes()).hexdigest(),
        "predeclared_full_species_pairs":780,
        "privacy":"NO ORIGINAL RECORDS OR USER IDENTIFIERS IN RESULT",
    }
    output.mkdir(parents=True,exist_ok=True)
    (output/"native_density_robustness.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
        encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "colour_outline_difference":result["colour_vs_outline"]["observed_difference_fraction_colour_minus_outline"],
        "referenced_pairs":result["colour_vs_outline"]["shared_species_pairs"],
        "rarefaction_reps":reps,
    }),flush=True)
    return result


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("table","receipt","models","integration","spatial","mapping",
                "observers","baseline","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--reps",type=int,default=30)
    p.add_argument("--null-reps",type=int,default=30)
    p.add_argument("--bootstrap",type=int,default=1000)
    p.add_argument("--scenario-reps",type=int,default=6)
    a=p.parse_args()
    execute(a.table,a.receipt,a.models,a.integration,a.spatial,a.mapping,
            a.observers,a.baseline,a.output,reps=a.reps,null_reps=a.null_reps,
            bootstrap=a.bootstrap,scenario_reps=a.scenario_reps)


if __name__=="__main__":
    main()
