"""H3d: climatic niche breadth versus individual trait diversity, and climate
moderation of *cross-domain residual trait covariance* within Cirsium taxa.

Three photographed traits (camera-angle, corolla chroma, apparent head aspect)
are continuous; three realised climate variables are BIO12, BIO01 and RSDS.
Data are source-verified confidential native S1 observations, reduced to
equal-weight geographic cell medians. Each taxon is an inferential unit.

H3d-1: 3 predeclared 1D trait log-SDs versus climate 3D log-volume after
same taxon sampling controls used by H3a; Holm correction for 3 comparisons.
H3d-2: 3 pairwise residual cross-products after jointly removing additive
within-taxon climate effects, compared between higher/lower *within-taxon*
BIO12 and RSDS cells: 6 tests, taxon-level sign-flip null, Holm correction.
This is observational environment-dependent photographic trait association,
not a statistical fitness interaction, genetic pleiotropy or selection.
No source taxon IDs, photo IDs, observers or coordinates are exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from analysis.v4.s1_climate_trait_hypervolume import (
    AXES, CLIMATE, PHENOTYPE, N_DRAWS, RAREFIED_CELLS,
    SEED, MIN_TAXA, _design, load_source, geo_aggregate,
    standardise_cells, draw_taxon_scores,
)

PAIRS = ((0,1),(0,2),(1,2))
EXPOSURES = (0,2)   # BIO12 annual precipitation; RSDS mean solar radiation
PERMUTATIONS = 999
BOOTSTRAPS = 999
MIN_LOW_HIGH = 3
MIN_TAXA_FOR_COVARIANCE = 25


def taxon_axis_breadth(cells: pd.DataFrame, scores: pd.DataFrame,
                       draws=N_DRAWS, permutations=399, seed=SEED+7):
    """Same 8-cell opportunity as H3a, per photographic dimension, Holm3."""
    rng=np.random.default_rng(seed)
    traits={}
    for taxon,part in cells.groupby("accepted_key",sort=True):
        if len(part)<RAREFIED_CELLS:
            continue
        values=part[["_z3","_z4","_z5"]].to_numpy(float)
        samples=[]
        for k in range(draws):
            ix=rng.choice(len(values),RAREFIED_CELLS,replace=False)
            # A small fixed variance floor avoids negative infinity for
            # species with very little photographed colour or shape variation.
            samples.append(.5*np.log(values[ix].var(axis=0,ddof=1)+0.025))
        traits[str(taxon)]=np.mean(samples,axis=0)
    if len(traits)!=len(scores) or set(traits)!=set(scores._taxon):
        raise ValueError("H3_AXIS_BREADTH_COHORT_NOT_MATCHED_TO_H3_VOLUME")
    x,active,omitted=_design(scores)
    null=x[:,[0,*range(2,x.shape[1])]]
    rng=np.random.default_rng(seed+1)
    sector=scores.sector.to_numpy(int)
    outputs=[]
    for j,name in enumerate(PHENOTYPE):
        y=np.array([traits[t][j] for t in scores._taxon],dtype=float)
        if y.std()<1e-8:
            outputs.append({"trait":name,"status":"unidentifiable_zero_trait_variance"})
            continue
        fitted=sm.OLS(y,x).fit(cov_type="HC3",use_t=False)
        coefficient=float(fitted.params[1])
        lo,hi=map(float,fitted.conf_int()[1])
        nuisance=sm.OLS(y,null).fit()
        baseline=nuisance.fittedvalues
        residual=nuisance.resid
        perm=np.empty(permutations)
        for k in range(permutations):
            shuffled=residual.copy()
            for sec in np.unique(sector):
                ids=np.where(sector==sec)[0]
                shuffled[ids]=rng.permutation(residual[ids])
            perm[k]=np.linalg.lstsq(x,baseline+shuffled,rcond=None)[0][1]
        p=(1+int(np.count_nonzero(np.abs(perm)>=abs(coefficient))))/(permutations+1)
        outputs.append({
            "trait":name,
            "status":"estimated",
            "n_taxa":len(scores),
            "log_trait_sd_slope_per_3d_climate_logvolume":coefficient,
            "HC3_pointwise_95CI":[lo,hi],
            "within_sector_residual_permutation_p":p,
            "rarefaction_draws":draws,
            "same_species_climate_volume_as_original_H3a":True
        })
    valid=[i for i,r in enumerate(outputs) if r["status"]=="estimated"]
    if valid:
        adjusted=multipletests(
            [outputs[i]["within_sector_residual_permutation_p"] for i in valid],
            method="holm")[1]
        for i,p in zip(valid,adjusted):
            outputs[i]["exploratory_Holm_3_traits"]=float(p)
    return {"status":"estimated" if len(valid)==3 else "incomplete",
            "n_original_predeclared_traits":3,
            "controls_retained":active,"controls_omitted":omitted,
            "taxon_units":len(scores),
            "results":outputs,
            "no_phenotypic_plasticity_or_adaptation_test":True}


def pooled_additive_residuals(cells: pd.DataFrame):
    """Remove taxon means and 3 linear main climate axes from each trait.

    Species are equal-weighted to avoid photo-rich lineages determining the
    residualisation. Environmental moderation uses residual cross products,
    not post-hoc selected signed-r associations.
    """
    original=cells.copy()
    taxa=original.accepted_key.astype(str)
    exposure=original[["_z0","_z1","_z2"]].to_numpy(float)
    shape=original[["_z3","_z4","_z5"]].to_numpy(float)
    clim_mu=original[["_z0","_z1","_z2"]].groupby(taxa).transform("mean").to_numpy(float)
    trait_mu=original[["_z3","_z4","_z5"]].groupby(taxa).transform("mean").to_numpy(float)
    xc=exposure-clim_mu
    yc=shape-trait_mu
    count=original.groupby("accepted_key").size()
    weights=original.accepted_key.map(1/count).to_numpy(float)
    if np.linalg.matrix_rank(xc)!=3:
        raise ValueError("THREE_CLIMATE_CONTRAST_COLUMNS_NOT_IDENTIFIABLE")
    wroot=np.sqrt(weights)
    beta=np.linalg.lstsq(xc*wroot[:,None],yc*wroot[:,None],rcond=None)[0]
    residual=yc-xc@beta
    if not np.isfinite(residual).all():
        raise ValueError("INVALID_PHOTO_TRAIT_RESIDUALS")
    return original,residual,beta


def within_taxon_covariance_moderation(cells: pd.DataFrame,
                                      permutations=PERMUTATIONS,
                                      bootstraps=BOOTSTRAPS,
                                      seed=SEED+11):
    """Six fixed tests: 3 cross-domain pairs × BIO12 or RSDS, Holm6.

    A positive difference indicates higher *residual covariance* in cells
    at greater within-taxon climate exposure, not a causal climate effect.
    """
    d,resid,beta=pooled_additive_residuals(cells)
    d=d.copy()
    for i in range(3):
        d["_res"+str(i)]=resid[:,i]
    rows=[]
    for exposure_index in EXPOSURES:
        exposure=CLIMATE[exposure_index]
        diffs={pair:[] for pair in PAIRS}
        taxon_eligible=0
        for name,g in d.groupby("accepted_key",sort=True):
            env=g["_z"+str(exposure_index)].to_numpy(float)
            midpoint=np.median(env)
            hi=env>midpoint
            lo=env<midpoint
            if int(hi.sum())<MIN_LOW_HIGH or int(lo.sum())<MIN_LOW_HIGH:
                continue
            r=g[["_res0","_res1","_res2"]].to_numpy(float)
            taxon_eligible+=1
            for i,j in PAIRS:
                cross=r[:,i]*r[:,j]
                diffs[(i,j)].append(float(cross[hi].mean()-cross[lo].mean()))
        for i,j in PAIRS:
            values=np.asarray(diffs[(i,j)],dtype=float)
            if len(values)<MIN_TAXA_FOR_COVARIANCE:
                rows.append({"status":"insufficient_taxa","climate":exposure,
                             "trait_a":PHENOTYPE[i],"trait_b":PHENOTYPE[j],
                             "n_taxa":len(values)})
                continue
            rng=np.random.default_rng(seed+exposure_index*100+i*10+j)
            estimate=float(values.mean())
            signed=(rng.integers(0,2,size=(permutations,len(values)))*2-1)
            null=(signed*values).mean(axis=1)
            p=(1+int((np.abs(null)>=abs(estimate)).sum()))/(permutations+1)
            draws=rng.integers(0,len(values),size=(bootstraps,len(values)))
            ci=np.quantile(values[draws].mean(axis=1),[.025,.975])
            rows.append({
                "status":"estimated","climate":exposure,
                "trait_a":PHENOTYPE[i],"trait_b":PHENOTYPE[j],
                "n_taxa":len(values),
                "higher_minus_lower_exposure_residual_covariance":estimate,
                "taxon_bootstrap_pointwise_95CI":[float(ci[0]),float(ci[1])],
                "taxon_sign_flip_permutation_p":float(p),
                "n_permutations":permutations,
                "taxon_sign_flip_not_phylogenetic_independence":True,
            })
    valid=[i for i,r in enumerate(rows) if r["status"]=="estimated"]
    if valid:
        adjusted=multipletests(
            [rows[i]["taxon_sign_flip_permutation_p"] for i in valid],
            method="holm")[1]
        for i,p in zip(valid,adjusted):
            rows[i]["exploratory_Holm_6_moderations"]=float(p)
    return {"status":"estimated" if len(valid)==6 else "incomplete",
            "n_tests_registered":6,
            "n_estimable_tests":len(valid),
            "n_cells_source":len(cells),
            "n_taxa_source":int(cells.accepted_key.nunique()),
            "n_eligible_taxa_by_exposure":{
                exposure:next((r["n_taxa"] for r in rows if r["climate"]==exposure),0)
                for exposure in (CLIMATE[0],CLIMATE[2])},
            "trait_additive_climate_coefficients_matrices":{
                name: [float(v) for v in beta[:,i]]
                for i,name in enumerate(PHENOTYPE)
            },
            "results":rows,
            "regional_and_phylogenetic_effects_not_causally_separated":True,
            "photo_measured_traits_not_fitness_or_genetic_covar":True,
    }


def execute(source,receipt,output):
    frame,provenance=load_source(source,receipt)
    results={}
    for degree in (1,5):
        cells,coverage=geo_aggregate(frame,degree)
        if coverage["source_taxa_with_cell_and_photo_support"]<MIN_TAXA:
            results[f"{degree}degree"]={"status":"insufficient_source_support",
                                        "source":coverage}
            continue
        cells,_=standardise_cells(cells)
        scores=draw_taxon_scores(cells)
        response={
            "status":"estimated",
            "source":coverage,
            "n_species":len(scores),
            "n_source_cells":len(cells),
            "trait_specific_diversity_vs_climate_breadth":
                taxon_axis_breadth(cells,scores),
            "environment_modulated_residual_trait_covariation":
                within_taxon_covariance_moderation(cells),
        }
        results[f"{degree}degree"]=response
    payload={
        "schema":"azami_s1_h3d_environment_moderated_photo_trait_covariance_v1",
        "status":"SOURCE_VERIFIED_EXPLORATORY_H3_NOT_ORIGINAL_HOLM",
        "original_source":provenance,
        "dimensions":{"climate":list(CLIMATE),
                      "phenotype":list(PHENOTYPE)},
        "results":results,
        "scope":{
            "one_degree_registered_primary_five_degree_secondary":True,
            "rareify_to_eight_cells_same_sampling_opportunity_as_H3a":True,
            "three_1d_breadth_tests_Holm3_and_six_moderation_tests_Holm6":True,
            "observational_conditional_trait_covariance_not_functional_interaction":True,
            "no_original_208_tests_refitted":True,
            "no_taxon_observation_geographic_or_photographer_ids_exported":True,
        },
    }
    dest=Path(output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(payload,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                    encoding="utf-8")
    print(json.dumps({"status":payload["status"],
                      "source":provenance,
                      "summary":{k:{
                          "taxa":v.get("n_species"),
                          "breadth":v.get("trait_specific_diversity_vs_climate_breadth"),
                          "moderation":v.get("environment_modulated_residual_trait_covariation")
                      } if v.get("status")=="estimated" else v for k,v in results.items()}
                      },ensure_ascii=False,allow_nan=False),flush=True)
    return payload


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    execute(args.source,args.receipt,args.out)


if __name__=="__main__":
    main()
