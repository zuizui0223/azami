"""Does the original S1 continuous head-image axis carry identifiable history?

Pagel-lambda profile on exactly the verified S1 species-angle cohort:
194 mapped taxa; 191 after three source-concept collisions; plus pre-existing
>=5 (171) and >=10 (143) image/taxon data-depth sensitivities.
Both species MEAN and MEDIAN photo angle, with/without the SAME S1 five
CHELSA exposures (BIO1, BIO4, BIO12, BIO15, RSDS).

All 52 GBOTB/V.PhyloMaker2 trees are grafting-placement scenarios, not
52 independent phylogenetic replicates. Testing signal is a PRIOR requirement
for interpreting conditional Brownian branch-change counts as lineage history.
Neither signal, lack of signal nor any photo-angle change identifies
gravity-referenced biological nodding gains or climate-dependent selection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES
from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.s1_continuous_phylogenetic_history import (
    SHA, audited_pruned_tree, shared_path_matrix,
)

SOURCE_RUN=39591120377
EXCLUDED_COLLISIONS={"Cirsium_parryi","Cirsium_setidens","Cirsium_yezoense"}
NUGGETS=(0.05,0.25)


def verify_source(path):
    raw=json.loads(Path(path).read_text(encoding="utf-8"))
    audit=raw.get("source_audit",{})
    if (audit.get("source_sha256")!=SHA
            or not audit.get("native_observation_source_authenticated")
            or not audit.get("nonindependent_52_tree_ensemble")
            or raw.get("summary",{}).get("n_tree_taxon_keys_with_original_continuous_orientation")!=194
            or raw.get("summary",{}).get("n_original_images_in_S1")!=170066
            or raw.get("summary",{}).get("S1_original_five_exposures")!=list(EXPOSURES)):
        raise ValueError("Original S1 continuous-angle source integrity/coverage failed")
    rows=raw.get("source_cohort_194_species",[])
    if len(rows)!=194:
        raise ValueError("194 source-authenticated species-angle vectors not available")
    names=[x["tip"] for x in rows]
    keys=[str(x["accepted_key"]) for x in rows]
    if len(set(names))!=194 or len(set(keys))!=194:
        raise ValueError("Non-unique S1 phylogenetic tip or WCVP accepted key")
    collisions={x["tip"] for x in rows if x["source_concept_collision"]}
    if collisions != EXCLUDED_COLLISIONS:
        raise ValueError("Original three unresolved source-name collisions changed")
    if sum(int(x["n_eligible_images"]) for x in rows)!=169912:
        raise ValueError("S1 194-tree-covered image count changed")
    for x in rows:
        angles=(float(x["image_mean_deg"]),float(x["image_median_deg"]))
        if any(not np.isfinite(y) or y<0 or y>180 for y in angles):
            raise ValueError("Image-vertical angle absent/outside S1 measurement ontology")
        climate=x.get("S1_five_climate_means")
        if not isinstance(climate,dict) or set(climate)!=set(EXPOSURES):
            raise ValueError("S1 five climate source means absent / different selection")
        vals=[float(climate[e]) for e in EXPOSURES]
        if not np.isfinite(vals).all() or vals[2]<=0:
            raise ValueError("Invalid five original CHELSA source values")
        if int(x["n_eligible_images"])<1:
            raise ValueError("Invalid source image count")
    return rows


def taxon_sets(rows):
    clean=[x for x in rows if not x["source_concept_collision"]]
    return {
      "original194":rows,
      "collision_free191":clean,
      "collision_free_at_least5": [x for x in clean if x["n_eligible_images"]>=5],
      "collision_free_at_least10":[x for x in clean if x["n_eligible_images"]>=10],
    }


def model_design(rows,adjust_climate):
    ordered=sorted(rows,key=lambda row:row["tip"])
    n=len(ordered)
    if n<7:
        raise ValueError("Insufficient species for five-environment correction")
    five=np.array([[float(r["S1_five_climate_means"][v]) for v in EXPOSURES]
                   for r in ordered])
    if not adjust_climate:
        return ordered,np.ones((n,1))
    sd=five.std(axis=0)
    if np.any(sd<=1e-12):
        raise ValueError("Constant S1 selected climate")
    X=np.column_stack([np.ones(n), (five-five.mean(axis=0))/sd])
    if np.linalg.matrix_rank(X)!=6 or np.linalg.cond(X)>1e6:
        raise ValueError("Five CHELSA environmental covariates are collinear")
    return ordered,X


class LambdaProfile:
    """GLS profile using one covariance eigensystem for all model evaluations."""
    def __init__(self,C,X,nugget=.05):
        C=np.asarray(C,float)
        X=np.asarray(X,float)
        n=len(C)
        if (C.shape!=(n,n) or X.shape[0]!=n or n<=X.shape[1]+2
                or nugget<0 or nugget>1
                or not np.isfinite(C).all() or not np.isfinite(X).all()):
            raise ValueError("Invalid source phylogeny / fixed model design")
        C=(C+C.T)/2
        diagonal=np.diag(C)
        if np.any(diagonal<=0):
            raise ValueError("No positive inherited phylogeny depth")
        self.depth=np.sqrt(diagonal)
        corr=C/self.depth[:,None]/self.depth[None,:]
        self.eig,self.rotation=np.linalg.eigh((corr+corr.T)/2)
        self.UX=self.rotation.T@(X/self.depth[:,None])
        self.X=X
        self.n=n
        self.p=X.shape[1]
        self.nugget=nugget
        self.logdet_diag=float(np.log(diagonal).sum())
        self.profile_Xcache={}
        if np.min(self.eig)<-1e-7:
            raise ValueError("Historical Brownian correlation matrix not PSD")

    def fit_at(self,lam,Uy):
        if lam<0 or lam>1:
            raise ValueError("Pagel lambda constrained to [0,1]")
        v=(1+self.nugget)+lam*(self.eig-1)
        if np.min(v)<=0:
            raise ValueError("Invalid covariance at fixed lambda")
        weighted_X=self.UX/np.sqrt(v[:,None])
        weighted_y=Uy/np.sqrt(v)
        beta=np.linalg.lstsq(weighted_X,weighted_y,rcond=None)[0]
        ssr=float(np.sum((weighted_y-weighted_X@beta)**2))
        if ssr<=1e-14:
            raise ValueError("Degenerate S1 source trait residual")
        sigma2=ssr/self.n
        loglike=-.5*(self.n*(np.log(2*np.pi)+1+np.log(sigma2))
                     +self.logdet_diag+float(np.log(v).sum()))
        return float(loglike),beta,float(sigma2)

    def fit(self,y):
        yy=np.asarray(y,float)
        if yy.shape!=(self.n,) or not np.isfinite(yy).all():
            raise ValueError("Expected one continuous original S1 angle per tip")
        Uy=self.rotation.T@(yy/self.depth)
        fit0=self.fit_at(0.,Uy)
        grid=np.linspace(0,1,21)
        candidates=[(float(x),self.fit_at(float(x),Uy)) for x in grid]
        anchor=max(range(len(candidates)),key=lambda i:candidates[i][1][0])
        lo=max(anchor-1,0);hi=min(anchor+1,len(grid)-1)
        opt=minimize_scalar(lambda z:-self.fit_at(float(z),Uy)[0],
             bounds=(float(grid[lo]),float(grid[hi])),method="bounded",
             options={"xatol":1e-5})
        if opt.success:
            candidates.append((float(opt.x),self.fit_at(float(opt.x),Uy)))
        lam,best=max(candidates,key=lambda x:x[1][0])
        lrt=max(0.,2*(best[0]-fit0[0]))
        return {
          "lambda_ML":float(lam),
          "log_likelihood_lambda0":fit0[0],
          "log_likelihood_best":best[0],
          "likelihood_ratio_against_lambda0":float(lrt),
          "deltaAIC_best_vs_lambda0":float(lrt-2),
          "ML_sigma2":best[2],
          "beta_at_lambda0":fit0[1].tolist(),
          "sigma2_at_lambda0":fit0[2],
          "model_rank":self.p,
          "n_taxa":self.n,
        }

    def bootstrap_lambda0(self,observed_y,reps=399,seed=20261021):
        """Parametric H0: independent tips after fixed X, refitting λ each time."""
        if reps<20:
            raise ValueError("At least 20 replicates needed to assess LRT tails")
        observed=self.fit(observed_y)
        rng=np.random.default_rng(seed)
        beta0=np.array(observed["beta_at_lambda0"])
        sigma0=np.sqrt(observed["sigma2_at_lambda0"])
        expected=self.X@beta0
        rootvar=self.depth * np.sqrt(1+self.nugget) * sigma0
        null=np.empty(reps)
        for i in range(reps):
            simulate=expected+rootvar*rng.normal(size=self.n)
            null[i]=self.fit(simulate)["likelihood_ratio_against_lambda0"]
        p=float((1+np.count_nonzero(null >=
                observed["likelihood_ratio_against_lambda0"]-1e-10))/(reps+1))
        return {
          "null":"parametric_bootstrap_phylogenetically_independent_tips_given_original_five_climates",
          "requested_draws":reps,"seed":seed,
          "observed_lrt":observed["likelihood_ratio_against_lambda0"],
          "boundary_calibrated_bootstrap_p":p,
          "null_lrt_q50_q95_q99":np.quantile(null,[.50,.95,.99]).tolist(),
          "observed_lambda_ML":observed["lambda_ML"],
          "scope":"fixed one V.PhyloMaker2/GBOTB topology, one fixed nugget, normal residuals; does not cover history across 52 placement hypotheses or observer-level photo variance",
        }


def covariance_for(tree,rows):
    ordered=sorted(rows,key=lambda x:x["tip"])
    names=[r["tip"] for r in ordered]
    pruned,_=audited_pruned_tree(tree,names)
    C=shared_path_matrix(pruned,names)[0]
    return C


def analysis(rows,trees,*,nuggets=NUGGETS,bootstrap_reps=399):
    sets=taxon_sets(rows)
    expected={"original194":194,"collision_free191":191,
              "collision_free_at_least5":171,"collision_free_at_least10":143}
    if {k:len(v) for k,v in sets.items()}!=expected:
        raise ValueError("S1 original 194/191/171/143 image-precision ladder changed")
    if len(trees)!=52:
        raise ValueError("Expected exact 52 frozen placement hypotheses")
    details=[]
    bootstraps={}
    for cohort_name,species in sets.items():
        ordered,X0=model_design(species,False)
        _,X5=model_design(species,True)
        values={
          "species_mean":np.array([r["image_mean_deg"] for r in ordered],float),
          "species_median":np.array([r["image_median_deg"] for r in ordered],float),
        }
        names=[r["tip"] for r in ordered]
        for index,(tree_name,tree) in enumerate(trees):
            C=covariance_for(tree,ordered)
            for nugget in nuggets:
                for design_label,design in (("intercept_only",X0),("S1_five_climate_adjusted",X5)):
                    fitter=LambdaProfile(C,design,nugget)
                    for phenotype,yy in values.items():
                        result=fitter.fit(yy)
                        details.append({
                          "cohort":cohort_name,"n_taxa":len(species),
                          "tree_scenario":tree_name,"tree_scenario_index":index,
                          "phenotype":phenotype,"model":design_label,
                          "nugget_fraction_of_tip_depth":float(nugget),
                          "lambda_ML":result["lambda_ML"],
                          "likelihood_ratio_against_lambda0":result["likelihood_ratio_against_lambda0"],
                          "deltaAIC_best_vs_lambda0":result["deltaAIC_best_vs_lambda0"],
                        })
                        # Primary diagnostic calibrated at the first frozen
                        # scenario for a clean concept (191) and species with
                        # >=10 original photos (143): no outcome-dependent picks.
                        if (index==0 and abs(nugget-.05)<1e-12 and phenotype=="species_mean"
                                and ((cohort_name=="collision_free191" and
                                      design_label in ("intercept_only","S1_five_climate_adjusted"))
                                  or (cohort_name=="collision_free_at_least10" and
                                      design_label=="S1_five_climate_adjusted"))):
                            label=f"{cohort_name}__{design_label}__{phenotype}"
                            bootstraps[label]=fitter.bootstrap_lambda0(
                                 yy,reps=bootstrap_reps,seed=20261021+len(bootstraps))
    def group_stats(group):
        lam=np.array([x["lambda_ML"] for x in group])
        lrt=np.array([x["likelihood_ratio_against_lambda0"] for x in group])
        return {
            "n_nonindependent_placement_scenarios":len(group),
            "lambda_ML_range":[float(lam.min()),float(lam.max())],
            "lambda_ML_median":float(np.median(lam)),
            "n_lambda_ML_at_least_0_1":int(np.count_nonzero(lam>=.1)),
            "n_deltaAIC_favours_lambda_more_than_2":int(np.count_nonzero(lrt>4)),
            "lrt_range":[float(lrt.min()),float(lrt.max())],
            "lambda_values_not_independent_evolutionary_replicates":True,
        }
    summary={
      "status":"EMPIRICAL_S1_IMAGE_ANGLE_PHYLOGENETIC_SIGNAL_SENSITIVITY",
      "primary_response":"image_vertical_relative_continuous_orientation_not_gravitational_nodding",
      "original_native_S1_eligible_image_observations":170066,
      "original_native_S1_accepted_taxon_keys":244,
      "n_frozen_grafted_phylogenetic_scenarios":52,
      "source_194_includes_three_name_collisions":True,
      "analysis_cohorts":expected,
      "same_S1_five_climate_exposures":list(EXPOSURES),
      "not_comparable_to_independent_nuclear_species_tree":True,
      "frozen_original_S1_observation_model_reestimated":False,
      "nugget_scales_tested":list(nuggets),
      "bootstrap_N":bootstrap_reps,
      "phylogenetic_signal_does_not_validate_independent_nodding_origins":True,
      "result_panels":{},
      "first_archived_scenario_bootstrap_lambda0":bootstraps,
    }
    for cohort in expected:
        for nugget in nuggets:
            for phenotype in ("species_mean","species_median"):
                for design in ("intercept_only","S1_five_climate_adjusted"):
                    name=f"{cohort}|{phenotype}|{design}|nugget{nugget}"
                    group=[x for x in details if
                           x["cohort"]==cohort and x["nugget_fraction_of_tip_depth"]==nugget
                           and x["phenotype"]==phenotype and x["model"]==design]
                    if len(group)!=52:
                        raise ValueError("S1 52-tree scenario coverage dropped")
                    summary["result_panels"][name]=group_stats(group)
    return {"summary":summary,"scenario_diagnostics":details}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",required=True)
    ap.add_argument("--out",required=True)
    ap.add_argument("--bootstrap-reps",type=int,default=399)
    args=ap.parse_args()
    if args.bootstrap_reps<20 or args.bootstrap_reps>5000:
        raise ValueError("Bootstrap simulation budget outside bounds")
    species=verify_source(args.source)
    trees=restore_trees()
    result=analysis(species,trees,bootstrap_reps=args.bootstrap_reps)
    result["source_file_sha256"]=__import__("hashlib").sha256(Path(args.source).read_bytes()).hexdigest()
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    panels=result["summary"]["result_panels"]
    print("S1_PHYLOGENETIC_SIGNAL_SUMMARY="+json.dumps({
      "status":result["summary"]["status"],
      "source_cohorts":result["summary"]["analysis_cohorts"],
      "result_panels":{k:v for k,v in panels.items()
            if ("collision_free191" in k or "at_least10" in k)
            and "species_mean" in k and "nugget0.05" in k},
      "first_scenario_bootstrap":result["summary"]["first_archived_scenario_bootstrap_lambda0"],
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
