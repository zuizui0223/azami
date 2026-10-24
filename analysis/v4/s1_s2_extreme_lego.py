"""V4 S1/S2 extreme-climate analysis + phenotypic LEGO falsification.

Unlike V3's 9-construct / 36-edge result, this V4 analysis has eight actual
image-observable constructs / 28 pairs. Test coupled vs reconfigurable visible
phenotype dimensions on source-matched native observations, then check
native/introduced S2 coverage for two focal species. Phenotypic reconfiguration
is neither independently inherited genetic parts nor adaptive function.

All input obs_id, photographer IDs and precise locations remain PRIVATE in
GitHub Actions; only aggregate JSON and pair matrices are exported.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

from analysis.v4.native_integration import FEATURES, PAIRS, association

BASE = ["chelsa_bio01","chelsa_bio04","chelsa_bio12",
        "chelsa_bio15","chelsa_rsds_mean"]
EXTRA = ["chelsa_bio05","chelsa_bio06","chelsa_bio13","chelsa_bio14",
         "chelsa_bio17","chelsa_bio18","topo_elevation"]
ALL_ENV = BASE + EXTRA
FEATURE_COLS = list(dict.fromkeys(c for features in FEATURES.values() for c in features))
MIN_PAIR_TAXA = 20
MIN_TAXON_N = 5
COHORT_ORIGINAL_S1 = 313350
COHORT_ORIGINAL_S2 = 297622
S2_NAMES = {"Cirsium arvense","Cirsium vulgare"}
S2_RENAMING = {
    "orientation_image_vertical_angle":["orientation_angle_degrees_median"],
    "corolla_lab_lightness":["corolla_lab_lightness_median","corolla_lightness_median"],
    "corolla_lab_chroma":["corolla_lab_chroma_median"],
    "corolla_hue_sin":["corolla_hue_sin_median"],
    "corolla_hue_cos":["corolla_hue_cos_median"],
    "capitulum_outline_aspect_ratio":["shape_aspect_ratio_median","capitulum_outline_aspect_ratio_median"],
    "capitulum_outline_circularity":["shape_circularity_median","capitulum_outline_circularity_median"],
    "capitulum_outline_solidity":["shape_solidity_median","capitulum_outline_solidity_median"],
    "capitulum_width_profile_cv":["shape_width_cv_median","capitulum_width_profile_cv_median"],
}
MODULES = {"angle":"presentation", "lightness":"colour",
           "chroma":"colour","hue":"colour",
           "aspect":"outline","circularity":"outline",
           "solidity":"outline","width_profile":"outline"}
SCENARIOS = {
  "annual_BIO12": {"exposure":"chelsa_bio12","thermal":"chelsa_bio01","elevation":False},
  "annual_BIO12_elevation": {"exposure":"chelsa_bio12","thermal":"chelsa_bio01","elevation":True},
  "warmest_quarter_BIO18": {"exposure":"chelsa_bio18","thermal":"chelsa_bio01","elevation":False},
  "warmest_quarter_BIO18_elevation": {"exposure":"chelsa_bio18","thermal":"chelsa_bio01","elevation":True},
  "driest_quarter_BIO17": {"exposure":"chelsa_bio17","thermal":"chelsa_bio01","elevation":True},
  "driest_month_BIO14": {"exposure":"chelsa_bio14","thermal":"chelsa_bio01","elevation":True},
  "wettest_month_BIO13": {"exposure":"chelsa_bio13","thermal":"chelsa_bio01","elevation":True},
  "heat_BIO05_warm_rain": {"exposure":"chelsa_bio18","thermal":"chelsa_bio05","elevation":True},
  "cold_BIO06_warm_rain": {"exposure":"chelsa_bio18","thermal":"chelsa_bio06","elevation":True},
}
MODULE_NAMES = list(FEATURES)
SEED = 20261024


def require(flag, code):
    if not flag:
        raise ValueError(code)


def read_tables(native_path, focal_path, external_path, measurements_path=None):
    native = pd.read_csv(native_path,low_memory=False,dtype={"obs_id":str})
    focal = pd.read_csv(focal_path,low_memory=False,dtype={"obs_id":str})
    env = pd.read_csv(external_path,low_memory=False,dtype={"obs_id":str})
    require(len(native)==COHORT_ORIGINAL_S1 and native.obs_id.is_unique,
            "WRONG_S1_ORIGINAL_SOURCE")
    require(len(focal)==COHORT_ORIGINAL_S2 and focal.obs_id.is_unique,
            "WRONG_S2_ORIGINAL_SOURCE")
    require(env.obs_id.is_unique and set(EXTRA)<=set(env),"UNVERIFIED_EXTERNAL_JOIN")
    require(set(focal.taxon_name)==S2_NAMES and
            set(focal.range_state)=={"native","introduced"},"S2_WCVP_IDENTITY")
    require(set(BASE)<=set(native) and
            set(FEATURE_COLS)<=set(native),"S1_V4_FEATURES_MISSING")
    require({"accepted_key","latitude","longitude"}<=set(native),"S1_AUTHORITY_OR_GEOGRAPHY_MISSING")
    require({"obs_id","latitude","longitude","user_id","range_state"}<=set(focal),
            "S2_PRIVATE_WCVP_SOURCE_INCOMPLETE")
    original = len(native),len(focal)
    native=native.merge(env[["obs_id",*EXTRA]],on="obs_id",how="left",validate="one_to_one")
    focal=focal.merge(env[["obs_id",*EXTRA]],on="obs_id",how="left",validate="one_to_one")
    require((len(native),len(focal))==original,"ENVIRONMENT_JOIN_CHANGED_GRAIN")
    if measurements_path is not None:
        measured=pd.read_csv(measurements_path,dtype={"obs_id":str},low_memory=False)
        require(measured.obs_id.is_unique,"ORIGINAL_S2_MEASUREMENT_DUPLICATE")
        cols=["obs_id"]
        matched={}
        for canonical, aliases in S2_RENAMING.items():
            available=[col for col in (canonical,*aliases) if col in measured]
            if available:
                matched[canonical]=available[0]
                cols.append(available[0])
        if len(cols)>1:
            m=measured[cols].rename(columns={value:key for key,value in matched.items()})
            focal=focal.merge(m,on="obs_id",how="left",validate="one_to_one")
        focal.attrs["measurement_columns_verified"]=matched
    return native,focal


def env_levels(frame, keys, columns):
    values=frame[columns].apply(pd.to_numeric,errors="coerce")
    mask=np.isfinite(values.to_numpy(float)).all(axis=1)
    data=values.loc[mask]
    group=keys.loc[mask].astype(str)
    # Species means computed on ALL observed environments before complete-row restriction.
    means=values.groupby(keys).transform("mean")
    within=(values-means).loc[mask]
    among=means.loc[mask].groupby(group).first()
    return {"raw":data,"among_equal_taxa":among,"within_taxon_demeaned":within}


def diagnostic_matrix(level, names):
    v=level[names].to_numpy(float)
    if len(v)<len(names)+3 or not np.isfinite(v).all():
        return {"status":"INSUFFICIENT_ENVIRONMENT_SUPPORT","n":len(v)}
    v=v-v.mean(axis=0)
    sd=v.std(axis=0,ddof=0)
    if np.any(sd<1e-10):
        return {"status":"CONSTANT_ENVIRONMENT_PREDICTOR","n":len(v)}
    z=v/sd
    c=np.corrcoef(z,rowvar=False)
    eigen=np.linalg.eigvalsh(c)
    if np.min(eigen)<1e-9:
        vif={name:None for name in names}
        status="NEAR_SINGULAR"
    else:
        vif=dict(zip(names,map(float,np.diag(np.linalg.inv(c)))))
        status="OK"
    pear=pd.DataFrame(c,index=names,columns=names)
    spearman=level[names].corr(method="spearman")
    return {"status":status,"n":len(v),"variables":names,
            "pearson":np.round(pear.to_numpy(float),7).tolist(),
            "spearman":np.round(spearman.to_numpy(float),7).tolist(),
            "vif":vif,
            "max_abs_pair_r":float(np.max(np.abs(c-np.eye(len(c))))) if len(c)>1 else 0.,
            "condition_index":None if eigen[0]<=0 else float(np.sqrt(eigen[-1]/eigen[0]))}


def all_diagnostics(frame, names):
    levels=env_levels(frame,frame.accepted_key,names)
    return {name:diagnostic_matrix(matrix,names) for name,matrix in levels.items()}


def max_vif(diagnostics,names):
    if any(d.get("status")!="OK" for d in diagnostics.values()):
        return float("inf")
    return max(d["vif"][v] for d in diagnostics.values() for v in names)


def outcome_blind_set(env_frame,scenario):
    """Select one hydric + one temperature exposure and RSDS, prune extras.

    Comparisons use SAME all-12-complete original observations. Rejected
    protected conflicts stop; never search trait p-values for a better set.
    """
    rain=scenario["exposure"]; temp=scenario["thermal"]
    protected=[rain,temp,"chelsa_rsds_mean"]
    extra=["chelsa_bio04","chelsa_bio15"]
    if scenario["elevation"]:
        extra.append("topo_elevation")
    selected=list(dict.fromkeys(protected+extra))
    removal_priority=["chelsa_bio04","chelsa_bio15","topo_elevation"]
    actions=[]
    for _ in range(len(selected)+2):
        levels=env_levels(env_frame,env_frame.accepted_key,selected)
        diagno={k:diagnostic_matrix(v,selected) for k,v in levels.items()}
        large_pairs=[]
        for name,d in diagno.items():
            if d.get("status")!="OK":
                continue
            corr=np.asarray(d["pearson"])
            for i,a in enumerate(selected):
                for j in range(i+1,len(selected)):
                    if abs(corr[i,j])>.8:
                        large_pairs.append((name,a,selected[j],float(corr[i,j])))
        implicated=set(x for _,l,r,_ in large_pairs for x in (l,r))
        removable=next((v for v in removal_priority if v in implicated and v in selected),None)
        if removable is None:
            bad=max_vif(diagno,selected)>5
            if not large_pairs and not bad:
                return {"status":"SELECTED","selected":selected,"history":actions,
                        "diagnostics":diagno}
            if bad:
                max_by={v:max((d["vif"] or {}).get(v,float("inf")) or float("inf")
                              for d in diagno.values()) for v in selected}
                rem=[v for v in removal_priority if v in selected]
                removable=max(rem,key=lambda v:max_by[v]) if rem else None
            if removable is None:
                return {"status":"UNRESOLVED_PROTECTED_CONFOUNDING",
                        "selected":selected,"history":actions,
                        "diagnostics":diagno,
                        "protected":protected}
        actions.append({"dropped":removable,
                        "reason":"pairwise_abs_r_gt_0.8" if large_pairs else "vif_gt_5",
                        "pairs_above_0.8":large_pairs[:18]})
        selected.remove(removable)
    raise AssertionError("VARIABLE_SELECTION_LOOP_NOT_TERMINATING")


def complete_lego_cohort(native):
    eligible=np.ones(len(native),dtype=bool)
    for column in FEATURE_COLS:
        eligible &= np.isfinite(pd.to_numeric(native[column],errors="coerce").to_numpy(float))
        flag=column+"__eligible"
        if flag in native:
            eligible &= native[flag].fillna(False).astype(bool).to_numpy()
        else:
            raise ValueError("MISSING_V4_FEATURE_ELIGIBILITY "+column)
    part=native.loc[eligible].copy()
    counts=part.groupby("accepted_key").size()
    part=part.loc[part.accepted_key.isin(counts[counts>=MIN_TAXON_N].index)].copy()
    require(part.accepted_key.nunique()>=MIN_PAIR_TAXA,"INSUFFICIENT_COMPLETE_V4_LEGO_TAXA")
    require(part.obs_id.is_unique,"DUPLICATE_S1_HEAD_OBSERVATION")
    return part


def standardize_traits(frame):
    out=frame[FEATURE_COLS].to_numpy(float).copy()
    means=out.mean(axis=0)
    scale=out.std(axis=0)
    h1=FEATURE_COLS.index("corolla_hue_sin")
    h2=FEATURE_COLS.index("corolla_hue_cos")
    hue_scale=float(np.sqrt((scale[h1]**2+scale[h2]**2)/2))
    scale[h1]=scale[h2]=hue_scale
    require(np.isfinite(scale).all() and (scale>1e-8).all(),"DEGENERATE_TRAIT_COORDINATE")
    return (out-means)/scale


def split_covariances(values,group):
    labels,inverse=np.unique(np.asarray(group),return_inverse=True)
    require(len(labels)>=MIN_PAIR_TAXA,"TOO_FEW_TAXA_FOR_LEGO")
    counts=np.bincount(inverse)
    n=len(values); d=values.shape[1]
    mean=np.zeros((len(labels),d))
    np.add.at(mean,inverse,values)
    mean/=counts[:,None]
    middle=mean-mean.mean(axis=0)
    among=middle.T@middle/len(labels)
    residual=values-mean[inverse]
    # One total weight per taxon, independent of number of source photographs.
    within=(residual/counts[inverse,None]).T@residual/len(labels)
    return within,among,labels,inverse,mean


def rv_rows(within,among):
    idx={name:[FEATURE_COLS.index(c) for c in members]
         for name,members in FEATURES.items()}
    rows=[]
    for left,right in PAIRS:
        w=association(within,idx[left],idx[right])
        a=association(among,idx[left],idx[right])
        rows.append({"left":left,"right":right,"left_module":MODULES[left],
                     "right_module":MODULES[right],
                     "within_rv":float(w[0]),"among_rv":float(a[0]),
                     "among_minus_within":float(a[0]-w[0]),
                     "within_signed_r":float(w[1]) if np.isfinite(w[1]) else None,
                     "among_signed_r":float(a[1]) if np.isfinite(a[1]) else None})
    return rows


def regression_design(frame,variables):
    x=frame[variables].to_numpy(float)
    group=frame.accepted_key.astype(str).to_numpy()
    labels,inverse=np.unique(group,return_inverse=True)
    counts=np.bincount(inverse)
    mean=np.zeros((len(labels),len(variables)))
    np.add.at(mean,inverse,x)
    mean/=counts[:,None]
    between=mean[inverse]
    within=x-between
    centre=x.mean(axis=0)
    sd=x.std(axis=0)
    require(np.all(sd>1e-9),"COLLINEAR_ZERO_SD")
    b=(between-centre)/sd
    w=within/sd
    design=np.column_stack([np.ones(len(x)),b,w])
    if np.linalg.matrix_rank(design)!=design.shape[1]:
        raise ValueError("COLLINEAR_MODEL_DESIGN")
    return design,group


def residualize(frame,values,variables):
    design,group=regression_design(frame,variables)
    # Equal species total weight; the same sample is used for raw vs adjusted.
    counts=frame.groupby("accepted_key").size()
    weights=frame.accepted_key.map(1/counts).to_numpy(float)
    root=np.sqrt(weights/weights.sum())
    beta=np.linalg.lstsq(design*root[:,None],values*root[:,None],rcond=None)[0]
    return values-design@beta,beta,group,design


def equal_n_stats(values,adjusted,group,replicates=500,seed=SEED):
    """Matched taxon counts; precompute centred features once (not every draw)."""
    require(replicates>=100,"NOT_ENOUGH_RESAMPLES")
    _, adjusted_a, labels, inv, means_adjusted = split_covariances(adjusted,group)
    _, raw_a, _, _, means_raw = split_covariances(values,group)
    rng=np.random.default_rng(seed)
    indices=[np.flatnonzero(inv==j) for j in range(len(labels))]
    centred={"raw":values-means_raw[inv], "adjusted":adjusted-means_adjusted[inv]}
    among={"raw":raw_a,"adjusted":adjusted_a}
    fixed_between={}
    for mode in ("raw","adjusted"):
        edges=rv_rows(among[mode],among[mode])
        fixed_between[mode]=np.array([e["among_rv"] for e in edges])
    draws={mode:[] for mode in ("raw","adjusted")}
    pair_counts={mode:[] for mode in ("raw","adjusted")}
    for _ in range(replicates):
        selected=np.array([rng.choice(ix) for ix in indices])
        for mode in draws:
            cov=np.cov(centred[mode][selected],rowvar=False,bias=True)
            edges=rv_rows(cov,among[mode])
            within=np.array([e["within_rv"] for e in edges])
            between=fixed_between[mode]
            draws[mode].append(float(np.median(between)-np.median(within)))
            pair_counts[mode].append(int((between>within).sum()))
    return {mode:{
        "median_differential_RV":float(np.median(draws[mode])),
        "resampling_q025":float(np.quantile(draws[mode],.025)),
        "resampling_q975":float(np.quantile(draws[mode],.975)),
        "fraction_difference_positive":float(np.mean(np.asarray(draws[mode])>0)),
        "median_pairs_among_stronger":float(np.median(pair_counts[mode])),
        "not_confidence_interval":True
    } for mode in draws}


def lego_module_summary(rows):
    module_pairs=[]
    for row in rows:
        module_pairs.append(
            {"pair":row["left"]+"/"+row["right"],
             "same_module":row["left_module"]==row["right_module"],
             "within":row["within_rv"],"among":row["among_rv"]})
    ans={}
    for scale in ("within","among"):
        same=np.array([r[scale] for r in module_pairs if r["same_module"]])
        cross=np.array([r[scale] for r in module_pairs if not r["same_module"]])
        ans[scale]={"n_within_module_pairs":len(same),
                    "n_cross_module_pairs":len(cross),
                    "median_same_module_RV":float(np.median(same)),
                    "median_cross_module_RV":float(np.median(cross)),
                    "cross_minus_within_module_RV":float(np.median(cross)-np.median(same)),
                    "not_independent_genetic_modules":True}
    return ans


def fit_environment_trait_gradients(frame,values,selected):
    # Separate observation- and taxon-structure effects; no discovery-driven
    # predictor inclusion or automatic merger with the old Holm family.
    design,groups=regression_design(frame,selected)
    counts=frame.groupby("accepted_key").size()
    weight=frame.accepted_key.map(1/counts).to_numpy(float)
    result={}
    for name,features in FEATURES.items():
        betas=[]
        raw_ps=[]
        for col in features:
            j=FEATURE_COLS.index(col)
            fit=sm.WLS(values[:,j],design,weights=weight).fit(
                cov_type="cluster",cov_kwds={"groups":groups},
                use_t=False)
            coefficients=[]
            for index,v in enumerate(selected):
                for suffix,pos in (("between",index+1),
                                   ("within",index+1+len(selected))):
                    coefficients.append({"environment":v,"scale":suffix,
                       "beta_std":float(fit.params[pos]),
                       "se_cluster":float(fit.bse[pos]),
                       "p_component":float(fit.pvalues[pos])})
                    if v==selected[0]:
                        raw_ps.append(float(fit.pvalues[pos]))
            betas.append({"trait_coordinate":col,"components":coefficients})
        result[name]={"n":len(frame),"n_taxa":len(counts),"joint_circular_hue":name=="hue",
                      "components":betas,
                      "p_values_are_unadjusted_exploratory":True}
    return result


def analyse_s1(native):
    trait_only=complete_lego_cohort(native)
    # Every annual/warm/dry comparison must use the exact same observations.
    # Retain the pre-climate trait cohort separately for transparent coverage.
    complete_climate=np.isfinite(trait_only[ALL_ENV].to_numpy(float)).all(axis=1)
    common=trait_only.loc[complete_climate].copy()
    counts=common.groupby("accepted_key").size()
    common=common.loc[common.accepted_key.isin(counts[counts>=MIN_TAXON_N].index)].copy()
    require(common.accepted_key.nunique()>=MIN_PAIR_TAXA,
            "INSUFFICIENT_JOINT_TRAIT_CLIMATE_COHORT")
    diagnostics=all_diagnostics(native,ALL_ENV)
    outcome_blind={}
    for name,scenario in SCENARIOS.items():
        outcome_blind[name]=outcome_blind_set(native[["accepted_key",*ALL_ENV]],scenario)
    vals=standardize_traits(common)
    # aza3 LEGO-inspired observable three-way configuration audit.
    # This is a phenotype test, never a functional genomic G123 inventory.
    from analysis.v4.s1_s2_lego_three_way import test_three_way
    raw_three_way=test_three_way(vals,common.accepted_key,499,SEED)
    raw_within,raw_among,*_=split_covariances(vals,common.accepted_key)
    original=rv_rows(raw_within,raw_among)
    from analysis.v4.s1_s2_lego_triple import analyze_triad
    raw_triad_scale=analyze_triad(vals, common.accepted_key.to_numpy(), replicates=400)
    model_results={}
    for name,choice in outcome_blind.items():
        if choice["status"]!="SELECTED":
            model_results[name]={"status":choice["status"],
                                 "selected":choice["selected"]};continue
        variables=choice["selected"]
        selected=common.dropna(subset=variables)
        if len(selected)!=len(common):
            # ALL original candidate columns were required for this single
            # matched V4 cohort. Abandon comparisons if coverage is inconsistent.
            raise ValueError("NEW_ENVIRONMENT_MODEL_DIFFERENT_PHENOTYPE_ROWS")
        try:
            residual, coefficients, groups, design=residualize(common,vals,variables)
            within,among,*_=split_covariances(residual,groups)
            edges=rv_rows(within,among)
            triad_scale=analyze_triad(residual,groups,replicates=400)
            nboot=equal_n_stats(vals,residual,groups,replicates=500)
            gradients=fit_environment_trait_gradients(common,vals,variables)
        except (ValueError,np.linalg.LinAlgError) as exc:
            model_results[name]={"status":"NOT_IDENTIFIABLE",
                                 "reason":str(exc)[:240]};continue
        model_results[name]={
          "status":"ESTIMATED_POSTHOC_EXPLORATORY",
          "predictors":variables,
          "sample_n":len(common),
          "taxa_n":int(common.accepted_key.nunique()),
          "v4_integration_raw":{"n_edges":len(original),
             "count_among_stronger":sum(r["among_minus_within"]>0 for r in original),
             "module_summary":lego_module_summary(original)},
          "v4_integration_adjusted":{
             "n_edges":len(edges),
             "count_among_stronger":sum(r["among_minus_within"]>0 for r in edges),
             "median_within_RV":float(np.median([v["within_rv"] for v in edges])),
             "median_among_RV":float(np.median([v["among_rv"] for v in edges])),
             "module_summary":lego_module_summary(edges)},
          "equal_n_sensitivity":nboot,
          "within_and_among_three_domain_configurations":triad_scale,
          "three_way_phenotype_after_environment":test_three_way(
              residual,groups,199,SEED+len(name)),
          "pairs":edges,
          "trait_environment_component_gradients":gradients,
          "warning":"No individual pair is a genetic module or an independent evolutionary transition; raw RV has sample-size/dimensionality bias."
        }
    return {
      "native_original_population":COHORT_ORIGINAL_S1,
      "all_eight_trait_before_climate_n":len(trait_only),
      "all_eight_trait_before_climate_taxa":int(trait_only.accepted_key.nunique()),
      "common_all_eight_trait_n":len(common),
      "common_all_eight_trait_taxa":int(common.accepted_key.nunique()),
      "common_complete12_environment_n":int(
          np.isfinite(common[ALL_ENV].to_numpy(float)).all(axis=1).sum()),
      "all_12_predictor_collinearity":diagnostics,
      "selection_by_scenario":outcome_blind,
      "raw_28_pairs":original,
      "raw_within_and_among_three_domain_configurations":raw_triad_scale,
      "raw_three_way_phenotype":raw_three_way,
      "scenario_results":model_results
    }


def aggregate_s2(focal):
    report={"source_wild_two_species_observations":len(focal),
            "climate_exposures":ALL_ENV,
            "original_measurement_column_availability":focal.attrs.get(
                "measurement_columns_verified",{}),
            "taxa":{}}
    for species,sp in focal.groupby("taxon_name"):
        output={}
        for state,g in sp.groupby("range_state"):
            variables={}
            for name in ALL_ENV:
                arr=g[name].to_numpy(float)
                arr=arr[np.isfinite(arr)]
                variables[name]={"n":len(arr),"median":float(np.median(arr)) if len(arr) else None,
                                 "q05":float(np.quantile(arr,.05)) if len(arr) else None,
                                 "q95":float(np.quantile(arr,.95)) if len(arr) else None}
            measured={col:int(pd.to_numeric(g[col],errors="coerce").notna().sum())
                      for col in FEATURE_COLS if col in g.columns}
            output[state]={"observations":len(g),
                           "observers":int(g.user_id.nunique()),
                           "climate":variables,
                           "trait_coverage":measured}
            available=[k for k in MODULE_NAMES
                       if all(c in g.columns for c in FEATURES[k])]
            paired=[]
            for i,left in enumerate(available):
                for right in available[i+1:]:
                    cs=FEATURES[left]+FEATURES[right]
                    dat=g.dropna(subset=cs)
                    if len(dat)<100 or dat.user_id.nunique()<10:
                        continue
                    x=dat[cs].to_numpy(float)
                    if np.min(x.std(axis=0))<1e-12:
                        continue
                    # Only independent original observation records, NOT genetic
                    # replicate plants. No categorical tag inferred from colour.
                    cov=np.cov(x,rowvar=False,bias=True)
                    nl=len(FEATURES[left]);nr=len(FEATURES[right])
                    lr=association(cov,list(range(nl)),list(range(nl,nl+nr)))
                    paired.append({"left":left,"right":right,
                                  "n":len(dat),"rv":float(lr[0])})
            output[state]["available_pairwise_visible_phenotypes"]=paired
        require(set(output)=={"native","introduced"},"S2_RANGE_MISSING")
        output["introduced_minus_native_climate_medians"]={
          v:(output["introduced"]["climate"][v]["median"]-
             output["native"]["climate"][v]["median"])
          if (output["introduced"]["climate"][v]["median"] is not None and
              output["native"]["climate"][v]["median"] is not None) else None
          for v in ALL_ENV}
        output["cannot_infer_time_change_or_heritable_genetic_modules"]=True
        report["taxa"][species]=output
    return report


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--native",type=Path,required=True)
    ap.add_argument("--focal",type=Path,required=True)
    ap.add_argument("--environment-ext",type=Path,required=True)
    ap.add_argument("--s2-measurements",type=Path)
    ap.add_argument("--out-public",type=Path,required=True)
    args=ap.parse_args()
    args.out_public.mkdir(parents=True,exist_ok=True)
    native,focal=read_tables(args.native,args.focal,args.environment_ext,
                             args.s2_measurements)
    if not set(BASE)<=set(focal):
        from analysis.v4.s1_s2_restore_original_five import recover_original_five
        focal,baseline_receipt=recover_original_five(
            focal, args.out_public.parent/"private"/"original_s2_five",
            Path("analysis/v4/environment_source_metadata.json"))
    else:
        baseline_receipt={"status":"FIVE_NATIVE_S2_SOURCE_CLIMATES_ALREADY_PRESENT"}
    source={"analysis_id":"azami_v4_S1_S2_extreme_climate_and_phenotypic_LEGO_20261024",
            "status":"SOURCE_VERIFIED_V4_EXTENSION_EXPLORATORY",
            "input_source":"V4 original encrypted source, WCVP-native/introduced, original-pixel rasters",
            "not_V3_reanalysis":True,
            "not_aza3_molecular_gene_or_functional_test":True,
            "S2_original_five_CHELSA_provenance":baseline_receipt,
            "S1":analyse_s1(native),"S2":aggregate_s2(focal),
            "limits":{"gravity_calibrated_orientation":False,
                      "true_phyllary_spines_and_stickiness":False,
                      "molecular_Lego_functional_resource_test":False,
                      "spatial_independence_or_adaptation_proven":False,
                      "phenotype_contrasts_not_independent_colonization_events":True,
                      "S1_and_S2_native_source_reused_not_independent":True,
                      "extra_climate_family_is_posthoc_20261024":True}}
    (args.out_public/"S1_S2_extreme_LEGO_results.json").write_text(
        json.dumps(source,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":source["status"],
                      "S1_native":source["S1"]["native_original_population"],
                      "S1_common_n":source["S1"]["common_all_eight_trait_n"],
                      "S1_common_taxa":source["S1"]["common_all_eight_trait_taxa"],
                      "S1_models_estimated":{k:v["status"] for k,v in
                           source["S1"]["scenario_results"].items()},
                      "S2_source":source["S2"]["source_wild_two_species_observations"]},
                     ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
