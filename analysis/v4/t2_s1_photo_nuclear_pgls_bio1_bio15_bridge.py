"""Compare S1 continuous CAMERA-head angle with EAzami BIO1/BIO15 on true nuclear trees.

PGLS uses phylogenetically correlated residuals in angle ~ climate; this is
NOT the separate nearest-tip clustering test and NOT EAzami's discrete
flora-posture U->D transition/environment reconstruction. Keep the two true
nuclear species trees separate: NA published 63 pruned to 42 and Japan38
176-nuclear-locus tree pruned to 12. The S1 original 251 measurements and
climate fields are exactly source-SHA-pinned. No U/D, plastid, ITS109, or
posthoc failed nuclear 110 candidate tree is allowed as primary data.

This is a post-result cross-study translation, NOT a preregistered independent
validation of EAzami H1 or a causal mechanistic VPD/wetting test.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
from scipy import linalg, optimize, stats
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as src

FIVE=("chelsa_bio01","chelsa_bio04","chelsa_bio12",
      "chelsa_bio15","chelsa_rsds_mean")
TWO=("chelsa_bio01","chelsa_bio15")
DESIGNS={
    "two_axis_BIO1_BIO15":TWO,
    "two_axis_plus_BIO12":(TWO[0],TWO[1],"chelsa_bio12"),
    "frozen_S1_five_axes":FIVE,
    "EAzami_vector_BIO15up_BIO1down":("BIO15_minus_BIO1",),
    "EAzami_vector_plus_BIO12":("BIO15_minus_BIO1","chelsa_bio12"),
}
RESPONSES=("mean_image_angle_deg","median_image_angle_deg")
TIP_NOISE=.05
MIN_DF=5
MIN_TIPS=8

def brownian_kernel(tree,tipmap):
    t=src.prune(tree,tipmap)
    tips=t.get_terminals()
    names=[tipmap[x.name] for x in tips]
    depths=t.depths()
    K=np.empty((len(tips),len(tips)))
    for i,tip in enumerate(tips):
        K[i,i]=float(depths[tip])
        for j in range(i):
            K[i,j]=K[j,i]=float(depths[t.common_ancestor(tip,tips[j])])
    if np.min(np.diag(K))<=0:
        raise ValueError("Nuclear tree lacks positive original tip depth")
    K=K/float(np.median(np.diag(K)))
    K=(K+K.T)/2
    if np.linalg.eigvalsh(K).min()<=-1e-7:
        raise ValueError("Invalid nuclear Brownian shared branch covariance")
    return names,K

def design_matrix(rows,fields):
    original=np.array([[float(r[c]) for c in FIVE] for r in rows])
    mu=original.mean(0)
    sd=original.std(0)
    if np.min(sd)<=1e-12 or not np.isfinite(original).all():
        raise ValueError("Original 5 CHELSA species axes invalid")
    z=(original-mu)/sd
    raw={name:z[:,i] for i,name in enumerate(FIVE)}
    # The vector is the same *direction* as the frozen EAzami sign expectation,
    # not its original CTMC/Brownian state-transition estimator.
    combo=raw["chelsa_bio15"]-raw["chelsa_bio01"]
    if combo.std()<=0:raise ValueError("BIO15 - BIO1 composition is degenerate")
    raw["BIO15_minus_BIO1"]=combo/combo.std()
    X=np.column_stack([np.ones(len(rows))]+[raw[f] for f in fields])
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("Cannot identify requested original source design")
    return X,{
      "within_panel_pearson_BIO1_BIO15":
          float(np.corrcoef(raw["chelsa_bio01"],raw["chelsa_bio15"])[0,1]),
      "original_physical_variables":[*FIVE],
      "regression_fields":[*fields],
      "predictor_unit":"one standard deviation WITHIN selected source nuclear panel",
      "EAzami_composite_is_posthoc_translation_not_CTMC_U_to_D":True,
    }

def covariance(K,lam):
    if lam=="IID":
        return np.eye(K.shape[0])
    value=float(lam)
    if not 0<=value<=1:
        raise ValueError("Pagel lambda must be in the 0..1 interval")
    # Conventional Pagel off-diagonal transform with original unequal
    # Japan nuclear branch depths retained, fixed 5% independent tip error.
    V=np.diag(np.diag(K))+value*(K-np.diag(np.diag(K)))
    V+=TIP_NOISE*np.eye(K.shape[0])
    return (V+V.T)/2

def fit_gls(y,X,K,lam):
    n,p=X.shape
    if n<p+MIN_DF or K.shape!=(n,n) or not np.isfinite(y).all():
        raise ValueError("Insufficient matched taxa for 5-axis PGLS")
    V=covariance(K,lam)
    L=linalg.cholesky(V,lower=True,check_finite=True)
    wx=linalg.solve_triangular(L,X,lower=True)
    wy=linalg.solve_triangular(L,y,lower=True)
    if np.linalg.matrix_rank(wx)!=p:raise ValueError("Singular source climate + phylogeny regression")
    bet=np.linalg.lstsq(wx,wy,rcond=None)[0]
    error=wy-wx@bet
    rss=float(error@error)
    if rss<=1e-12:raise ValueError("Perfect source phenotype fit")
    sigma2=rss/(n-p)
    inv=np.linalg.inv(wx.T@wx)
    se=np.sqrt(np.diag(inv)*sigma2)
    critical=stats.t.ppf(.975,n-p)
    pp=2*stats.t.sf(np.abs(bet/se),n-p)
    ci=np.column_stack([bet-critical*se,bet+critical*se])
    ll=-.5*(n*(np.log(2*np.pi)+1+np.log(rss/n))+
            2*np.log(np.diag(L)).sum())
    return {
      "n_taxa":int(n),"residual_degrees_freedom":int(n-p),
      "coefficients_degrees_per_SD":bet.tolist(),
      "t_approx_two_sided_p_conditional_on_covariance":pp.tolist(),
      "conditional_t_95CI_degrees_per_SD":ci.tolist(),
      "maximized_gaussian_loglikelihood_fixed_covariance_shape":float(ll),
      "conditional_inference_does_not_propagate_lambda_or_tree_uncertainty":True,
    }

def profile_ml_lambda(y,X,K):
    likelihood=lambda x:fit_gls(y,X,K,float(x))[
        "maximized_gaussian_loglikelihood_fixed_covariance_shape"]
    grid=np.linspace(0,1,21)
    xs=[(float(t),likelihood(t)) for t in grid]
    i=int(np.argmax([t[1] for t in xs]))
    lower=float(grid[max(i-1,0)])
    upper=float(grid[min(i+1,20)])
    if lower<upper:
        sol=optimize.minimize_scalar(lambda x:-likelihood(x),
            bounds=(lower,upper),method="bounded",
            options={"xatol":1e-6})
        if sol.success:xs.append((float(sol.x),-float(sol.fun)))
    lam,ll=max(xs,key=lambda v:v[1])
    # Numerical boundary close to zero is still no measurable residual
    # phylogenetic correlation; do not invent a nonzero rate.
    lam=0. if lam<1e-5 else 1. if lam>1-1e-5 else lam
    return lam,ll

def analyse_tip_panel(tree,tipmap,by,*,photo_min):
    filt={tip:name for tip,name in tipmap.items()
          if int(by[name]["original_photo_count"])>=photo_min}
    if len(filt)<MIN_TIPS:
        return {"status":"NOT_EVALUABLE_TOO_FEW_SOURCE_NUCLEAR_TIPS",
                "n_species":len(filt),"required_n_species":MIN_TIPS}
    names,K=brownian_kernel(tree,filt)
    rows=[by[n] for n in names]
    result={
      "status":"SOURCE_MATCHED_TRUE_NUCLEAR_PGLS_COMPLETED",
      "n_taxa":len(names),
      "n_original_S1_photo_angles":sum(int(r["original_photo_count"]) for r in rows),
      "source_camera_angle_response_not_botanical_U_D":True,
      "model_tables":{},
    }
    for response in RESPONSES:
        y=np.array([float(r[response]) for r in rows])
        models={}
        for model,fields in DESIGNS.items():
            X,diagnostics=design_matrix(rows,fields)
            if len(y)-X.shape[1]<MIN_DF:
                models[model]={"status":"NOT_EVALUABLE_UNDER_PARAMETER_GATE","n_taxa":len(y)}
                continue
            lam,ll=profile_ml_lambda(y,X,K)
            fits={}
            for key,weight in (("IID_no_phylogeny","IID"),
                               ("fixed_Brownian_PGLS",1.0),
                               ("ML_Pagel_lambda_PGLS",lam)):
                fits[key]=fit_gls(y,X,K,weight)
            models[model]={
              "status":"evaluated",
              "response":response,
              "original_frozen_contrasts":list(fields),
              "source_diagnostics":diagnostics,
              "Pagel_lambda_ML_residual_from_this_design":float(lam),
              "fixed_tip_noise_fraction_of_median_nuclear_tree_depth":TIP_NOISE,
              "models":fits
            }
        result["model_tables"][response]=models
    return result

def bh_adjust(pvalues):
    p=np.array(pvalues,dtype=float)
    if p.size==0:return []
    if np.any(~np.isfinite(p)) or (p<0).any() or (p>1).any():
        raise ValueError("Cannot FDR-adjust invalid PGLS p-values")
    order=np.argsort(p);out=np.empty(len(p))
    ranked=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out[order]=np.minimum(1,ranked)
    return out.tolist()

def render_plot(summary,path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yvals=[];labels=[];pp=[];ci=[]
    labelscl={"chelsa_bio01":"BIO1 · annual temperature",
              "chelsa_bio15":"BIO15 · rain seasonality"}
    for pane,report in summary["nuclear_panels"].items():
        primary=report["min_photos_1"]["model_tables"]["mean_image_angle_deg"]
        for model in ("two_axis_BIO1_BIO15","frozen_S1_five_axes"):
            entry=primary[model]
            fit=entry["models"]["ML_Pagel_lambda_PGLS"]
            for j,field in enumerate(entry["original_frozen_contrasts"]):
                if field not in labelscl:continue
                yvals.append(len(yvals))
                label=("NA42" if pane.startswith("NA42") else "Japan12")
                labels.append(f"{label} | {('2-axis' if model.startswith('two') else '5-axis')} | {labelscl[field]}")
                pp.append(fit["coefficients_degrees_per_SD"][j+1])
                ci.append(fit["conditional_t_95CI_degrees_per_SD"][j+1])
    fig,ax=plt.subplots(figsize=(9,5.2))
    for i,(v,(low,high)) in enumerate(zip(pp,ci)):
        ax.plot([low,high],[i,i],color="0.4",lw=1.5)
        ax.plot(v,i,marker="o",color="C0" if i<4 else "C1")
    ax.axvline(0,color="black",lw=.75,ls="--")
    ax.set_yticks(yvals,labels)
    ax.invert_yaxis()
    ax.set_xlabel("PGLS angle coefficient (image degrees per within-panel 1 SD climate)")
    ax.set_title("Original S1 continuous photo angle | true multilocus nuclear trees\n"
                 "Conditional 95% intervals; family adjusted p-values in JSON")
    fig.tight_layout()
    fig.savefig(path,dpi=160)
    plt.close(fig)

def run(original251,*,na_raw,j_raw,j_map_raw,j_manifest_raw,output=None):
    by,panels,source=src.read_sources(original251,na_raw,j_raw,j_map_raw,j_manifest_raw)
    out={}
    for pane,(tree,tipmap) in panels.items():
        out[pane]={
          f"min_photos_{cut}":analyse_tip_panel(tree,tipmap,by,photo_min=cut)
          for cut in (1,10)
        }
    report={
      "schema":"V4_T2_S1_original_photo_continuous_true_nuclear_PGLS_EAzami_BIO1_BIO15_20261024",
      "status":"SOURCE_VERIFIED_TWO_SEPARATE_TRUE_MULTILOCUS_NUCLEAR_PGLS",
      "original_S1_all_251_names_and_170066_photos_not_refitted_as_mixed_model":True,
      "outcome":"S1_EXIF_IMAGE_VERTICAL_continuous_head_angle_degrees_NOT_gravitational_nodding",
      "authentic_nuclear_panels":2,"full_s1_original_species":251,
      "not_EAzami_U_D_transition_score_or_replication":True,
      "not_phylogenetically_grafted_or_plastid":True,
      "EAzami_frozen_H1_reference":{"discrete_source_flora_U_D_9taxa":"5U_4D",
        "original_BIO15_D_minus_U_SDs":[1.320,1.330],
        "original_BIO1_D_minus_U_SDs":[-0.975,-0.967],
        "original_H1_transition_composite_primary_12":"16/792",
        "original_H1_transition_composite_strict9":"4/126",
        "original_PGLS_BIO15_p":[.05054,.05239],
        "original_PGLS_BIO1_p":[.09604,.09793],
        "original_verdict":"unresolved_single_axis_phylogenetic_PGLS",
        "original_BIO15_BIO1_posthoc_conditional_geographic_models_not_independently_significant":True},
      "covariance_policy":{"Pagel_lambda_ML":"ML profile fitted within each mean/median angle~climate design",
        "lambda0_not_equal_IID_if_tree_tips_nonultrametric":True,
        "lambda1":"fixed Brownian shared-branch covariance",
        "independent_noise_fraction":TIP_NOISE,
        "conditional_student_t_intervals_not_full_phylogenetic_uncertainty":True},
      "source_nuclear_provenance":source,"nuclear_panels":out,
      "inferential_status":"posthoc_cross_study_hypothesis_translation_not_independent_confirmation",
    }
    for model,keys in (("two_axis_BIO1_BIO15",TWO),("frozen_S1_five_axes",FIVE)):
        pts=[]
        for pane,dat in out.items():
            entry=dat["min_photos_1"]["model_tables"]["mean_image_angle_deg"][model]
            fit=entry["models"]["ML_Pagel_lambda_PGLS"]
            for i,label in enumerate(keys):
                pts.append((pane,label,fit["t_approx_two_sided_p_conditional_on_covariance"][i+1]))
        adj=bh_adjust([z[2] for z in pts])
        for (pane,label,p),q in zip(pts,adj):
            out[pane]["min_photos_1"]["model_tables"]["mean_image_angle_deg"][model].setdefault(
                "primary_BH_FDR_across_both_regions_and_axes",{})[label]=float(q)
    if output is not None:
        output=Path(output);output.mkdir(parents=True,exist_ok=True)
        (output/"T2_REAL_NUCLEAR_S1_PHOTO_PGLS_BIO1_BIO15.json").write_text(
           json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
        render_plot(report,output/"T2_REAL_NUCLEAR_S1_ANGLE_BIO1_BIO15_PGLS.png")
    return report

def main():
    cli=argparse.ArgumentParser()
    cli.add_argument("--original251",required=True,type=Path)
    cli.add_argument("--outdir",required=True,type=Path)
    cli.add_argument("--NA-tree",type=Path)
    cli.add_argument("--Japan-tree",type=Path)
    cli.add_argument("--Japan-map",type=Path)
    cli.add_argument("--Japan-manifest",type=Path)
    args=cli.parse_args()
    original=(args.NA_tree,args.Japan_tree,args.Japan_map,args.Japan_manifest)
    if any(x is not None for x in original) and not all(x is not None for x in original):
        raise ValueError("Supply all immutable genuine nuclear sources, or none")
    raws=([x.read_bytes() for x in original] if all(x is not None for x in original)
         else [src.upstream(src.NA_REPO,src.NA_REF,src.NA_TREE),
               src.upstream(src.J_REPO,src.J_REF,src.J_TREE),
               src.upstream(src.J_REPO,src.J_REF,src.J_MAP),
               src.upstream(src.J_REPO,src.J_REF,src.J_MANIFEST)])
    result=run(args.original251,na_raw=raws[0],j_raw=raws[1],
       j_map_raw=raws[2],j_manifest_raw=raws[3],output=args.outdir)
    for pane,stats in result["nuclear_panels"].items():
        print("S1_NUCLEAR_PGLS_BIO1_BIO15="+json.dumps({
           "panel":pane,"n":stats["min_photos_1"]["n_taxa"],
           "two_axis":stats["min_photos_1"]["model_tables"]["mean_image_angle_deg"]["two_axis_BIO1_BIO15"]["models"]["ML_Pagel_lambda_PGLS"],
           "lambda":stats["min_photos_1"]["model_tables"]["mean_image_angle_deg"]["two_axis_BIO1_BIO15"]["Pagel_lambda_ML_residual_from_this_design"]
        },sort_keys=True))
if __name__=="__main__":main()
