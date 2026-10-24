"""Original S1 photographic orientation x newly reconstructed S1 GenBank ITS tree.

Use the OWNER'S measured continuous photograph-derived head-presentation
angle (mean and median over 170,066 original eligible image observations)
from the frozen S1 251 original-source-name registry, NOT any flora-derived
categorical upright/nodding label, and NOT the historical retired MobileNet
classifier that learned from species-description labels.

This runs an actual continuous comparative model on all 109 original-S1
taxon names with NCBI accession-verified ITS sequence tips. The predictor
set is exactly the original FIVE S1 environmental axes BIO1, BIO4, BIO12,
BIO15 and RSDS. Study outcome is 0..180 image-vertical angle, not a
gravity-calibrated tilt, flower peduncle bending, rainfall-protective fitness,
or a tree-informed inferred independent U->D evolutionary event.

The only nuclear marker in this standalone phylogeny is the multicopy ITS
region. Root is unverified; midpoint and raw Newick display roots are
separate documented statistical sensitivities, NEVER accepted ancestry.
"""
from __future__ import annotations
import argparse
import csv
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from Bio import Phylo
from scipy import linalg, optimize, stats

from analysis.v4.t2_na63_actual_tree_rain_gls import brownian_shared_ancestry

EXPECTED_MASTER=251
EXPECTED_WCVP=244
EXPECTED_PHOTOS=170066
EXPECTED_ITS_TIPS=109
EXPECTED_ITS_REF_SITES=645
EXPECTED_ITS_INFORMATIVE=103
EXPOSURES=("chelsa_bio01","chelsa_bio04","chelsa_bio12",
           "chelsa_bio15","chelsa_rsds_mean")
RESPONSES=("mean_image_angle_deg","median_image_angle_deg")
# Immutable original public S1 v4 source-name split, and new genuine GenBank
# accession-verified ITS ML results from independent GitHub Actions.
S1_ALL251_CSV_SHA256="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
ACCESSION_CSV_SHA256="1c778aa75b723f6af6f57f7f741ca0bccd65d3cbc14a6745531eb4596703f6f1"
ITS_ML_TREE_SHA256="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"
PHOTO_ARTIFACT=12421842342
ACCESSION_ARTIFACT=12483849333
ITS_TREE_ARTIFACT=12484282607
ROOT_SCENARIOS=("midpoint_ITS_unrooted","original_unrooted_Newick_display")
SUBCOHORTS=("full109","without_four_original_WCVP_concept_collisions")
SEED=20261023
# Exploratory measurement-error robustness from original S1 5%-tip-nugget
# convention, not an estimated photo/gravity calibration error.
BROWNIAN_TIP_NUGGET_FRAC=.05


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def load_evidence(s1path,accessionpath,treepath,receiptpath,preppath):
    if (sha256(s1path)!=S1_ALL251_CSV_SHA256 or
        sha256(accessionpath)!=ACCESSION_CSV_SHA256 or
        sha256(treepath)!=ITS_ML_TREE_SHA256):
        raise ValueError("Original S1 photos or independently queried nuclear ITS evidence SHA changed")
    master=rows(s1path);accessions=rows(accessionpath)
    if (len(master)!=EXPECTED_MASTER
        or len({r["accepted_key"] for r in master})!=EXPECTED_WCVP
        or len({(r["accepted_key"],r["original_source_name"]) for r in master})!=EXPECTED_MASTER
        or sum(int(r["original_photo_count"]) for r in master)!=EXPECTED_PHOTOS
        or sum(r["WCVP_key_collides_source_names"]=="True" for r in master)!=14):
        raise ValueError("Original own S1 photo cohort 244/251/170066 or 7 split keys changed")
    original={r["original_source_name"]:r for r in master}
    if (len(accessions)!=EXPECTED_ITS_TIPS or
        len({r["original_source_name"] for r in accessions})!=EXPECTED_ITS_TIPS):
        raise ValueError("Original 109 nuclear ITS species linked to S1 not present")
    receipt=json.loads(Path(receiptpath).read_text(encoding="utf-8"))
    prep=json.loads(Path(preppath).read_text(encoding="utf-8"))
    if (not receipt["all251_queries_completed_without_network_failures"]
        or receipt["source_names_actually_queried"]!=EXPECTED_MASTER
        or receipt["accepted_exact_species_nuclear_ITS_accessions"]!=EXPECTED_ITS_TIPS
        or receipt["output_accession_ledger_SHA256"]!=ACCESSION_CSV_SHA256
        or prep["sequence_quality"]["selected_real_ITS_alignment_columns"]!=EXPECTED_ITS_REF_SITES
        or prep["sequence_quality"]["selected_parsimony_informative_sites"]!=EXPECTED_ITS_INFORMATIVE
        or prep["sequence_quality"]["n_genbank_exact_species_tips"]!=EXPECTED_ITS_TIPS):
        raise ValueError("GenBank original 251 discovery and real-site ITS phylogeny no longer source-gated")
    keep=set()
    for r in accessions:
        name=r["original_source_name"].replace(" ","_")
        if (name not in original or name in keep or
            r["source_original_species"].replace(" ","_")!=name or
            r["GenBank_original_organism"].replace(" ","_")!=name or
            r["S1_species_to_GB_record_exact_name"]!="True" or
            r["ITS_phylogeny_single_copy_orthology_confirmed"]!="False" or
            r["accepted_key"]!=original[name]["accepted_key"] or
            r["original_photo_count"]!=original[name]["original_photo_count"]):
            raise ValueError("ITS source name is not the exact original S1 photographic taxon")
        keep.add(name)
    tree=Phylo.read(str(treepath),"newick")
    tips=[x.name for x in tree.get_terminals()]
    if (len(tips)!=EXPECTED_ITS_TIPS or len(set(tips))!=EXPECTED_ITS_TIPS
        or set(tips)!=keep or
        any(x.branch_length is None or not np.isfinite(x.branch_length) or x.branch_length<0
            for x in tree.find_clades() if x is not tree.root)):
        raise ValueError("IQ-TREE ML 109 source ITS tree does not match real GenBank accession tip IDs")
    if len({original[x]["accepted_key"] for x in tips})!=EXPECTED_ITS_TIPS:
        raise ValueError("Duplicated WCVP concept among genuine ITS tree tips")
    for r in master:
        if (r["botanical_gravity_posture_inferred_from_photo_angle"]!="False" or
            not 0 <= float(r["mean_image_angle_deg"]) <= 180 or
            not 0 <= float(r["median_image_angle_deg"]) <= 180 or
            any(not math.isfinite(float(r[k])) for k in EXPOSURES) or
            float(r["chelsa_bio12"])<=0):
            raise ValueError("Original S1 photo-angle/5 CHELSA sources are corrupted")
    return master,accessions,tree,receipt,prep


def source_design(master_by_tip,ordered):
    """Equal taxon weight, original S1 5-climate covariates, original photo angles."""
    n=len(ordered)
    if not 12<=n<=EXPECTED_ITS_TIPS:
        raise ValueError("S1 genus original molecular tree has too few real taxa")
    raw=np.array([[float(master_by_tip[name][v]) for v in EXPOSURES]
                  for name in ordered],dtype=float)
    if not np.isfinite(raw).all():
        raise ValueError("Non-finite source original CHELSA")
    sd=raw.std(axis=0)
    if (sd<=1e-12).any():
        raise ValueError("Original environmental axes are collinear or invariant")
    centered=(raw-raw.mean(axis=0))/sd
    # Keep rain in its original physical scale (photographic degrees / 100mm).
    centered[:,2]=(raw[:,2]-raw[:,2].mean())/100
    X=np.column_stack((np.ones(n),centered))
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Five fixed S1 climate predictors are not identifiable")
    return X, raw


def fit_gls(y,X,C,lam,*,structure):
    n,p=X.shape
    if n!=len(y) or C.shape!=(n,n) or not 0<=lam<=1:
        raise ValueError("Non-matching authentic nuclear source covariance or slope index")
    Sigma=(lam*C+(1-lam)*np.diag(np.diag(C)))
    L=linalg.cholesky(Sigma,lower=True,check_finite=True)
    x=linalg.solve_triangular(L,X,lower=True,check_finite=True)
    yy=linalg.solve_triangular(L,y,lower=True,check_finite=True)
    coef=np.linalg.lstsq(x,yy,rcond=None)[0]
    residual=yy-x@coef
    rss=float(residual@residual)
    if rss<=0:raise ValueError("Zero residual ecological angular variance")
    cov_beta=np.linalg.inv(x.T@x)*(rss/(n-p))
    se=np.sqrt(np.diag(cov_beta))
    j=3
    degrees=n-p
    ci=stats.t.ppf(.975,degrees)*se[j]
    logL=float(-.5*(n*np.log(2*np.pi)+2*np.log(np.diag(L)).sum()
                   +n*np.log(rss/n)+n))
    return {
        "original_nuclear_tip_count":n,
        "original_nuclear_ITS_only":True,
        "photo_derived_angle_is_dependent_variable":True,
        "source_yearly_rain_BIO12_beta_image_degrees_per100mm":float(coef[j]),
        "source_rain_beta_se_image_degrees_per100mm":float(se[j]),
        "BIO12_95pct_CI_image_degrees_per100mm":[float(coef[j]-ci),float(coef[j]+ci)],
        "BIO12_two_sided_t_p_conditional_on_fixed_lambda":float(
            2*stats.t.sf(abs(coef[j]/se[j]),degrees)),
        "BIO12_positive_direction_one_sided_t_p_conditional_on_fixed_lambda":float(
            stats.t.sf(coef[j]/se[j],degrees)),
        "loglikelihood":logL,
        "Pagel_lambda":float(lam),
        "residual_covariance":structure,
        "observational_model_n_not_independent_images":True,
        "lambda_estimation_uncertainty_not_in_Wald_CI":True,
    }


def fit_models(y,X,brownian):
    n=len(y)
    no_relation=fit_gls(y,X,np.eye(n),0,structure="IID_equal_species_weight")
    no_covariance=fit_gls(y,X,brownian,0,structure="Pagel_lambda0_nonultrametric_tip_specific_variance")
    bm=fit_gls(y,X,brownian,1,structure="Brownian_ITS_shared_ancestry")
    # Full Brownian on an ITS gene tree with ~1e-6 terminal edges can be
    # overconfident about images measured with camera-relative technical noise.
    # Keep the unregularized result, but separately audit exactly 5% diagonal
    # measurement nugget; do not tune the amount using rainfall association.
    C_meas=brownian + np.eye(n)*(
        BROWNIAN_TIP_NUGGET_FRAC*np.median(np.diag(brownian)))
    bm_measurement=fit_gls(
        y,X,C_meas,1,structure="Brownian_ITS_plus_5pct_image_measurement_nugget")
    def loss(val):
        try:return -fit_gls(y,X,brownian,float(val),structure="fitted_Pagel_lambda")["loglikelihood"]
        except (linalg.LinAlgError,ValueError):return float("inf")
    opt=optimize.minimize_scalar(loss,bounds=(0,1),method="bounded",
                                options={"xatol":1e-5})
    if not opt.success or not np.isfinite(opt.fun):
        raise ValueError("Cannot profile original nuclear ITS Pagel covariance")
    best=min((0.,float(opt.x),1.),key=loss)
    ml=fit_gls(y,X,brownian,best,structure="estimated_Pagel_lambda")
    return {
        "IID_equal_species_weight":no_relation,
        "Pagel_lambda0_diagonal_nonultrametric":no_covariance,
        "fixed_Brownian_lambda1":bm,
        "fixed_Brownian_plus_5pct_photo_measurement_nugget":bm_measurement,
        "estimated_Pagel_lambda":ml,
        "estimated_lambda_as_additional_parameter_in_AIC":float(
            2*(X.shape[1]+2)-2*ml["loglikelihood"]),
        "IID_AIC":float(2*(X.shape[1]+1)-2*no_relation["loglikelihood"]),
        "fixed_Brownian_AIC":float(2*(X.shape[1]+1)-2*bm["loglikelihood"]),
        "photos_do_not_measure_absolute_gravity_direction":True,
        "technical_tip_nugget_is_5pct_of_median_root_tip_depth":BROWNIAN_TIP_NUGGET_FRAC,
        "technical_nugget_was_not_tuned_to_observed_climate_signal":True,
    }


def evaluate_tree(tree,source,*,drop_original_WCVP_alias_collisions=False,
                  rooting="midpoint_ITS_unrooted"):
    t=deepcopy(tree)
    if rooting not in ROOT_SCENARIOS:
        raise ValueError("Not an auditable alternative original nuclear ITS root")
    if rooting=="midpoint_ITS_unrooted":
        t.root_at_midpoint()
    if drop_original_WCVP_alias_collisions:
        for leaf in list(t.get_terminals()):
            if source[leaf.name]["WCVP_key_collides_source_names"]=="True":
                t.prune(leaf)
    tips=[leaf.name for leaf in t.get_terminals()]
    n=len(tips)
    if (n!= (105 if drop_original_WCVP_alias_collisions else EXPECTED_ITS_TIPS)
        or len(set(tips))!=n):
        raise ValueError("Actual original 109/105 source subset not reproducible")
    C=brownian_shared_ancestry(t,tips)
    X,climate=source_design(source,tips)
    result={
        "root_strategy":rooting,
        "root_not_biologically_verified":True,
        "removed_any_original_WCVP_key_with_multiple_names":drop_original_WCVP_alias_collisions,
        "taxon_count":n,
        "CHELSA_five_original_variables":list(EXPOSURES),
        "different_tip_branch_depths_at_lambda0":True,
        "model_does_not_use_botanical_flora_U_D_or_species_label_classifier":True,
        "model_DOES_use_original_170066_image_derived_continuous_angle":True,
        "image_angle_0up_90horizontal_180down_is_camera_frame_not_gravity":True,
        "min_covariance_eigenvalue":float(np.linalg.eigvalsh(C)[0]),
        "angle_models":{},
    }
    for response in RESPONSES:
        y=np.array([float(source[tip][response]) for tip in tips])
        result["angle_models"][response]={
            "source_photo_angle_mean":float(y.mean()),
            "source_photo_angle_min_max":[float(y.min()),float(y.max())],
            "source_n_photo_observations":sum(int(source[tip]["original_photo_count"]) for tip in tips),
            "rain_gradient":fit_models(y,X,C)
        }
    return result


def reconstruct(s1path,accessionpath,treepath,receiptpath,preppath):
    master,acc,tree,receipt,prep=load_evidence(s1path,accessionpath,treepath,receiptpath,preppath)
    names={r["original_source_name"]:r for r in master}
    tips={x.name for x in tree.get_terminals()}
    groups={}
    for stage,drop in (("full109",False),("without_four_original_WCVP_concept_collisions",True)):
        for rooting in ROOT_SCENARIOS:
            groups[stage+"__"+rooting]=evaluate_tree(
                tree,names,drop_original_WCVP_alias_collisions=drop,rooting=rooting)
    matrix=[{
        "original_S1_source_name":row["original_source_name"],
        "accepted_WCVP_taxon_key":row["accepted_key"],
        "photographic_orientation_mean_deg":float(row["mean_image_angle_deg"]),
        "photographic_orientation_median_deg":float(row["median_image_angle_deg"]),
        "photographic_measurements":int(row["original_photo_count"]),
        "photo_axis_not_gravity":True,
        "exact_GenBank_nuclear_ITS_tree_tip":(
            row["original_source_name"] if row["original_source_name"] in tips else ""),
        "S1_original_name_trees_from_GenBank_accessions":row["original_source_name"] in tips,
        "source_name_WCVP_accepted_taxon_concept_conflict":row["WCVP_key_collides_source_names"]=="True",
        **{v:float(row[v]) for v in EXPOSURES},
    } for row in master]
    run={
        "schema":"azami_v4_T2_S1_photo_continuous_original5_BIO12_own109_GenBank_ITS_PGLS_v1",
        "status":"ACTUAL_SOURCE_S1_PHOTO_ANGLE_AND_INDEPENDENTLY_RECONSTRUCTED_GENBANK_ITS_PHYLOGENY",
        "biological_input_S1_photo_or_flora":"S1_owner_measured_continuous_photo_angle_ONLY",
        "previous_flora_binary_3D_36U_analyses_excluded_from_primary":True,
        "original_frozen_S1_170066_eligible_observations":EXPECTED_PHOTOS,
        "original_frozen_S1_WCVP_accepted_keys":EXPECTED_WCVP,
        "original_S1_source_name_operational_units":EXPECTED_MASTER,
        "real_GenBank_nuclear_ITS_one_original_S1_name_tip":EXPECTED_ITS_TIPS,
        "genuinely_new_ITS_phylogeny_from_MAFFT_IQTREE":True,
        "nuclear_ITS_alignment_original_sites":EXPECTED_ITS_REF_SITES,
        "nuclear_ITS_parsimony_informative_sites":EXPECTED_ITS_INFORMATIVE,
        "original_S1_five_covariates_unchanged":list(EXPOSURES),
        "environmental_BIO12_niche_is_current_photo_sampling_not_paleoclimate":True,
        "source_SHA256":{
            "S1_251_original_photo_angle_five_climate":sha256(s1path),
            "GenBank_exact_source_original109_accessions":sha256(accessionpath),
            "GenBank_new_ITS_109_ML_tree":sha256(treepath),
            "NCBI_original_251_discovery_receipt":sha256(receiptpath),
            "original_ITS_alignment_preparation":sha256(preppath),
        },
        "model_selection_not_confirmatory":True,
        "reconstructed_ITS_gene_tree_not_single_copy_nuclear_species_tree":True,
        "independent_evolutionary_origins_of_botanical_nodding_not_counted":True,
        "no_photo_orientation_to_gravity_U_D_conversion":True,
        "photo_sample_size_not_used_to_select_clades_or_weight_tips":True,
        "taxon_identity_alias_scope_only_source_names":True,
        "groups":groups,
    }
    if len(matrix)!=EXPECTED_MASTER or sum(x["S1_original_name_trees_from_GenBank_accessions"] for x in matrix)!=EXPECTED_ITS_TIPS:
        raise ValueError("Photo phenotype species cannot vanish because NCBI ITS unavailable")
    return run,matrix


def write_csv(path,rows):
    with Path(path).open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),extrasaction="raise")
        w.writeheader();w.writerows(rows)


def plot_phylogeny(tree,master,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    t=deepcopy(tree)
    t.root_at_midpoint()
    tips=t.get_terminals()
    depths=t.depths()
    positions={tip:float(i) for i,tip in enumerate(tips)}
    for node in t.find_clades(order="postorder"):
        if not node.is_terminal():
            positions[node]=sum(positions[c] for c in node.clades)/len(node.clades)
    height=max(float(depths[k]) for k in tips)
    if height<=0:raise ValueError("No genuine ITS branch lengths")
    fig,ax=plt.subplots(figsize=(17,29),dpi=150)
    pos=lambda n:4.2*float(depths[n])/height
    for n in t.find_clades(order="preorder"):
        x=pos(n)
        if not n.is_terminal():
            ys=[positions[c] for c in n.clades]
            ax.plot([x,x],[min(ys),max(ys)],c="#98a5ae",lw=.5,zorder=1)
        for child in n.clades:
            ax.plot([x,pos(child)],[positions[child],positions[child]],c="#98a5ae",lw=.5,zorder=1)
    photo_x=9.8;photo_w=2.7
    rain_x=14.35;rain_w=2.3
    mx=max(float(master[t.name]["chelsa_bio12"]) for t in tips)
    for leaf in tips:
        row=master[leaf.name];y=positions[leaf]
        angle=float(row["photographic_orientation_mean_deg"])
        rain=float(row["chelsa_bio12"])
        ax.text(4.4,y,leaf.name.replace("Cirsium_","C. "),fontsize=6.4,
                color="#253645",va="center")
        ax.barh(y,photo_w*angle/180,height=.55,left=photo_x,color="#577f98",
                alpha=.85)
        ax.text(photo_x+photo_w+.18,y,f"{angle:.1f}°",fontsize=6.3,va="center",color="#24536a")
        ax.barh(y,rain_w*rain/mx,height=.55,left=rain_x,color="#aa785c",alpha=.78)
        ax.text(rain_x+rain_w+.15,y,f"{rain:,.0f}",fontsize=6.3,va="center")
    ax.text(0,-1.3,"New 109-tip GenBank nuclear ITS tree × original S1 photo angle",
            fontsize=13,fontweight="bold")
    ax.text(photo_x,-.25,"S1 photo mean angle (image °)",fontsize=10,fontweight="bold")
    ax.text(rain_x,-.25,"Original S1 BIO12 (mm)",fontsize=10,fontweight="bold")
    ax.set_ylim(len(tips)+.5,-2)
    ax.set_xlim(-.25,rain_x+rain_w+1.35)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)
    fig.text(.08,.03,"109 source-matched ITS accession tips; S1 photo angle is EXIF-image-relative (NOT gravity/nodding). Nuclear ITS is multicopy and tree midpoint is arbitrary.",
             fontsize=8)
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(str(out),bbox_inches="tight",dpi=150)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s1-251",required=True,type=Path)
    p.add_argument("--its-accessions",required=True,type=Path)
    p.add_argument("--its-tree",required=True,type=Path)
    p.add_argument("--source-receipt",required=True,type=Path)
    p.add_argument("--alignment-receipt",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    p.add_argument("--no-plot",action="store_true")
    a=p.parse_args()
    report,rows=reconstruct(a.s1_251,a.its_accessions,a.its_tree,a.source_receipt,a.alignment_receipt)
    a.outdir.mkdir(parents=True,exist_ok=True)
    write_csv(a.outdir/"S1_ALL251_GENBANK_109_ACTUAL_PHOTO_ANGLE_CHELSA_FULL_ROSTER.csv",rows)
    path=a.outdir/"S1_PHOTO_109_NEW_GENBANK_ITS_PHYLO_BIO12_5_COVARIATE_RESULTS.json"
    path.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                    encoding="utf-8")
    if not a.no_plot:
        orig=Phylo.read(str(a.its_tree),"newick")
        lookup={r["original_S1_source_name"]:r for r in rows}
        plot_phylogeny(orig,lookup,a.outdir/"S1_PHOTO109_GENBANK_ITS_TREE_HEAD_ANGLE_AND_BIO12.png")
    x=report["groups"]["full109__midpoint_ITS_unrooted"]["angle_models"]["mean_image_angle_deg"]["rain_gradient"]
    print("S1_ACTUAL_PHOTO109_ORIGINAL_FIVE_CHELSA_NEW_GENBANK_ITS="+json.dumps({
        "n":109,
        "IID_BIO12_deg_per100mm":x["IID_equal_species_weight"]["source_yearly_rain_BIO12_beta_image_degrees_per100mm"],
        "Brownian_BIO12_deg_per100mm":x["fixed_Brownian_lambda1"]["source_yearly_rain_BIO12_beta_image_degrees_per100mm"],
        "estimated_lambda":x["estimated_Pagel_lambda"]["Pagel_lambda"],
        "model_selected_BIO12_deg_per100mm":x["estimated_Pagel_lambda"]["source_yearly_rain_BIO12_beta_image_degrees_per100mm"],
        "n_photo_observations_matched":report["groups"]["full109__midpoint_ITS_unrooted"]["angle_models"]["mean_image_angle_deg"]["source_n_photo_observations"],
        "photo_angles_are_gravity_nodding":False
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
