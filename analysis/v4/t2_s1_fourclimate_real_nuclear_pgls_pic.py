"""T2 continuation: BIO18/BIO05/RSDS/elevation from original photo sites on true nuclear trees.

Use the 2026-10-24 S1 FULL ORIGINAL PHOTO cohort empirically sampled
CHELSA BIO18+BIO05, original CHELSA RSDS pixels and WorldClim 30-second DEM.
Match original 251 name units to 42 Siniscalchi NA genuinely multilocus
nuclear tips and 12 EAzami Japan38 176-nuclear-locus tips. Four-axis PGLS and
Felsenstein standardized independent contrasts are DIFFERENT ESTIMANDS.
Neither counts morphological nodding U/D events or genuine historical
climate/nodding gains. Never use chloroplast, ITS-only or original-Q20 failed trees.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from Bio import Phylo
from scipy import stats
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as dna
from analysis.v4 import t2_s1_photo_nuclear_pgls_bio1_bio15_bridge as gls

S1_ORIGINAL_SHA=dna.S1_SHA
NEW_AGGREGATE_SHA="823ecb8442eb03598253cb33f3e16bac86538c6c2652f48ca2dc0886d0e19626"
COVARIATES=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
NEW_SOURCE_PHOTOS=169020
REGIONS=("NA42_genuine_published_nuclear","Japan12_genuine_176_loci_nuclear")
NAMES=251
N_PIC_CORRELATION_MIN=8

def new_taxon_source(path, original):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=NEW_AGGREGATE_SHA:
        raise ValueError("Original same-photo four climate aggregate source SHA changed")
    rr=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    if len(rr)!=NAMES or len({r["original_source_name"] for r in rr})!=NAMES:
        raise ValueError("S1 251 original source name four-axis source missing or duplicate")
    if sum(int(r["n_photo_complete_four_axis"]) for r in rr)!=NEW_SOURCE_PHOTOS:
        raise ValueError("S1 169020 four-axis measured photo counts changed")
    x={}
    for row in rr:
        name=row["original_source_name"]
        if name not in original or row["accepted_key"]!=original[name]["accepted_key"]:
            raise ValueError("New 4-axis original S1 species concept mismatched")
        if not 1<=int(row["n_photo_complete_four_axis"])<=int(original[name]["original_photo_count"]):
            raise ValueError("Original image counts cannot increase after new climate extraction")
        if not 0<=float(row["mean_image_angle_deg"])<=180:
            raise ValueError("Original S1 response is no longer 0-180 camera-relative angle")
        if not all(np.isfinite(float(row[k])) for k in COVARIATES):
            raise ValueError("New four-axis sources must be complete at taxon mean")
        if float(row["chelsa_bio18"])<0 or not -80<=float(row["chelsa_bio05"])<=70:
            raise ValueError("Incorrect original BIO18 precipitation or BIO05 temperature units")
        x[name]=row
    return x

def design(rows,fields=COVARIATES):
    values=np.array([[float(row[k]) for k in fields] for row in rows])
    sd=values.std(0)
    if not np.isfinite(values).all() or np.min(sd)<=1e-9:
        raise ValueError("Environmental S1 taxonomy subset climate degenerate")
    X=np.column_stack([np.ones(len(rows)),(values-values.mean(0))/sd])
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Cannot separate original source four climates under true nuclear panel")
    corr=np.corrcoef(values,rowvar=False)
    invcorr=np.linalg.inv(corr)
    return X,{"means":values.mean(0).tolist(),"SDs":sd.tolist(),
              "pairwise_pearson_r":corr.tolist(),
              "VIF":[float(v) for v in np.diag(invcorr)],
              "high_VIF_gt5_any":bool((np.diag(invcorr)>5).any())}

def independent_contrasts(tree,tipmap,agg):
    """Standardized Brownian contrasts, NOT independently observed U->D gains."""
    root=dna.prune(tree,tipmap)
    cols=["mean_image_angle_deg",*COVARIATES]
    pairs=[]
    def walk(node):
        if node.is_terminal():
            name=tipmap[node.name]
            return np.array([float(agg[name][x]) for x in cols],float),0.
        if len(node.clades)!=2:
            raise ValueError("True nuclear source subtree is not bifurcating after pruning")
        A,ca=walk(node.clades[0]);B,cb=walk(node.clades[1])
        va=float(node.clades[0].branch_length or 0)+ca
        vb=float(node.clades[1].branch_length or 0)+cb
        if min(va,vb)<=0:
            raise ValueError("Insufficient real nuclear branch lengths for Felsenstein contrasts")
        contr=(A-B)/np.sqrt(va+vb)
        pairs.append(contr)
        v=1/(1/va+1/vb)
        ancestor=(A/va+B/vb)*v
        return ancestor,float(v)
    walk(root.root)
    out=np.array(pairs,float)
    if out.shape!=(len(tipmap)-1,len(cols)):
        raise ValueError("Ancestral arithmetic did not produce n-1 independent contrasts")
    if not np.isfinite(out).all():
        raise ValueError("Original nuclear branchwise contrasts invalid")
    return out

def pic_report(contr):
    n=len(contr)
    if n<N_PIC_CORRELATION_MIN:
        return {"status":"NOT_EVALUABLE_INSUFFICIENT_GENUINE_NUCLEAR_BRANCHES","n":n}
    results={}
    for i,axis in enumerate(COVARIATES,1):
        if np.std(contr[:,i])<1e-12 or np.std(contr[:,0])<1e-12:
            raise ValueError("S1 original photo and environment contrast variance vanished")
        r,p=stats.pearsonr(contr[:,0],contr[:,i])
        results[axis]={"PIC_continuous_image_angle_environment_Pearson_r":float(r),
          "approx_PIC_pairwise_correlation_two_sided_p":float(p),
          "not_independent_true_historical_U_to_D_gain_events":True}
    return {"status":"FELSENSTEIN_REAL_NUCLEAR_SOURCE_CONTRASTS","n":n,"axes":results,
            "tip_photo_sampling_and_phylogeny_uncertainty_not_propagated":True}

def analyze_panel(tree,tipmap,agg,minimum_photos):
    eligible={tip:name for tip,name in tipmap.items()
              if int(agg[name]["n_photo_complete_four_axis"])>=minimum_photos}
    if len(eligible)<9:
        return {"status":"NOT_EVALUABLE_TOO_FEW_SOURCE_PHOTOS_NUCLEAR_TIPS",
                "n_original_source_taxa":len(eligible)}
    names,K=gls.brownian_kernel(tree,eligible)
    rows=[agg[name] for name in names]
    X,cor=design(rows)
    y=np.array([float(row["mean_image_angle_deg"]) for row in rows])
    lam,ll=gls.profile_ml_lambda(y,X,K)
    fits={str(k):gls.fit_gls(y,X,K,v) for k,v in
          (("IID_equal_species","IID"),("estimated_Pagel_lambda",lam),("fixed_Brownian",1.0))}
    contrasts=independent_contrasts(tree,eligible,agg)
    return {"status":"S1_REAL_FOUR_CLIMATE_MULTILOCUS_NUCLEAR_PGLS_PIC",
      "n_authentic_source_taxa":len(names),
      "n_original_measured_photo_observations":sum(int(row["n_photo_complete_four_axis"]) for row in rows),
      "model":"S1_image_head_angle_mean ~ z(BIO18)+z(BIO05)+z(RSDS)+z(elevation_m)",
      "panel_standardization":cor,
      "Pagel_lambda_MLE_residual":float(lam),
      "Pagel_lambda_profile_logL":float(ll),
      "PGLS":fits,
      "PIC":pic_report(contrasts),
      "angle_is_EXIF_image_vertical_not_gravitational_nodding":True,
      "current_present_climate_not_ancestral_transition_time_climate":True}

def run(original251,new_aggregate,outdir,na_raw,j_raw,j_map,j_manifest):
    original,panels,proof=dna.read_sources(original251,na_raw,j_raw,j_map,j_manifest)
    agg=new_taxon_source(new_aggregate,original)
    reports={}
    for tag,(tree,mapping) in panels.items():
        reports[tag]={"min_photos_1":analyze_panel(tree,mapping,agg,1),
                      "min_photos_10":analyze_panel(tree,mapping,agg,10)}
    primary=[]
    for tag in REGIONS:
        dat=reports[tag]["min_photos_1"]
        if dat["status"]!="S1_REAL_FOUR_CLIMATE_MULTILOCUS_NUCLEAR_PGLS_PIC":
            continue
        for idx,label in enumerate(COVARIATES):
            primary.append((tag,label,
              dat["PGLS"]["estimated_Pagel_lambda"]["t_approx_two_sided_p_conditional_on_covariance"][idx+1],
              dat["PIC"]["axes"][label]["approx_PIC_pairwise_correlation_two_sided_p"]))
    q_gls=gls.bh_adjust([x[2] for x in primary])
    q_pic=gls.bh_adjust([x[3] for x in primary])
    for (tag,label,p1,p2),q1,q2 in zip(primary,q_gls,q_pic):
        reports[tag]["min_photos_1"].setdefault("BH_q_genuine_2panel_fouraxis_PGLS",{})[label]=q1
        reports[tag]["min_photos_1"]["PIC"]["axes"][label]["BH_q_across_true_two_nuclear_panels"]=q2
    out={
      "schema":"azami_V4_BIO18_BIO05_RSDS_ELEVATION_SOURCE_TRUE_NUCLEAR_S1_PGLS_PIC_20261024_v1",
      "source_original_251_S1_image_data_sha256":S1_ORIGINAL_SHA,
      "source_169020_new_four_climate_photo_species_aggregate_sha256":NEW_AGGREGATE_SHA,
      "response":"SOURCE_EXIF_RELATIVE_S1_IMAGE_HEAD_ANGLE_0_180_CONTINUOUS",
      "original_S1_native_image_universe":170066,
      "new_four_environment_matched_original_photo_observations":NEW_SOURCE_PHOTOS,
      "new_four_axes":list(COVARIATES),
      "four_exposures_selected_by_user_after_earlier_results":True,
      "authentic_true_multilocus_nuclear_species_trees":2,
      "independent_geographic_phylogenies_NOT_one_joint_54_or_161_tree":True,
      "no_plastid_or_flora_U_D_or_failed_111_Q20_tree":True,
      "signed_continuous_angle_change_and_climatic_contrast_do_NOT_count_historical_nodding_gains":True,
      "fixed_nuclear_tree_and_unestimated_camera_measurement_errors":True,
      "proof":proof,"nuclear_panels":reports}
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    (outdir/"T2_REAL_NUCLEAR_4AXIS_S1_ANGLE_PGLS_PIC.json").write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(10,5.2))
        labels=[]
        for tag in REGIONS:
            dat=reports[tag]["min_photos_1"]
            if dat["status"]!="S1_REAL_FOUR_CLIMATE_MULTILOCUS_NUCLEAR_PGLS_PIC":continue
            fit=dat["PGLS"]["estimated_Pagel_lambda"]
            short="NA42" if tag.startswith("NA") else "Japan12"
            for i,name in enumerate(COVARIATES,1):
                row=len(labels)
                labels.append(short+" | "+name)
                beta=fit["coefficients_degrees_per_SD"][i]
                lo,hi=fit["conditional_t_95CI_degrees_per_SD"][i]
                ax.plot([lo,hi],[row,row],linewidth=1.5)
                ax.plot(beta,row,"o")
        ax.set_yticks(np.arange(len(labels)),labels)
        ax.axvline(0,linestyle="--",linewidth=.8)
        ax.invert_yaxis()
        ax.set_xlabel("True nuclear PGLS S1 camera-head angle coefficient (° per 1 SD)")
        ax.set_title("Original S1 BIO18 / BIO5 / RSDS / elevation in true nuclear trees\n"+
                     "95% intervals conditional on topology, estimated lambda and 5% tip noise")
        fig.tight_layout();fig.savefig(outdir/"T2_REAL_NUCLEAR_S1_4AXIS_PGLS.png",dpi=160)
        plt.close(fig)
    except ImportError:
        pass
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original251",required=True,type=Path)
    p.add_argument("--new-four-species",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    source=[dna.upstream(dna.NA_REPO,dna.NA_REF,dna.NA_TREE),
        dna.upstream(dna.J_REPO,dna.J_REF,dna.J_TREE),
        dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MAP),
        dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MANIFEST)]
    r=run(args.original251,args.new_four_species,args.out,*source)
    for tag,cohort in r["nuclear_panels"].items():
        d=cohort["min_photos_1"]
        print("T2_REAL_BIO18_BIO05_RSDS_ELEVATION_NUCLEAR="+json.dumps({
          "panel":tag,"n":d["n_authentic_source_taxa"],
          "lambda":d["Pagel_lambda_MLE_residual"],
          "PGLS_four_climate":d["PGLS"]["estimated_Pagel_lambda"]["coefficients_degrees_per_SD"],
          "PGLS_95CI":d["PGLS"]["estimated_Pagel_lambda"]["conditional_t_95CI_degrees_per_SD"],
          "PIC":{k:v["PIC_continuous_image_angle_environment_Pearson_r"]
                 for k,v in d["PIC"]["axes"].items()}
        },sort_keys=True))
if __name__=="__main__":main()
