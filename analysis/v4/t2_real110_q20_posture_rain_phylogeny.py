"""T2: own 110 original 108-locus Q20 nuclear ML × botanical head U/D × rain.

One real posthoc GTR+G4 molecular tree from 110 S1-matched public BioSamples,
108 120k-read nuclear reference-copy candidate loci, original 111-specimen
89/111 locus masking and 29,310 observed Q20 nucleotide coordinates.
The original111 gate FAILED due sparse C. tioganum (119 calls); this source
uses a disclosed posthoc110 technical sensitivity, not an accepted, curated
multi-copy-proof EAzami-type nuclear species tree.

DO NOT use camera-vertical photographic angle as gravity-rooted head U/D.
38/110 independently flora-coded botanically (D6 U32); 72 remain unknown.
All scored source names and current CHELSA BIO12/BIO15/BIO1 estimates are
checked at EXACT S1 original-name level with copy/sex/voucher uncertainty.

Root-free nearest-U vs D in original 110 topology + two arbitrary-root BM/
Pagel-lambda sensitivity on the same 38 measured taxa. Each missing posture
is unscored, not assumed U; no ancestral gain or rainfall adaptation inferred.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from copy import deepcopy
from pathlib import Path
import numpy as np
from Bio import Phylo
from scipy import linalg, optimize, stats

TREE_SHA="541f21801d79a7dda872a020835e6838ea0fd002668324a157d24bca16b9fa08"
SAMPLE_SHA="fa35f48c164f956aa131089156f5cb577a3fd0267fb7dce32f1bf1f8bc365f10"
AUTHORITY_SHA="8d127ccb93136d0ea8bcdc8f7a15c9bd712130e3de2d0d739209ec65b9488063"
SOURCE_RECEIPT_SHA="9046fc3a9bb68b8861ac6763c570210308740a46384cb91d00c503e7782f91ef"
SCOPED_N=110
D=6
U=32
UNKNOWN=72
EXCLUDED="Cirsium tioganum"

def read_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def fingerprint(path,digest):
    actual=hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if actual!=digest: raise ValueError("Actual public-author source SHA is not the pinned original: "+Path(path).name)
    return actual

def source_data(tree_path,samples_path,authority_path,receipt_path):
    for path,expected in ((tree_path,TREE_SHA),(samples_path,SAMPLE_SHA),
                          (authority_path,AUTHORITY_SHA),(receipt_path,SOURCE_RECEIPT_SHA)):
        fingerprint(path,expected)
    rec=json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    if (rec["posthoc_sensitivity_taxa"]!=110 or rec["original_111_source_108_retained_loci"]!=108
        or rec["observed_Q20_reference_positions_reused_unchanged"]!=29310
        or rec["original111_technical_gate_passed"] is not False
        or rec["original_111_reference_position_selection_recomputed_on110"] is not False
        or rec["true_Mk_nodding_transition_or_rainfall_adaptation_supported"] is not False):
        raise ValueError("Original failed 111 nuclear source provenance changed")
    samples=read_csv(samples_path)
    flora=read_csv(authority_path)
    if (len(samples)!=110 or len(flora)!=251 or
        len({z["source_BioSample"] for z in samples})!=110 or
        len({z["original_source_name"] for z in flora})!=251):
        raise ValueError("110 true original BioSamples / 251 original biological names not conserved")
    if any(z["independently_verified_one_copy_orthology"]!="False" for z in samples):
        raise ValueError("Original public candidate loci have never passed copy-orthology validation")
    bytax={z["original_source_name"]:z for z in flora}
    tree=Phylo.read(str(tree_path),"newick")
    tips=[z.name for z in tree.get_terminals()]
    bytip={z["source_tip"]:z for z in samples}
    if len(tips)!=110 or set(tips)!=set(bytip) or len(set(tips))!=110:
        raise ValueError("Genuine public original Q20 ML species labels drift")
    coded={}
    for tip in tips:
        lab=bytip[tip]
        name=lab["taxon_name"]
        if name not in bytax or not tip.startswith("WCVP") or tip!="WCVP"+lab["source_tip"].removeprefix("WCVP"):
            raise ValueError("One posthoc110 tip has no exact biological original sample")
        f=bytax[name]
        if int(tip[4:])!=int(f["accepted_key"]) or f["original_source_name"]!=name:
            raise ValueError("Actual nuclear BioSample WCVP does not match original botanical concept")
        if (f["botanical_state_not_exif_angle_threshold"]!="True"
            or f["independent_botanical_source_same_SRA_voucher_verified"]!="False"):
            raise ValueError("Only independent flora posture, not photographic camera angle or claimed same voucher")
        if f["botanical_U_D_state"] not in ("U","D","?"):
            raise ValueError("Invalid original floristic head posture")
        coded[tip]={**lab,**f}
    counts=Counter(z["botanical_U_D_state"] for z in coded.values())
    if counts!={"D":D,"U":U,"?":UNKNOWN}:
        raise ValueError("Original full110 independent botanical states changed")
    if any(float(z["chelsa_bio12"])<=0 for z in coded.values()):
        raise ValueError("Missing real annual rain source")
    return tree,coded,rec

def covariance(tree,tipnames):
    nodes={x.name:x for x in tree.get_terminals()}
    if set(nodes)!=set(tipnames):
        raise ValueError("True DNA source/phenotype taxon tips must match")
    depth=tree.depths()
    C=np.empty((len(tipnames),len(tipnames)))
    for i,a in enumerate(tipnames):
        for j,b in enumerate(tipnames[:i+1]):
            C[i,j]=C[j,i]=float(depth[tree.common_ancestor(nodes[a],nodes[b])])
    if np.min(np.linalg.eigvalsh(C))<=0:
        raise ValueError("Candidate nuclear tree BM shared history covariance singular")
    return C

def gaussian_gls(y,X,C,lam,kind):
    n,p=X.shape
    S=lam*C+(1-lam)*np.diag(np.diag(C))
    L=linalg.cholesky(S,lower=True,check_finite=True)
    xx=linalg.solve_triangular(L,X,lower=True)
    yy=linalg.solve_triangular(L,y,lower=True)
    beta=np.linalg.lstsq(xx,yy,rcond=None)[0]
    d=yy-xx@beta
    rss=float(d@d)
    if rss<=0: raise ValueError("No residual rain variation")
    var=np.linalg.inv(xx.T@xx)*rss/(n-p)
    se=np.sqrt(np.diag(var))
    likelihood=-.5*(n*(1+np.log(2*np.pi)+np.log(rss/n))
                    +2*np.log(np.diag(L)).sum())
    k=stats.t.ppf(.975,n-p)
    return {"coefficient_D_minus_U_mm":float(beta[1]),
            "95CI_D_minus_U_mm":[float(beta[1]-k*se[1]),float(beta[1]+k*se[1])],
            "two_sided_conditional_p":float(2*stats.t.sf(abs(beta[1]/se[1]),n-p)),
            "one_sided_D_wetter_conditional_p":float(stats.t.sf(beta[1]/se[1],n-p)),
            "se_mm":float(se[1]),"loglikelihood":float(likelihood),
            "phylogeny_statistical_covariance_kind":kind,
            "Pagel_lambda":float(lam),
            "conditional_Wald_CI_does_not_propagate_tree_or_lambda_uncertainty":True}

def fit_panel(tree,coded,exclude_names=frozenset()):
    t=deepcopy(tree)
    for tip in list(t.get_terminals()):
        if coded[tip.name]["botanical_U_D_state"]=="?" or coded[tip.name]["taxon_name"] in exclude_names:
            t.prune(tip)
    names=[x.name for x in t.get_terminals()]
    dat=[coded[x] for x in names]
    d=np.array([x["botanical_U_D_state"]=="D" for x in dat],float)
    rain=np.array([float(x["chelsa_bio12"]) for x in dat])
    temp=np.array([float(x["chelsa_bio01"]) for x in dat])
    season=np.array([float(x["chelsa_bio15"]) for x in dat])
    if int(d.sum())<5 or int((1-d).sum())<25: raise ValueError("Too few original florist-confirmed two-state nuclear tips")
    X0=np.column_stack((np.ones(len(d)),d))
    X=np.column_stack((X0,(temp-temp.mean())/temp.std(),
                       (season-season.mean())/season.std()))
    C=covariance(t,names)
    def fit(X):
        # IID true equal-variance comparison, not lambda=0 on non-ultrametric tree.
        iid=gaussian_gls(rain,X,np.eye(len(rain)),0,"IID_unrelated_equal_tip_variance")
        bm=gaussian_gls(rain,X,C,1,"Brownian_original108_candidate_nuclear_loci")
        diag=gaussian_gls(rain,X,C,0,"Pagel0_nonultrametric_tip_variance")
        out=optimize.minimize_scalar(lambda x:-gaussian_gls(rain,X,C,float(x),"profile")["loglikelihood"],
                                     bounds=(0,1),method="bounded")
        if not out.success:raise ValueError("Pagel nuclear PGLS likelihood profile failed")
        lam=min([0.,1.,float(out.x)],
                key=lambda z:-gaussian_gls(rain,X,C,z,"profile")["loglikelihood"])
        estimate=gaussian_gls(rain,X,C,lam,"estimated_Pagel_correlation_nuclear")
        return {"IID_equal_weight":iid,"Brownian_fixed_tree":bm,
                "lambda0_diagonal_not_IID":diag,"Pagel_MLE":estimate}
    return {
        "n_morphology_scored_source_tips":len(d),
        "D":int(d.sum()),"U":int((1-d).sum()),
        "D_minus_U_raw_mean_mm":float(rain[d==1].mean()-rain[d==0].mean()),
        "D_mean_mm":float(rain[d==1].mean()),
        "U_mean_mm":float(rain[d==0].mean()),
        "unadjusted":fit(X0),
        "BIO1_BIO15_adjusted":fit(X),
        "root_unverified_technical_nuclear_reference_tree":True,
        "taxon_concept_sra_same_voucher_not_verified":True,
    }

def nearest_U(tree,coded):
    Dtips=[t for t in tree.get_terminals()
           if coded[t.name]["botanical_U_D_state"]=="D"]
    Utips=[t for t in tree.get_terminals()
           if coded[t.name]["botanical_U_D_state"]=="U"]
    pairs=[]
    for x in Dtips:
        best=min(Utips,key=lambda y:(tree.distance(x,y),y.name))
        d=coded[x.name];u=coded[best.name]
        pairs.append({
            "original_nodding_source_species":d["taxon_name"],
            "nearest_flora_erect_source_species":u["taxon_name"],
            "D_minus_nearest_U_BIO12_mm":float(d["chelsa_bio12"])-float(u["chelsa_bio12"]),
            "patristic_nuclear_gene_tree_distance":float(tree.distance(x,best)),
            "matched_original_D_BIO12_mm":float(d["chelsa_bio12"]),
            "nearest_erect_may_be_reused":True,
            "nearest_species_is_not_independently_verified_voucher_pair":True,
        })
    if len(pairs)!=6:return None
    return {
        "pairs":sorted(pairs,key=lambda x:x["original_nodding_source_species"]),
        "n_nodding_with_nearest_erect_wetter":sum(x["D_minus_nearest_U_BIO12_mm"]>0 for x in pairs),
        "n_nodding_with_nearest_erect_drier":sum(x["D_minus_nearest_U_BIO12_mm"]<0 for x in pairs),
        "nearest_tip_comparisons_not_independent_gain_events":True,
    }

def execute(treepath,samplespath,authoritypath,receiptpath):
    tree,coded,source_receipt=source_data(treepath,samplespath,authoritypath,receiptpath)
    closest=nearest_U(tree,coded)
    rootings={}
    for kind in ("source_unrooted_Newick_display","midpoint_root_sensitivity"):
        t=deepcopy(tree)
        if kind=="midpoint_root_sensitivity":t.root_at_midpoint()
        main=fit_panel(t,coded)
        # The extreme photograph-niche BIO12 estimate is 4025mm for one
        # botanical D C. norikurense, based on only ONE S1 observation.
        # This post-result vulnerability is disclosed; no pretended prereg.
        exclude=fit_panel(t,coded,exclude_names={"Cirsium norikurense"})
        rootings[kind]={"all38_biologically_scored":main,
                        "posthoc_without_C_norikurense_single_site_climate":exclude}
    return {
        "schema":"azami_T2_real110_original108_nuclear_gene_tree_botanical_U_D_rain_history_v1",
        "status":"EXECUTED_AUTHENTIC_REAL110_MULTILOCUS_Q20_TECHNICAL_TREE_POSTURE_RAIN_SCOPE",
        "source_sha256":{
            "IQTREE_110_true_public_nuclear":TREE_SHA,
            "source110_original_SRA_BioSamples":SAMPLE_SHA,
            "independent_botanical_posture_BIO12_BIO1_BIO15":AUTHORITY_SHA,
            "original_111_failed_Q20_gate_and_110_posthoc_receipt":SOURCE_RECEIPT_SHA},
        "n_original_nuclear_source_tip":SCOPED_N,
        "n_original_nuclear_reference_candidate_loci":108,
        "n_real_original111_Q20_masked_nucleotide_sites":29310,
        "n_observed_D_original_flora":D,
        "n_observed_U_original_flora":U,
        "n_unresolved_botanical_phenotypes_kept_unknown":UNKNOWN,
        "original111_quality_gate_was_not_passed":True,
        "posthoc_source_quality_exclusion":EXCLUDED,
        "original_108_loci_single_copy_orthology_accepted":False,
        "author_nuclear_species_tree_validated":False,
        "root_independent_nearest_known_erect_nuclear_relatives":closest,
        "alternative_arbitrary_root_placement_sensitivities":rootings,
        "photo_CAMERA_angle_used_to_infer_D_U":False,
        "species_root_and_internal_branch_gain_events_inferred":False,
        "rainfall_selection_mediating_pollen_fitness_measured":False,
        "non_random_flora_coded_subset_38_may_bias_climate":True,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("tree","samples","botanical","receipt","out"):
        p.add_argument("--"+key,type=Path,required=True)
    q=p.parse_args()
    result=execute(q.tree,q.samples,q.botanical,q.receipt)
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,
                                indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_OWN_REAL110_SOURCE_NUCLEAR_POSTURE_RAIN="+json.dumps({
        "D":result["n_observed_D_original_flora"],
        "U":result["n_observed_U_original_flora"],
        "nearest_D_wetter":result["root_independent_nearest_known_erect_nuclear_relatives"]["n_nodding_with_nearest_erect_wetter"],
        "nearest_D_drier":result["root_independent_nearest_known_erect_nuclear_relatives"]["n_nodding_with_nearest_erect_drier"],
        "root_hypothesis_adjusted_BM":{k:r["all38_biologically_scored"]["BIO1_BIO15_adjusted"]["Brownian_fixed_tree"]["coefficient_D_minus_U_mm"]
               for k,r in result["alternative_arbitrary_root_placement_sensitivities"].items()},
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":main()
