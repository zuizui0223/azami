"""T2: OWN genuine S1 photo head angle × real 108-candidate-locus / ITS gene trees.

Uses actual original source S1 251 names (244 accepted keys, 170066 measured
image-relative head angles), public Moreyra source 110 unique SRA BioSamples
with 108 SAME-original111-mask Q20 reference-copy nuclear candidate loci
(29310 positions), and new 109 exact-S1 original GenBank ITS gene-tree taxa.

CRITICAL comparison: 110 / 109 whole tree vs exactly SAME 58 original species
on both nuclear gene-tree reconstructions. The change in the BIO12 signal may
reflect taxon coverage rather than an ITS-vs-108-locus evolutionary mechanism.
PHOTO values are measured by S1 authors; flora states and alleged true
gravity-direction head nodding are NEVER used to classify or filter taxa.
Pagel's lambda is a correlated-tip phenomenological signal, not a BM/OU
evolutionary generator or adaptive radiation/trait-gain count.
Unrooted trees midpoint rooted only for conditional covariance calculation.
110-tip tree is posthoc technical source (111 primary gate FAILED). Multi-
locus copy orthology, homeology/reticulation, topology posterior unverified.
"""
from __future__ import annotations
import argparse
import csv
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy import linalg, optimize, stats
from Bio import Phylo

N_ORIGINAL_S1=251
N_ACCEPTED_WCVP=244
N_ORIGINAL_PHOTOS=170066
N_MULTILOCUS=110
N_ITS=109
N_EXACT_OVERLAP=58
N_Q20_LOCI=108
N_Q20_ORIGINAL_SITES=29310
EXPECTED_S1_SHA="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
EXPECTED_ML108_SHA="541f21801d79a7dda872a020835e6838ea0fd002668324a157d24bca16b9fa08"
EXPECTED_ITS_SHA="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"
ORIGINAL_EXCLUDED="Cirsium tioganum"
FIVE=("chelsa_bio01","chelsa_bio04","chelsa_bio12",
      "chelsa_bio15","chelsa_rsds_mean")
PHENOTYPES=("log_BIO12","S1_photo_mean_image_head_angle","S1_photo_median_image_head_angle")
NOISES=(0.0,0.05)
SEED=20261024


def sha256(source):
    return hashlib.sha256(Path(source).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as handle:
        return list(csv.DictReader(handle))


def verify_sources(original,source110,source_receipt,ML110,ITS109):
    if sha256(original)!=EXPECTED_S1_SHA or sha256(ML110)!=EXPECTED_ML108_SHA or sha256(ITS109)!=EXPECTED_ITS_SHA:
        raise ValueError("Original S1 photograph or true 110/109 reconstructed nucleotide tree SHA drift")
    master=read_csv(original)
    if (len(master)!=N_ORIGINAL_S1 or
        len({x["accepted_key"] for x in master})!=N_ACCEPTED_WCVP or
        len({x["original_source_name"] for x in master})!=N_ORIGINAL_S1 or
        sum(int(x["original_photo_count"]) for x in master)!=N_ORIGINAL_PHOTOS or
        any(x["botanical_gravity_posture_inferred_from_photo_angle"]!="False" for x in master)):
        raise ValueError("Original S1 251 photographic taxon concepts or 170066 images changed")
    byname={x["original_source_name"]:x for x in master}
    samples=read_csv(source110)
    rec=json.loads(Path(source_receipt).read_text())
    if (len(samples)!=N_MULTILOCUS or
        len({x["source_tip"] for x in samples})!=N_MULTILOCUS or
        len({x["source_BioSample"] for x in samples})!=N_MULTILOCUS or
        rec["posthoc_sensitivity_taxa"]!=110 or
        rec["original_all111_Q20_successful_real_SRA_BioSamples"]!=111 or
        rec["original_111_source_108_retained_loci"]!=N_Q20_LOCI or
        rec["observed_Q20_reference_positions_reused_unchanged"]!=N_Q20_ORIGINAL_SITES or
        rec["original111_technical_gate_passed"] is not False or
        rec["original_111_reference_position_selection_recomputed_on110"] is not False or
        rec["copy_number_paralogy_polyploidy_review_completed"] is not False):
        raise ValueError("Real 111 FAILED nuclear Q20 source gate must never become an accepted 110 species phylogeny")
    mapping={}
    for r in samples:
        key=r["taxon_name"].replace(" ","_")
        if (key not in byname or "WCVP"+byname[key]["accepted_key"]!=r["source_tip"]
            or r["independently_verified_one_copy_orthology"]!="False"
            or r["original_sample_selection_after_phenotype"]!="False"
            or key==ORIGINAL_EXCLUDED):
            raise ValueError("108 nuclear locus original S1 species or copy annotation not strictly sourced")
        mapping[r["source_tip"]]=key
    if len(set(mapping.values()))!=N_MULTILOCUS:
        raise ValueError("Two 110 public biological samples collapsed into one S1 source species")
    nu=Phylo.read(str(ML110),"newick")
    its=Phylo.read(str(ITS109),"newick")
    nuc_tips=[t.name for t in nu.get_terminals()]
    its_tips=[t.name for t in its.get_terminals()]
    if (len(nuc_tips)!=N_MULTILOCUS or set(nuc_tips)!=set(mapping) or
        len(its_tips)!=N_ITS or len(set(its_tips))!=N_ITS or
        not set(its_tips)<={r["original_source_name"] for r in master}):
        raise ValueError("Original two sequence-inferred trees do not match original S1 taxa")
    same=set(mapping.values())&set(its_tips)
    if len(same)!=N_EXACT_OVERLAP:
        raise ValueError("Original S1 source taxonomy no longer has exactly 58 same-species sequence tips")
    return byname,nu,its,mapping,same,rec


def original_true_covariance(tree,tip_to_source,*,keep=None):
    t=deepcopy(tree)
    if keep is not None:
        for tip in list(t.get_terminals()):
            if tip_to_source[tip.name] not in keep:
                t.prune(tip)
    t.root_at_midpoint()
    tips=t.get_terminals()
    original_names=[tip_to_source[x.name] for x in tips]
    depths=t.depths()
    tipdepth=np.array([float(depths[x]) for x in tips])
    if not len(set(original_names))==len(tips) or np.min(tipdepth)<=0:
        raise ValueError("True input gene tree has duplicated tips or zero root-to-tip length")
    med=float(np.median(tipdepth))
    M=np.empty((len(tips),len(tips)))
    for i,a in enumerate(tips):
        M[i,i]=tipdepth[i]/med
        for j in range(i):
            M[i,j]=M[j,i]=float(depths[t.common_ancestor(a,tips[j])])/med
    if not np.all(np.isfinite(M)) or np.min(np.linalg.eigvalsh(M))<=0:
        raise ValueError("True nuclear source midpoint Brownian kernel not PD")
    return original_names,M


def phenotype(name,byname,trait):
    if trait=="log_BIO12":
        return np.log(np.array([float(byname[n]["chelsa_bio12"]) for n in name]))
    field={"S1_photo_mean_image_head_angle":"mean_image_angle_deg",
           "S1_photo_median_image_head_angle":"median_image_angle_deg"}.get(trait)
    if not field:raise ValueError("No other externally-coded/fictional biological phenotype accepted")
    return np.array([float(byname[n][field]) for n in name])


def equal_species_five_rain_design(original_names,byname):
    raw=np.array([[float(byname[n][col]) for col in FIVE] for n in original_names])
    if not np.isfinite(raw).all() or np.any(raw[:,2]<=0):
        raise ValueError("Original annual CHELSA BIO12 source invalid")
    SD=raw.std(axis=0)
    if min(SD)<=0:raise ValueError("No five original environmental covariates estimable")
    X=(raw-raw.mean(axis=0))/SD
    X[:,2]=(raw[:,2]-raw[:,2].mean())/100.0
    return np.column_stack((np.ones(len(original_names)),X))


def slope(y,X,C):
    if len(y)<X.shape[1]+5 or X.shape[0]!=len(y):
        raise ValueError("Exact same original S1 sequence taxa needed for rain association")
    L=linalg.cholesky(C,lower=True,check_finite=True)
    x=linalg.solve_triangular(L,X,lower=True)
    z=linalg.solve_triangular(L,y,lower=True)
    coeff=np.linalg.lstsq(x,z,rcond=None)[0]
    res=z-x@coeff
    vcov=np.linalg.inv(x.T@x)*float(res@res)/(len(y)-X.shape[1])
    se=float(np.sqrt(vcov[3,3]))
    sd=stats.t.ppf(.975,len(y)-X.shape[1])
    return {"BIO12_original_S1_photo_angle_degrees_per100mm":float(coeff[3]),
            "conditional_95CI_degrees_per100mm":[float(coeff[3]-sd*se),
                                                 float(coeff[3]+sd*se)],
            "conditional_two_sided_t_p":float(2*stats.t.sf(abs(coeff[3]/se),len(y)-X.shape[1])),
            "same_original_S1_five_CHELSA_adjustment":list(FIVE)}


def pagel_spectral_kernel(C,noise):
    if not 0<=noise<=.2:raise ValueError("Photo technical nugget outside source-declared sensitivity")
    K=C+np.eye(len(C))*noise*np.median(np.diag(C))
    D=np.sqrt(np.diag(K))
    R=K/np.outer(D,D)
    np.fill_diagonal(R,1.0)
    eig,Q=np.linalg.eigh(R)
    if np.min(eig)<-1e-8:raise ValueError("Phylo covariance is invalid")
    return eig,Q,K


def profile_lambda(y,eig,Q,*,xatol=.0005):
    n=len(y)
    q=Q.T@np.ones(n)
    vals=Q.T@y
    def ll(lam):
        d=(1-lam)+lam*eig
        if np.min(d)<=1e-12:return float("-inf")
        weight=1/d
        mu=float(np.sum(weight*q*vals)/np.sum(weight*q*q))
        rss=float(np.sum(weight*(vals-mu*q)**2))
        return -.5*(n*(1+np.log(2*np.pi)+np.log(rss/n))+np.sum(np.log(d)))
    fit=optimize.minimize_scalar(lambda x:-ll(float(x)),bounds=(0,1),
                                method="bounded",options={"xatol":xatol})
    if not fit.success:
        raise ValueError("Pagel original nuclear correlated-tip optimizer did not converge")
    best=max([0.,1.,float(fit.x)],key=ll)
    return {"ML_lambda":float(best),"likelihood_ratio_over_unrelated_IID":float(
        max(0,2*(ll(best)-ll(0)))),
        "loglikelihood_IID":float(ll(0)),"loglikelihood_ML":float(ll(best))}


def permute_phylo_signal(y,eig,Q,n=999,seed=SEED):
    original=profile_lambda(y,eig,Q)
    obs=original["likelihood_ratio_over_unrelated_IID"]
    rng=np.random.default_rng(seed)
    more=0
    for _ in range(n):
        stat=profile_lambda(rng.permutation(y),eig,Q,xatol=.003)["likelihood_ratio_over_unrelated_IID"]
        more+=stat>=obs-1e-9
    return {"permutations":n,"as_or_more_extreme":int(more),
            "plus_one_source_tip_label_permutation_p":float((1+more)/(n+1)),
            "fixed_real_gene_tree_not_topology_posterior":True,
            "original_taxa_permutation_is_not_geographic_climate_null":True}


def cohort(t,name_map,byname,keep=None,*,permutations=999):
    names,M=original_true_covariance(t,name_map,keep=keep)
    results={}
    for noise in (0.0,0.05):
        eig,Q,C=pagel_spectral_kernel(M,noise)
        diag={}
        for trait in PHENOTYPES:
            y=phenotype(names,byname,trait)
            result=profile_lambda(y,eig,Q)
            if noise==0.05 and trait=="log_BIO12":
                result["permutation_diagnostic"]=permute_phylo_signal(
                    y,eig,Q,n=permutations,seed=SEED+len(names))
            diag[trait]=result
        results[("five_pct_source_tip_noise" if noise==.05 else "zero_tip_noise")]=diag
    X=equal_species_five_rain_design(names,byname)
    y=phenotype(names,byname,"S1_photo_mean_image_head_angle")
    noise=.05*np.median(np.diag(M))
    regression={
        "IID_equal_species_error":slope(y,X,np.eye(len(y))),
        "fixed_BM_phylogenetic_error_with_five_pct_noise":slope(
            y,X,M+noise*np.eye(len(y)))
    }
    return {
        "n_original_S1_taxa":len(names),
        "n_original_photo_angle_measurements":sum(int(byname[n]["original_photo_count"]) for n in names),
        "n_original_S1_taxa_in_exact_list":len(set(names)),
        "median_original_BIO12_mm":float(np.median([float(byname[n]["chelsa_bio12"]) for n in names])),
        "signal":results,
        "same_source_photo_angle_five_climate_rain_regression":regression,
        "photo_measurements_are_camera_vertical_not_true_gravity_head_tilt":True,
        "multilocus_phylogeny_paralogs_not_audited":True,
        "tree_root_midpoint_assumed_not_verified":True,
        "not_BM_OU_fitness_adaptive_evolution_or_U_to_D_gain_count":True,
    }


def run(original,source110,receipt,ML110,ITS109,*,permutations=999):
    byname,nu,its,mapping,same,rec=verify_sources(original,source110,receipt,ML110,ITS109)
    nu_map=mapping
    its_map={x.name:x.name for x in its.get_terminals()}
    groups={
        "108_candidate_locus_real110":cohort(nu,nu_map,byname,permutations=permutations),
        "nuclear_ITS_real109":cohort(its,its_map,byname,permutations=permutations),
        "108_candidate_locus_same58":cohort(nu,nu_map,byname,same,permutations=permutations),
        "nuclear_ITS_same58":cohort(its,its_map,byname,same,permutations=permutations),
    }
    if ([groups[x]["n_original_S1_taxa"] for x in groups]!=[110,109,58,58] or
        groups["108_candidate_locus_real110"]["n_original_photo_angle_measurements"]!=132762 or
        groups["nuclear_ITS_real109"]["n_original_photo_angle_measurements"]!=162146):
        raise ValueError("Sequence matched S1 photo measurement denominators changed")
    # Identical biological 58 original taxa should yield numerically IDENTICAL
    # independent-IID regression slopes even if the gene trees differ.
    regress=groups["108_candidate_locus_same58"]["same_source_photo_angle_five_climate_rain_regression"]["IID_equal_species_error"]
    reg_its=groups["nuclear_ITS_same58"]["same_source_photo_angle_five_climate_rain_regression"]["IID_equal_species_error"]
    if abs(regress["BIO12_original_S1_photo_angle_degrees_per100mm"]-reg_its["BIO12_original_S1_photo_angle_degrees_per100mm"])>1e-8:
        raise ValueError("Common 58 species did not preserve genuinely same climate and photographic phenotype")
    return {
        "schema":"azami_T2_true_110_108_multilocus_vs_109_ITS_same_58_S1_photo_angle_and_BIO12_v1",
        "status":"EXECUTED_TWO_DIFFERENT_AUTHENTIC_NUCLEAR_GENEALOGIES_PHOTO_CLIMATE_PHYLOGENY_COVERAGE_TEST",
        "S1_original_taxonomy_units":251,"S1_accepted_WCVP_concepts":244,
        "S1_full_original_photo_measurements":170066,
        "original_110_posthoc_source_gate_failed":True,
        "original108_nuclear_candidate_loci_and_Q20_positions":[108,29310],
        "original_same_58_biological_name_units":58,
        "same_58_source_name_list":sorted(same),
        "same_photo_outcome_and_same_climate_used_on_common58_IID_test_passed":True,
        "genuinely_DNA_inferred_gene_trees_not_accepted_true_species_tree":True,
        "both_original_ITS_and_Q20_tree_root_midpoint_for_BM":True,
        "trait_source_is_original_S1_OWN_photo_angle_not_external_flora_U_D":True,
        "BIO12_is_present_photo_observation_niche_not_historical_rain_at_origin":True,
        "photo_angular_measurement_error_fraction_is_assumption_not_calibration":True,
        "geographic_and_sequence_sample_selection_bias_not_controlled":True,
        "independent_phylogenetically_replicated_adaptive_morphological_gain_established":False,
        "source_sha256":{
            "original_S1_photographic_measurements":EXPECTED_S1_SHA,
            "original_NA_EA110_real_SRA_gene_candidate_ML":EXPECTED_ML108_SHA,
            "original_ITS109_ML":EXPECTED_ITS_SHA,
            "original110_SRA_BioSample_crosswalk":sha256(source110),
            "original111_Q20_failed_110_posthoc_receipt":sha256(receipt),
        },
        "results":groups,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--S1-original-251",type=Path,required=True)
    ap.add_argument("--nuclear-source110",type=Path,required=True)
    ap.add_argument("--original111-failed-posthoc110-receipt",type=Path,required=True)
    ap.add_argument("--real-108-locus-110-tree",type=Path,required=True)
    ap.add_argument("--GenBank-ITS109-tree",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--permutations",type=int,default=999)
    opt=ap.parse_args()
    if opt.permutations<99 or opt.permutations>9999:
        raise ValueError("Source tip-label permutation count outside declared range")
    result=run(opt.S1_original_251,opt.nuclear_source110,
               opt.original111_failed_posthoc110_receipt,
               opt.real_108_locus_110_tree,opt.GenBank_ITS109_tree,
               permutations=opt.permutations)
    opt.out.parent.mkdir(parents=True,exist_ok=True)
    opt.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,
                                  sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_OWN_S1_PHOTO_MULTILOCUS110_VS_ITS109_SAME58="+json.dumps({
        n:{"n":x["n_original_S1_taxa"],
           "logBIO12_lambda":x["signal"]["five_pct_source_tip_noise"]["log_BIO12"]["ML_lambda"],
           "photo_angle_lambda":x["signal"]["five_pct_source_tip_noise"]["S1_photo_mean_image_head_angle"]["ML_lambda"],
           "BIO12_phylo_permutation_p":x["signal"]["five_pct_source_tip_noise"]["log_BIO12"]["permutation_diagnostic"]["plus_one_source_tip_label_permutation_p"],
           "IID_angle_rain_beta":x["same_source_photo_angle_five_climate_rain_regression"]["IID_equal_species_error"]["BIO12_original_S1_photo_angle_degrees_per100mm"],
           "BM_angle_rain_beta":x["same_source_photo_angle_five_climate_rain_regression"]["fixed_BM_phylogenetic_error_with_five_pct_noise"]["BIO12_original_S1_photo_angle_degrees_per100mm"]}
        for n,x in result["results"].items()},
        ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__": main()
