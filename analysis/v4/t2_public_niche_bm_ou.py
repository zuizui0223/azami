"""Source-backed Cirsium climatic-niche BM versus one-optimum OU screen.

Public original sources, not a reconstruction of Moreyra et al. 2025:
- EAzami independent Comp1061 tree (one maximum likelihood + six
  AU-nonrejected alternative optimized branch-length topologies)
- EAzami FDT4 frozen original GBIF/thinned/CHELSA taxon means (Japan/Taiwan).

These are present-day distributional climatic niche centroids: not genetic
physiological tolerances or ancestors' environments. Tree branch lengths are
substitutions/site, not time. Conditional small-n AICc comparison alone cannot
show an adaptive peak, historical niche conservatism, selective regime shift,
ecological speciation, or global Cirsium radiation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
from io import StringIO
import json
from pathlib import Path

from Bio import Phylo
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

AXES=("chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15")
NUGGETS=(0.0,0.05,0.25)
EXPECTED_COUNTS={"japan7":7,"east_asia9":9,"east_asia12":12}
SOURCE_COUNTS={"Cirsium alpicola":23,"Cirsium brevicaule":35,
"Cirsium gyojanum":6,"Cirsium irumtiense":44,
"Cirsium japonicum var. australe":23,
"Cirsium japonicum var. japonicum":389,
"Cirsium kamtschaticum":81,"Cirsium kawakamii":14,
"Cirsium morii":9,"Cirsium suffultum":81,
"Cirsium tatakaense":9,"Cirsium yezoense":28}

def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def _tip_name(name):
    return name.replace(" ","_").replace(".","")

def read_rooted_trees(primary:Path,alternate:Path):
    out=[]
    raw=[("independent_Comp1061_ML",Path(primary).read_text())]
    for i,line in enumerate(Path(alternate).read_text().splitlines()[:6],start=1):
        if "]" not in line or not line.startswith("["):
            raise ValueError("Original source AU optimized tree annotation missing")
        raw.append((f"AU_nonrejected_{i}",line.split("]",1)[1].strip()))
    if len(raw)!=7:
        raise ValueError("Expected original ML and six AU trees")
    for label,newick in raw:
        tree=Phylo.read(StringIO(newick),"newick")
        if not any(x.name=="OUTGROUP_saff" for x in tree.get_terminals()):
            raise ValueError("Real EAzami Compositae outgroup absent")
        tree.root_with_outgroup({"name":"OUTGROUP_saff"})
        if len({t.name for t in tree.get_terminals()})!=len(tree.get_terminals()):
            raise ValueError("Duplicate nuclear tree tips")
        out.append((label,tree))
    return out

def shared_paths(tree,names):
    terminals={t.name:t for t in tree.get_terminals()}
    missing=set(names)-set(terminals)
    if missing:
        raise ValueError("Target nuclear tips missing: "+str(sorted(missing)))
    target=[terminals[name] for name in names]
    root=tree.common_ancestor(*target)
    depth=np.array([tree.distance(root,t) for t in target],float)
    if np.any(depth<=0):
        raise ValueError("Nonpositive tip depth in measured tree")
    n=len(names)
    mrca=np.zeros((n,n),float)
    for i,a in enumerate(target):
        for j,b in enumerate(target):
            if i==j:
                mrca[i,j]=depth[i]
            else:
                mrca[i,j]=tree.distance(root,tree.common_ancestor(a,b))
    if np.linalg.eigvalsh(mrca).min()<=1e-10:
        raise ValueError("Shared-path BM covariance not positive definite")
    return depth,mrca

def _covariance(depth,mrca,alpha_dimless,nugget):
    ref=np.median(depth)
    if alpha_dimless<=1e-8:
        v=mrca.copy()
    else:
        alpha=alpha_dimless/ref
        difference=depth[:,None]+depth[None,:]-2*mrca
        v=np.exp(-alpha*difference)*(-np.expm1(-2*alpha*mrca))/(2*alpha)
    return v/ref+np.eye(len(depth))*nugget

def profiled_gaussian_negative_ll(y,depth,mrca,alpha_dimless,nugget):
    v=_covariance(depth,mrca,alpha_dimless,nugget)
    try:
        chol=cho_factor(v,lower=True)
        n=len(y)
        u=np.ones(n)
        v1=cho_solve(chol,u)
        vy=cho_solve(chol,y)
        mu=float((u@vy)/(u@v1))
        resid=y-mu
        sigma2=float(resid@cho_solve(chol,resid)/n)
        if sigma2<=0 or not np.isfinite(sigma2):
            return np.inf
        logdet=2*np.log(np.diag(chol[0])).sum()
        return float(.5*(n*np.log(2*np.pi)+logdet+n*(1+np.log(sigma2))))
    except np.linalg.LinAlgError:
        return np.inf

def fit_BM_OU(y,depth,mrca,nugget,grid_points=150):
    y=np.asarray(y,float)
    if not np.isfinite(y).all() or len(y)<7 or np.std(y)<=0:
        raise ValueError("Unidentified small-n constant niche")
    y=(y-y.mean())/y.std(ddof=1)
    n=len(y)
    fun=lambda alpha:profiled_gaussian_negative_ll(y,depth,mrca,alpha,nugget)
    nll_bm=fun(0.)
    grid=np.r_[0.,np.geomspace(1e-5,250,grid_points)]
    vals=np.array([fun(a) for a in grid])
    best=[(float(vals[i]),float(grid[i])) for i in np.argsort(vals)[:8]]
    for i in np.argsort(vals)[:8]:
        if i<=1 or i>=len(grid)-1:
            continue
        opt=minimize_scalar(lambda la:fun(np.exp(la)),
            bounds=(np.log(grid[i-1]),np.log(grid[i+1])),
            method="bounded")
        if opt.success:
            best.append((float(opt.fun),float(np.exp(opt.x))))
    nll_ou,alpha=min(best,key=lambda t:t[0])
    if nll_ou>nll_bm+1e-7:
        raise ValueError("OU must be nested in BM when attraction alpha=0")
    def aicc(ll,k):
        return 2*ll+2*k+2*k*(k+1)/(n-k-1)
    difference=aicc(nll_ou,3)-aicc(nll_bm,2)
    return {
        "n_taxa":n,
        "tip_nugget_relative_to_median_root_tip":nugget,
        "BM_logLik":float(-nll_bm),
        "OU_logLik":float(-nll_ou),
        "BM_AICc":float(aicc(nll_bm,2)),
        "OU_AICc":float(aicc(nll_ou,3)),
        "OU_minus_BM_AICc":float(difference),
        "AICc_preferred":"OU" if difference<0 else "BM",
        "OU_alpha_times_median_root_tip_substitution_distance":float(alpha),
        "OU_alpha_at_grid_upper_bound":bool(alpha>=249),
    }

def analyse(centroids:Path,primary:Path,alternate:Path):
    data=pd.read_csv(centroids)
    expected={"taxon","n_observations",*AXES}
    if set(data)!=expected or data.taxon.duplicated().any():
        raise ValueError("Source-verified niche table does not match contract")
    counts=dict(zip(data.taxon,data.n_observations))
    if counts!=SOURCE_COUNTS:
        raise ValueError("Frozen EAzami JP/TW occurrence counts changed")
    if not np.isfinite(data[list(AXES)].to_numpy(float)).all():
        raise ValueError("Absent CHELSA source cells")
    trees=read_rooted_trees(primary,alternate)
    panels={}
    groupings={
        "japan7":[x for x in counts if x!="Cirsium japonicum var. australe" and
            x not in ("Cirsium kawakamii","Cirsium morii","Cirsium tatakaense")
            and counts[x]>=10],
        "east_asia9":[x for x in counts if counts[x]>=10],
        "east_asia12":[x for x in counts if counts[x]>=5],
    }
    for label,taxa in groupings.items():
        taxa=sorted(taxa)
        if len(taxa)!=EXPECTED_COUNTS[label]:
            raise ValueError("Original taxon treatment changed "+label)
        y=data.set_index("taxon").loc[taxa]
        out=[]
        for tree_name,tree in trees:
            lengths,shared=shared_paths(tree,[_tip_name(x) for x in taxa])
            for axis in AXES:
                for nugget in NUGGETS:
                    item=fit_BM_OU(y[axis].to_numpy(float),
                                   lengths,shared,nugget)
                    out.append({"tree_scenario":tree_name,"axis":axis,**item})
        panels[label]={
            "n_taxa":len(taxa),"taxon_names":taxa,
            "tree_scenarios":len(trees),
            "fits":out,
            "nAICc_preferred_OU":int(sum(x["AICc_preferred"]=="OU" for x in out)),
            "nAICc_preferred_BM":int(sum(x["AICc_preferred"]=="BM" for x in out)),
        }
    return {
        "schema_version":"azami_v4_public_EAzami_niche_BM_OU_20261022_v1",
        "status":"EXECUTED_PUBLIC_CIRSIUM_EAST_ASIA_BM_VS_OU_LIMITED_NICHE_SCREEN",
        "inputs":{
            "original_repository":"zuizui0223/EAzami",
            "source_CHELSA":"frozen FDT4 Japan/Taiwan thinned GBIF 1981-2010 climate observations",
            "original_Japan_actions_artifact":9516784077,
            "original_Taiwan_actions_artifact":9520157997,
            "taxon_CHELSA_centroid_file_sha256":_sha256(centroids),
            "primary_comp1061_tree_file_sha256":_sha256(primary),
            "nine_AU_candidate_tree_file_sha256":_sha256(alternate),
            "six_AU_not_independent_phylogenetic_replicates":True,
        },
        "model":"single-optimum stationary-target OU with fixed starting optimum, nested BM at alpha=0",
        "branch_units":"substitutions per site, NOT years; alpha reported in dimensionless median-root-to-tip scale",
        "selection":"AICc with kBM=2 and kOU=3; Gaussian maximum likelihood conditional on frozen source niches, tree and nugget",
        "trait_measurement":"species present-day climatic occupancy centroid, not measured physiological niche or plant morphology",
        "panels":panels,
        "claim_stop":{
            "not_all_public_Cirsium_phylogenies":True,
            "not_final_Moreyra_2025_350_locus_tree":True,
            "not_global_niche_radiation":True,
            "not_a_multi_peak_OU_regime_mapping":True,
            "not_dated_evolutionary_rates":True,
            "not_evidence_of_selection_or_adaptation":True,
            "tip_niche_sampling_uncertainty_unmodeled":True
        }
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--centroids",type=Path,required=True)
    a.add_argument("--primary-tree",type=Path,required=True)
    a.add_argument("--alternative-trees",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    result=analyse(args.centroids,args.primary_tree,args.alternative_trees)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_BM_OU_SUMMARY="+json.dumps({
        p:{"fits":len(x["fits"]),"n_taxa":x["n_taxa"],
            "BM_preferred":x["nAICc_preferred_BM"],
            "OU_preferred":x["nAICc_preferred_OU"]}
        for p,x in result["panels"].items()}))
if __name__=="__main__":
    main()
