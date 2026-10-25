"""Source-type-aware, multivariate Brownian edge/niche tests for S1/S2.

PHOTO continuous variation is NOT a binary botanical CTMC transition.
Only independent source-validated categorical labels may enter a CTMC.
The method is the Luu conditional ancestral-mean operator used by EAzami;
with present niche centroids this is *not* historical transition-time climate.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np
from Bio import Phylo
from scipy.linalg import solve

FOUR = ("chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean", "elevation_m")
MIN_TIPS = 8


def selection_registry(registry: dict) -> dict:
    """Explicit measured/technical-only/missing distinction; no invented states."""
    endpoints = registry["endpoints"]
    if len(endpoints) != 27 or len({x["id"] for x in endpoints}) != 27:
        raise ValueError("The original S1/S2 27-endpoint registry drifted")
    measured = {x["id"] for x in endpoints if x["source_v2_status"] == "measured"}
    missing = {x["id"] for x in endpoints if x["source_v2_status"] != "measured"}
    if (len(measured), len(missing)) != (22, 5):
        raise ValueError("Original measured 22 / unavailable 5 status changed")
    groups = {}
    technical = set()
    for row in registry["dependent_measurement_groups"]:
        group = row["group"]
        members = [x["id"] for x in endpoints if x["mapping_group"] == group]
        if len(members) != row["n_endpoints"]:
            raise ValueError("Group size does not agree with 27-endpoint registry: " + group)
        if "technical_validation_only" in row["role"]:
            technical.update(members)
        else:
            active = tuple(x for x in members if x in measured)
            if active:
                groups[group] = active
    if len(technical & measured) != 4:
        raise ValueError("Original technical-only texture/specularity endpoints changed")
    if missing & technical:
        raise ValueError("Unmeasured outcomes cannot enter technical QC")
    inferential = measured - technical
    # Colour hue is a single paired circular endpoint in testing, never two tests.
    if not {"corolla_hue_sin", "corolla_hue_cos"} <= set(groups["flower_colour"]):
        raise ValueError("Circular hue pairing disappeared")
    return {"measured": sorted(measured), "unmeasured": sorted(missing),
            "technical_only": sorted(technical), "inferential": sorted(inferential),
            "groups": groups}


def tests_from_registry(reg: dict) -> dict[str, tuple[str, ...]]:
    """Separate endpoint and block estimands; do not independently test hue components."""
    units = {}
    for group, members in reg["groups"].items():
        if group == "flower_colour":
            units["endpoint:corolla_lab_lightness"] = ("corolla_lab_lightness",)
            units["endpoint:corolla_lab_chroma"] = ("corolla_lab_chroma",)
            units["block:corolla_hue_circular"] = ("corolla_hue_sin", "corolla_hue_cos")
            units["block:flower_colour_joint"] = members
        else:
            if len(members)>1:
                units["block:"+group] = members
            for v in members:
                units["endpoint:"+v] = (v,)
    if set(x for es in units.values() for x in es) != set(reg["inferential"]):
        raise ValueError("Observed and tested endpoint registry do not match")
    if "endpoint:corolla_hue_sin" in units or "endpoint:corolla_hue_cos" in units:
        raise ValueError("Circular hue cannot be tested as independent components")
    return units


def primary_trait_registry(inventory:dict)->dict:
    """Frozen S1/S2: three domains, eight constructs, nine coordinates.

    The 27-entry acquisition inventory is provenance, not the primary family.
    Bract, involucre and texture variables do not enter this analysis.
    """
    from analysis.v4.native_integration import FEATURES
    full=selection_registry(inventory)
    columns={c for values in FEATURES.values() for c in values}
    if len(FEATURES)!=8 or len(columns)!=9 or not columns<=set(full['measured']):
        raise ValueError('Frozen eight primary constructs or measurements changed')
    groups={k:full['groups'][k] for k in ('presentation','flower_colour','whole_head_shape')}
    if {c for values in groups.values() for c in values}!=columns:
        raise ValueError('Three-domain inventory differs from frozen S1/S2')
    return {'measured':sorted(columns),'inferential':sorted(columns),
            'groups':groups,'unmeasured':[],'technical_only':[]}


def primary_trait_tests(reg:dict)->dict[str,tuple[str,...]]:
    from analysis.v4.native_integration import FEATURES
    units={'construct:'+key:tuple(values) for key,values in FEATURES.items()}
    if {c for values in units.values() for c in values}!=set(reg['inferential']):
        raise ValueError('Primary trait scope mismatch')
    return units


def branch_operator(tree, tip_order: list[str]):
    """Return Luu-inverted Brownian child-minus-parent operator and edge masks."""
    if len(tip_order) < MIN_TIPS or len(tip_order) != len(set(tip_order)):
        raise ValueError("At least eight unique biologically matched taxa required")
    tips = {x.name: x for x in tree.get_terminals()}
    if set(tips) != set(tip_order):
        raise ValueError("Tree tip and original photographic source taxa do not match")
    nodes = list(tree.find_clades(order="preorder"))
    indices = {x: i for i, x in enumerate(nodes)}
    all_n = len(nodes)
    L = np.zeros((all_n, all_n))
    pairs = []
    for p in nodes:
        for c in p.clades:
            t = c.branch_length
            if t is None or not np.isfinite(t) or t < 0:
                raise ValueError("Invalid or missing source molecular branch length")
            w = 1.0 / max(float(t), 1e-10)
            i,j=indices[p],indices[c]
            L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
            pairs.append((i,j,not c.is_terminal()))
    fix = np.array([indices[tips[n]] for n in tip_order])
    free = np.array([i for i,n in enumerate(nodes) if not n.is_terminal()])
    if len(fix)+len(free)!=len(nodes):
        raise ValueError("Unexpected graph node classification")
    F = np.zeros((len(nodes),len(tip_order)))
    F[fix,np.arange(len(tip_order))]=1
    F[free]=solve(L[np.ix_(free,free)],-L[np.ix_(free,fix)],assume_a="pos")
    D = np.stack([F[j]-F[i] for i,j,_ in pairs])
    internal = np.array([b for _,_,b in pairs],dtype=bool)
    if not np.isfinite(D).all() or not np.allclose(D@np.ones(len(tip_order)),0,atol=1e-7):
        raise ValueError("Brownian operator failed constant-trait conservation")
    return D, internal


def standardize_matrix(matrix: np.ndarray) -> np.ndarray:
    x = np.asarray(matrix, dtype=float)
    if x.ndim==1:x=x[:,None]
    if x.ndim!=2 or x.shape[0]<MIN_TIPS or not np.isfinite(x).all():
        raise ValueError("Unreliable tip trait/environment matrix")
    sd=x.std(axis=0,ddof=1)
    if np.any(sd<1e-12):
        raise ValueError("One measured trait/environment component is constant")
    return (x-x.mean(axis=0))/sd


def _alignment(delta_t:np.ndarray,delta_e:np.ndarray,*,null_zero_is_zero:bool=False):
    """Signed scalar or norm-valued multivariate branchwise climate score."""
    den=np.linalg.norm(delta_t,axis=1).sum()
    if den < 1e-12:
        if null_zero_is_zero:
            # Exhaustive/MC null may place all internal ancestors at the same
            # value (e.g., balanced complementary tips); this is legitimate
            # zero alignment, not a numerical/selection failure.
            kind="signed_scalar" if delta_t.shape[1]==1 else "multivariate_nondirectional_norm"
            return np.zeros(delta_e.shape[1]),kind
        raise ValueError("No identifiable reconstructed trait change")
    matrix = delta_t.T @ delta_e / den
    if matrix.shape[0]==1:
        return matrix[0],"signed_scalar"
    return np.linalg.norm(matrix,axis=0),"multivariate_nondirectional_norm"


def _permutation_scores(D,E,T,nperm,seed):
    rng=np.random.default_rng(seed)
    denv=D @ standardize_matrix(E)
    temp=standardize_matrix(T)
    delta=D @ temp
    observed,kind=_alignment(delta,denv)
    ns=temp.shape[0]
    outcomes=np.empty((nperm,len(FOUR)),dtype=float)
    # Keep vectors paired (hue, colour and contour coordinates) under the null.
    for i in range(nperm):
        shuffled=temp[rng.permutation(ns)]
        outcomes[i],_= _alignment(D@shuffled,denv,null_zero_is_zero=True)
    return observed,outcomes,kind


def alignment_test(D,E,T,*,nperm:int=999,seed:int=20261024):
    """Tip-shuffle conditional rank, not biological replicate p-value."""
    if nperm<999:raise ValueError("At least 999 tip shuffles needed")
    if D.shape[1]!=len(E) or len(E)!=len(T) or np.shape(E)[1]!=len(FOUR):
        raise ValueError("Tree, trait and four environmental taxon cohorts differ")
    observed,null,kind=_permutation_scores(D,E,T,nperm,seed)
    values={}
    for i,climate in enumerate(FOUR):
        o=float(observed[i]); n=null[:,i]
        count=int(np.sum((np.abs(n)>=abs(o)-1e-12) if kind=="signed_scalar"
                         else (n>=o-1e-12)))
        values[climate]={"score":o, "conditional_tip_shuffle_fraction":(count+1)/(nperm+1),
                         "extreme_shuffles":count,"score_type":kind}
    return values


def co_shift_test(D,A,B,*,nperm:int=999,seed:int=20261024):
    """Cross-module alignment of magnitude, preserving entire module jointly."""
    if nperm<999:raise ValueError("At least 999 shuffles needed")
    a=standardize_matrix(A);b=standardize_matrix(B)
    da=np.linalg.norm(D@a,axis=1)
    db=np.linalg.norm(D@b,axis=1)
    def cosine(x,y):
        den=np.linalg.norm(x)*np.linalg.norm(y)
        return float((x@y)/den) if den>1e-12 else 0.0
    observed=cosine(da,db)
    rng=np.random.default_rng(seed)
    extremes=0
    for _ in range(nperm):
        if cosine(da,np.linalg.norm(D@b[rng.permutation(len(b))],axis=1))>=observed-1e-12:
            extremes+=1
    return {"co_shift_cosine":observed,"extreme_shuffles":extremes,
            "conditional_tip_shuffle_fraction":(extremes+1)/(nperm+1),
            "meaning":"shared reconstructed branchwise change magnitude, not shared genetic switches"}


def holm_family(pairs:list[tuple[str,float]]):
    """Strong-FWER correction, conservative for dependent endpoint/module tests."""
    order=sorted(pairs,key=lambda kv:kv[1])
    out={}; maxprev=0.0;m=len(order)
    for i,(name,p) in enumerate(order):
        maxprev=max(maxprev,(m-i)*float(p))
        out[name]=min(1.,maxprev)
    return out


def permutation_resolution(nperm:int, family_size:int, alpha:float=0.05):
    """Attainable first Holm threshold, not Monte Carlo precision assurance."""
    if (not isinstance(nperm,int) or nperm<1 or not isinstance(family_size,int)
            or family_size<1 or not 0<alpha<1):
        raise ValueError("Positive integer permutations/family size and valid alpha required")
    floor=family_size/(nperm+1)
    return {"permutations":nperm,"family_size":family_size,"alpha":alpha,
            "minimum_attainable_first_holm":min(1.,floor),
            "can_reject_at_alpha":floor<=alpha,
            "Monte_Carlo_precision_guaranteed":False}


def stable_seed(text:str)->int:
    return int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest()[:8],"big")
