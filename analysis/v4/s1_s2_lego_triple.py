"""V4 capitulum triple-configuration phenotypes, NOT genetic LEGO resources.

The anto/aza3 higher-order logic motivates direct measurement of joint
presentational, floral-colour and whole-head outline configurations. This
module quantifies an observed *phenotypic* 2x2x2 distribution relative to the
maximum-entropy distribution with exactly the same marginal and pairwise
probabilities. It cannot infer functional gene-resource co-retention, inherited
developmental modularity or repeated independent trait origins.

The triad is specified outcome-blind: image-relative presentation angle,
image Lab chroma and 2D outline aspect ratio. Thresholds are within-taxon
mean-centering or across-taxon median centering and are never optimized to
maximize a climate or fitness association. Cluster bootstrap is descriptive.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

TRIAD=("orientation_image_vertical_angle","corolla_lab_chroma",
       "capitulum_outline_aspect_ratio")
CODES=("low_low_low","low_low_high","low_high_low","low_high_high",
       "high_low_low","high_low_high","high_high_low","high_high_high")


def valid_probabilities(p):
    a=np.asarray(p,float)
    if a.shape != (8,) or not np.isfinite(a).all() or min(a)<-1e-12:
        raise ValueError("INVALID_THREE_MODULE_CONFIGURATION_PROBABILITIES")
    if abs(a.sum()-1)>1e-8:
        raise ValueError("THREE_MODULE_PROBABILITIES_DONT_SUM_TO_ONE")
    return np.maximum(a,0)/sum(np.maximum(a,0))


def frequencies_from_pairwise(t, singles, pairs):
    a,b,c=singles
    ab,ac,bc=pairs
    return np.array([
        1-a-b-c+ab+ac+bc-t, # 000
        c-ac-bc+t,          # 001
        b-ab-bc+t,          # 010
        bc-t,               # 011
        a-ab-ac+t,          # 100
        ac-t,               # 101
        ab-t,               # 110
        t,                  # 111
    ],float)


def maxent_from_observed(cells):
    p=valid_probabilities(cells)
    bits=np.array([[(i>>2)&1,(i>>1)&1,i&1] for i in range(8)],int)
    singles=np.sum(p[:,None]*bits,axis=0).tolist()
    pairs=[float(np.sum(p*(bits[:,i]*bits[:,j]))) for i,j in ((0,1),(0,2),(1,2))]
    a,b,c=singles
    ab,ac,bc=pairs
    lower=max(0.,ab+ac-a,ab+bc-b,ac+bc-c)
    upper=min(ab,ac,bc,1-a-b-c+ab+ac+bc)
    lower=max(0.,min(1.,lower))
    upper=max(0.,min(1.,upper))
    if lower>upper+1e-9:
        raise ValueError("IMPOSSIBLE_THREE_WAY_BOUNDS")
    p111=float(p[7])
    if p111<lower-1e-8 or p111>upper+1e-8:
        raise ValueError("OBSERVED_TRIPLE_OUTSIDE_PAIRWISE_MARGINAL_BOUNDS")
    def entropy(t):
        v=frequencies_from_pairwise(t,singles,pairs)
        q=np.clip(v,0,1)
        positive=q[q>0]
        return -float(np.sum(positive*np.log(positive)))
    if upper-lower<1e-9:
        expected=(lower+upper)/2
        status="PAIRWISE_DETERMINES_TRIPLE"
    else:
        solution=minimize_scalar(lambda t:-entropy(t),bounds=(lower,upper),
                                 method="bounded",options={"xatol":1e-12})
        if not solution.success or not np.isfinite(solution.x):
            raise ValueError("MAXENT_PAIRWISE_REFERENCE_DID_NOT_CONVERGE")
        expected=float(solution.x)
        status="TRIPLE_NOT_DETERMINED_BY_PAIRWISE_MARGINALS"
    return {
        "status":status,
        "p111_observed":p111,
        "p111_maximum_entropy_pairwise":expected,
        "p111_pairwise_feasible_low":lower,
        "p111_pairwise_feasible_high":upper,
        "observed_minus_pairwise_maxent":p111-expected,
        "fraction_of_feasible_triple_range":None if upper-lower<1e-9
                                            else float((p111-lower)/(upper-lower)),
        "p111_not_genetic_joint_resource_frequency":True,
        "all_8_observed_state_probabilities":list(map(float,p)),
        "all_8_maxent_reference_state_probabilities":list(map(
            float,np.maximum(0,frequencies_from_pairwise(expected,singles,pairs)))),
        "single_module_positive_probabilities":singles,
        "pairwise_joint_positive_probabilities":pairs,
        "n_nonempty_triple_configurations":int(np.sum(p>0)),
    }


def states_from_three_values(v):
    a=np.asarray(v,float)
    if a.ndim!=2 or a.shape[1]!=3 or not np.isfinite(a).all():
        raise ValueError("THREE_MODULE_TRAITS_REQUIRED")
    return (a[:,0]>0).astype(int)*4+(a[:,1]>0).astype(int)*2+(a[:,2]>0).astype(int)


def source_states(values,group,columns):
    from analysis.v4.s1_s2_extreme_lego import MIN_PAIR_TAXA,FEATURE_COLS
    if tuple(columns)!=TRIAD:
        raise ValueError("TRIAD_DEFINITION_DRIFT")
    indices=[FEATURE_COLS.index(k) for k in TRIAD]
    arr=np.asarray(values,float)[:,indices]
    labels,inv=np.unique(np.asarray(group),return_inverse=True)
    if len(labels)<MIN_PAIR_TAXA or len(arr)<len(labels)*5:
        raise ValueError("INSUFFICIENT_TAXON_SUPPORT_FOR_TRIAD")
    counts=np.bincount(inv)
    means=np.zeros((len(labels),3))
    np.add.at(means,inv,arr)
    means=means/counts[:,None]
    states_within=states_from_three_values(arr-means[inv])
    # The matching expected population is one species, not one photograph.
    state_species=np.zeros((len(labels),8),float)
    np.add.at(state_species,(inv,states_within),1)
    state_species/=counts[:,None]
    states_among=states_from_three_values(means-np.median(means,axis=0))
    return labels,state_species,states_among


def analyze_triad(values,group,columns=TRIAD,replicates=400,seed=20261024):
    if replicates<100 or replicates>5000:
        raise ValueError("INVALID_BIOLOGICAL_CLUSTER_REPLICATES")
    labels,within_by_taxon,among_states=source_states(values,group,columns)
    n=len(labels)
    within=np.mean(within_by_taxon,axis=0)
    among=np.bincount(among_states,minlength=8).astype(float)/n
    base={"within":maxent_from_observed(within),
          "among":maxent_from_observed(among)}
    rng=np.random.default_rng(seed)
    draws={"within":[],"among":[]}
    for _ in range(replicates):
        idx=rng.integers(0,n,size=n)
        pwithin=within_by_taxon[idx].mean(axis=0)
        pamong=np.bincount(among_states[idx],minlength=8).astype(float)/n
        for name,distribution in (("within",pwithin),("among",pamong)):
            try:
                m=maxent_from_observed(distribution)
                draws[name].append(m["observed_minus_pairwise_maxent"])
            except ValueError:
                draws[name].append(np.nan)
    out={}
    for name in base:
        arr=np.asarray(draws[name],float)
        good=arr[np.isfinite(arr)]
        out[name]={
            **base[name],
            "n_source_taxa":n,
            "taxon_bootstrap_replicates":replicates,
            "bootstrap_finite_replicates":len(good),
            "observed_minus_maxent_bootstrap_range95":
                list(map(float,np.quantile(good,[.025,.975]))) if len(good)>=.9*replicates else None,
            "bootstrap_is_descriptive_not_independent_evolutionary_replicates":True,
        }
    return {
        "model":"three-way visible phenotype configuration beyond all single/pairwise marginals",
        "triad":list(TRIAD),
        "binary_rule":{"within":"taxon-centred image-trait coordinate >0",
                       "among":"taxon-mean image trait > equal-taxon median"},
        "unselected_all_8_configuration_states":list(CODES),
        "within_and_among":out,
        "not_functional_resource_genetic_LEGO_test":True,
        "environmental_conditioning_does_not_prove_adaptation":True,
    }
