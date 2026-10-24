"""V4 phenotypic three-way head configurations under a pairwise-only null.

A direct *visible phenotype* test inspired by aza3 LEGO's need to distinguish
joint module configuration from marginal/pairwise information. It is NOT a
genomic co-retention (G123) assay, and same photo phenotypes are NOT plant
genetic parts. The three endpoints are fixed: image angle, colour chroma,
head-outline aspect. Other features continue in the full 28-pair analysis.

Each taxon contributes equal weight, and data are dichotomized by its own
observed median for this bounded descriptive test. Pairwise-only IPF is a
statistical log-linear reference, not a biological evolutionary null.
"""
from __future__ import annotations
import numpy as np
from analysis.v4.s1_s2_extreme_lego import FEATURE_COLS

TRIPLE = ("orientation_image_vertical_angle",
          "corolla_lab_chroma",
          "capitulum_outline_aspect_ratio")
IDX = [FEATURE_COLS.index(t) for t in TRIPLE]
MIN_OBS_PER_TAXON = 12


def taxon_tables(values,groups):
    v = np.asarray(values,float)[:,IDX]
    labels,inv=np.unique(np.asarray(groups),return_inverse=True)
    tables=[]; skipped=0; used_n=0
    for i in range(len(labels)):
        part=v[inv==i]
        if len(part)<MIN_OBS_PER_TAXON:
            skipped+=1
            continue
        median=np.median(part,axis=0)
        bits=(part>median).astype(int)
        if not all(np.unique(bits[:,j]).size==2 for j in range(3)):
            skipped+=1
            continue
        code=4*bits[:,0]+2*bits[:,1]+bits[:,2]
        counts=np.bincount(code,minlength=8).astype(float).reshape(2,2,2)
        tables.append(counts/counts.sum())
        used_n+=len(part)
    if len(tables)<15:
        return None,{"status":"INSUFFICIENT_INDEPENDENT_TAXON_CONFIGURATION_SUPPORT",
                     "taxa_eligible":len(tables),"taxa_excluded":skipped}
    return np.asarray(tables),{"taxa_eligible":len(tables),"taxa_excluded":skipped,
                               "source_observations_eligible":used_n}


def pairwise_max_entropy(table,iterations=1000,tol=1e-11):
    observed=np.asarray(table,float)
    if observed.shape!=(2,2,2) or np.min(observed)<0 or observed.sum()<=0:
        raise ValueError("INVALID_THREE_WAY_SOURCE_MATRIX")
    observed=observed/observed.sum()
    expected=np.ones((2,2,2),float)/8
    for step in range(iterations):
        prior=expected.copy()
        for axis in (2,1,0):
            target=observed.sum(axis=axis,keepdims=True)
            estimate=expected.sum(axis=axis,keepdims=True)
            expected *= np.divide(target,estimate,out=np.zeros_like(target),
                                  where=estimate>0)
        if np.max(np.abs(expected-prior))<tol:
            break
    if abs(expected.sum()-1)>1e-8:
        raise ValueError("PAIRWISE_IPF_MARGINAL_MASS_CHANGED")
    for axis in (0,1,2):
        if np.max(np.abs(observed.sum(axis=axis)-
                         expected.sum(axis=axis)))>1e-7:
            raise ValueError("PAIRWISE_IPF_DID_NOT_PRESERVE_TWO_WAY_MARGINS")
    return expected


def three_way_excess(table):
    p=np.asarray(table,float)
    q=pairwise_max_entropy(p)
    # Even vs odd parity, a contrast of the 3-way tensor beyond the
    # pairwise-only max-entropy table, not necessarily an independent trait.
    sign=np.array([[[1.,-1.],[-1.,1.]],
                   [[-1.,1.],[1.,-1.]]])
    parity=float(np.sum(sign*(p/p.sum()-q)))
    eps=1e-12
    dev=float(2*np.sum(p/p.sum()*np.log((p/p.sum()+eps)/(q+eps))))
    return parity,dev,q


def test_three_way(values,groups,reps=499,seed=20261024):
    """Taxon-cluster resampling, not a p-value for inherited LEGO parts."""
    tables,coverage=taxon_tables(values,groups)
    if tables is None:
        return coverage
    observed=np.mean(tables,axis=0)
    parity,divergence,expected=three_way_excess(observed)
    n_taxa=len(tables)
    rng=np.random.default_rng(seed)
    draws=[]
    for _ in range(reps):
        subset=tables[rng.integers(n_taxa,size=n_taxa)].mean(axis=0)
        try:
            moment,_,_=three_way_excess(subset)
            draws.append(moment)
        except ValueError:
            continue
    if len(draws)<reps*.95:
        return {**coverage,"status":"INSUFFICIENT_FINITE_TAXON_BOOTSTRAPS",
                "bootstrap_valid":len(draws)}
    occupied=(tables>0).sum(axis=(1,2,3))
    return {**coverage,"status":"OBSERVED_VISIBLE_THREE_WAY_CONFIGURATION",
            "component_ids":TRIPLE,
            "binary_rule":"above each source taxon's own observed median",
            "taxa_with_all_8_visible_configurations":int(np.sum(occupied==8)),
            "median_configurations_per_taxon":float(np.median(occupied)),
            "observed_taxon_equal_state_fraction":np.round(observed.reshape(-1),8).tolist(),
            "max_entropy_pairwise_only_fraction":np.round(expected.reshape(-1),8).tolist(),
            "three_way_excess_parity":parity,
            "pairwise_only_reference_KL_twice":divergence,
            "cluster_resampling_95pct_parity_interval":np.quantile(draws,[.025,.975]).tolist(),
            "cluster_resampling_fraction_positive":float(np.mean(np.asarray(draws)>0)),
            "bootstrap_reps":len(draws),
            "not_a_test_of_genetic_resource_co_retention":True,
            "photo_camera_colour_and_head_geometry_may_be_correlated":True}
