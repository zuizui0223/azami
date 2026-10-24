"""Independent external construct-check for actual S1 photographic orientation.

Compare S1 *continuous* image-angle medians against 46 source-authored U/D
botanical states. No thresholding/relabeling of remaining S1 194 taxa; test
only whether externally scored nodding taxa rank higher on a photographic
axis. Source strata NA42/EA4 fixed; this is NOT gravity validation or new
independent field data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from analysis.v4.t2_source_blocked_five_climate import (
    authenticate_input, RECEIPT, NA, EA, cohort, exact_source_maps
)


def rank_auc(scores, y):
    """Mann-Whitney ROC AUC, including half credit for ties."""
    s=np.asarray(scores,float)
    labels=np.asarray(y,int)
    if s.ndim!=1 or s.shape!=labels.shape or not np.all(np.isfinite(s)):
        raise ValueError("Invalid measured continuous angles")
    if not set(np.unique(labels)) <= {0,1}:
        raise ValueError("External morphology state must be binary")
    n1=int(labels.sum());n0=len(labels)-n1
    if n1==0 or n0==0:
        raise ValueError("No botanical-state contrast in source")
    return float((rankdata(s,method="average")@labels -
                 n1*(n1+1)/2)/(n1*n0))


def source_blocked_exact_auc(rows):
    """Exact source-prevalence permutation test for direction consistency."""
    rows=sorted(rows,key=lambda z:z["source_taxon"])
    c=cohort(rows)
    values=np.array([float(r["s1_image_orientation_median"]) for r in rows])
    if np.any((values<0)|(values>180)):
        raise ValueError("S1 angle outside stated image-axis range")
    labels,nmaps=exact_source_maps(c)
    strata={}
    per_source_arrays=[]
    for key in sorted(set(c["lane"])):
        pos=np.flatnonzero(c["lane"]==key)
        z=c["y"][pos]
        if z.sum()<=0 or z.sum()==len(z):
            raise ValueError("Source lacks one of the independent head states")
        ranking=rankdata(values[pos],method="average")
        k=int(z.sum());m=len(z)-k
        null=((labels[:,pos]@ranking)-k*(k+1)/2)/(k*m)
        observed=rank_auc(values[pos],z)
        category="North_American_published_author" if key==0 else "East_Asian_EAzami_curated"
        strata[category]={
           "n_taxa":len(z),"n_nodding":k,"n_erect":m,
           "observed_AUC":observed,
           "exact_null_2sided_around_half_p":float(np.mean(np.abs(null-.5)>=abs(observed-.5)-1e-12)),
           "exact_null_upper_p":float(np.mean(null>=observed-1e-12)),
           "source_specific_total_possible_assignments":int(np.unique(labels[:,pos],axis=0).shape[0]),
        }
        per_source_arrays.append((null,observed,k*m))
    pooled=rank_auc(values,c["y"])
    numerator=sum(x[0]*x[2] for x in per_source_arrays)
    denominator=sum(x[2] for x in per_source_arrays)
    null=numerator/denominator
    observed=sum(x[1]*x[2] for x in per_source_arrays)/denominator
    source_equal_null=np.mean([x[0] for x in per_source_arrays],axis=0)
    source_equal_obs=float(np.mean([x[1] for x in per_source_arrays]))
    return {
        "status":"EMPIRICAL_EXTERNAL_SOURCE_MORPHOLOGY_IMAGE_ANGLE_RANK_CHECK",
        "n_S1_angle_taxa":len(rows),
        "n_botanical_nodding":int(c["y"].sum()),
        "n_S1_eligible_image_taxon_keys_total":244,
        "n_S1_tree_taxa_with_continuous_angle":194,
        "image_angle_measurement":"EXIF_image_vertical_not_gravity",
        "external_binary_states_not_derived_from_image_angle":True,
        "source_prevalence_exact_assignments":nmaps,
        "pooled_mixture_AUC_descriptive_not_source_adjusted":pooled,
        "within_source_pair_weighted_AUC":float(observed),
        "within_source_pair_weighted_source_blocked_two_sided_p":float(np.mean(np.abs(null-.5)>=abs(observed-.5)-1e-12)),
        "within_source_pair_weighted_upper_p":float(np.mean(null>=observed-1e-12)),
        "within_source_equal_weight_AUC":source_equal_obs,
        "within_source_equal_weight_two_sided_p":float(np.mean(np.abs(source_equal_null-.5)>=abs(source_equal_obs-.5)-1e-12)),
        "source_strata":strata,
        "direction_convention":"AUC>0.5 means independent botanical nodding tips have larger IMAGE-vertical angles within the original source lane",
        "cannot_convert_194_image_angle_to_true_nodding_independent_gain_counts":True,
        "sample_selection_and_external_taxonomic_source_not_independent_validation_of_rainfall_causation":True,
    }


def run(source,output):
    raw=json.loads(Path(source).read_text(encoding="utf-8"))
    receipt=json.loads(RECEIPT.read_text(encoding="utf-8"))
    rows=authenticate_input(raw,receipt)
    report=source_blocked_exact_auc(rows)
    report["original_source_run_id"]=receipt["workflow_run_id"]
    report["original_source_artifact_digest"]=receipt["workflow_artifact_digest"]
    report["independent_morphological_coding_sources"]=[NA,EA]
    dest=Path(output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                    encoding="utf-8")
    print("S1_ANGLE_EXTERNAL_U_D_CALIBRATION_SUMMARY="+json.dumps(report,
                         ensure_ascii=False,sort_keys=True),flush=True)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--joint",required=True)
    ap.add_argument("--out",required=True)
    o=ap.parse_args()
    run(o.joint,o.out)


if __name__=="__main__":
    main()
