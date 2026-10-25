"""Exact descriptive module-label allocation reference on observed V4 28 RV edges.

There are exactly 8 image constructs divided into 1 presentation, 3 colour and
4 outline coordinates. Enumerating 8*C(7,3) = 280 possible assignments
assesses whether the registered grouping shows unusually strong *within-block*
phenotypic association compared with arbitrary equal-sized group labels.
It is NOT a phylogenetic null, gene/developmental LEGO proof, independent
statistical replication, or an independent permutation of photo measurements.
Colour and outline coordinates share source pixels/mathematical definitions;
a strong label effect can be entirely technical.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import numpy as np

CONSTRUCTS = ("angle","lightness","chroma","hue","aspect","circularity",
              "solidity","width_profile")
GROUPS = {"angle":"presentation",
          "lightness":"colour","chroma":"colour","hue":"colour",
          "aspect":"outline","circularity":"outline","solidity":"outline",
          "width_profile":"outline"}
EXPECTED_PAIRS=28
EXPECTED_NULLS=280


def all_partitions():
    for special in CONSTRUCTS:
        rest=[v for v in CONSTRUCTS if v!=special]
        for colour in itertools.combinations(rest,3):
            yield {name:(0 if name==special else 1 if name in colour else 2)
                   for name in CONSTRUCTS}


def contrast(rows,mapping,scale):
    same,cross=[],[]
    for item in rows:
        left,right=item["left"],item["right"]
        if left not in mapping or right not in mapping:
            raise ValueError("UNREGISTERED_V4_CONSTRUCT")
        (same if mapping[left]==mapping[right] else cross).append(float(item[scale]))
    if len(same)!=9 or len(cross)!=19:
        raise ValueError("V4_GROUP_SIZE_INCONSISTENT")
    return float(np.median(same)-np.median(cross)),float(np.median(same)),float(np.median(cross))


def evaluate(rows):
    canonical={tuple(sorted((r["left"],r["right"]))) for r in rows}
    actual={tuple(sorted(p)) for p in itertools.combinations(CONSTRUCTS,2)}
    if len(rows)!=EXPECTED_PAIRS or canonical!=actual:
        raise ValueError("V4_28_PAIR_SOURCE_IDENTITY_MISMATCH")
    original={name:(0 if group=="presentation" else 1 if group=="colour" else 2)
              for name,group in GROUPS.items()}
    partitions=list(all_partitions())
    if len(partitions)!=EXPECTED_NULLS:
        raise ValueError("WRONG_FIXED_SIZE_LABEL_NULL_SPACE")
    result={}
    for scale in ("within_rv","among_rv"):
        observed,same,cross=contrast(rows,original,scale)
        null=np.array([contrast(rows,p,scale)[0] for p in partitions],float)
        p=float(np.mean(null>=observed-1e-12))
        result[scale]={"observed_within_module_minus_cross_module_median_RV":observed,
                       "observed_within_module_median_RV":same,
                       "observed_cross_module_median_RV":cross,
                       "exact_one_sided_module_label_reference_p":p,
                       "null_280_label_allocations":len(partitions),
                       "null_reference_quantiles_q025_q50_q975":
                           list(map(float,np.quantile(null,[.025,.5,.975]))),
                       "explicitly_not_an_inferential_evolutionary_p_value":True}
    return result


def report(aggregate):
    if aggregate.get("S1",{}).get("common_all_eight_trait_taxa",0)<20:
        raise ValueError("ORIGINAL_V4_COHORT_NOT_ESTABLISHED")
    s=aggregate["S1"]
    results={"raw":evaluate(s["raw_28_pairs"])}
    for name,data in s["scenario_results"].items():
        if data["status"]=="ESTIMATED_POSTHOC_EXPLORATORY":
            results[name]=evaluate(data["pairs"])
    return {
        "analysis_id":"azami_v4_s1_fixed_size_module_label_reference_20261024",
        "original_cohort_observations":s["common_all_eight_trait_n"],
        "original_cohort_taxa":s["common_all_eight_trait_taxa"],
        "registered_block_sizes":{"presentation":1,"colour":3,"outline":4},
        "pairwise_tests_derived_only_from_existing_V4_aggregate":28,
        "modules_mathematically_shared_image_coordinates_not_independent_gene_lego":True,
        "not_evolutionary_phylogenetic_randomization":True,
        "not_higher_order_genetic_resource_test":True,
        "descriptive_reference":results,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--aggregate",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    data=report(json.loads(a.aggregate.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"taxa":data["original_cohort_taxa"],
                      "raw":data["descriptive_reference"]["raw"]}))
if __name__=="__main__":
    main()
