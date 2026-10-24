#!/usr/bin/env python3
"""Exact label null for eight V4 image-defined capitulum constructs, not genetic LEGO.

Input MUST be frozen V4 S1/S2 original-pixel result ZIP's aggregate JSON.
Enumerate all 8!/(1!3!4!) = 280 allocations of one presentation, three colour,
four outline traits; hold the 28 observed RV edges fixed. One-sided label P
tests whether prespecified module labels concentrate strong RV edges more than
equally sized arbitrary labels. It does not correct shared image algorithms.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

INPUT_SHA="f1fc39fa829318ce7289b27857255992f00b941094faf16d3b43812a77c04209"
MODULE_COUNTS={"presentation":1,"colour":3,"outline":4}


def exact(pairs):
    if len(pairs)!=28:
        raise ValueError("EXPECTED_28_V4_PAIRS")
    names=sorted({r[k] for r in pairs for k in ("left","right")})
    if len(names)!=8 or len({(r["left"],r["right"]) for r in pairs})!=28:
        raise ValueError("INCONSISTENT_V4_TRAIT_PAIRS")
    real={}
    for row in pairs:
        for end in ("left","right"):
            if end in row and row[end] in real and real[row[end]]!=row[end+"_module"]:
                raise ValueError("CONFLICTING_MODULE_ASSIGNMENT")
            real[row[end]]=row[end+"_module"]
    if {k:list(real.values()).count(k) for k in MODULE_COUNTS}!=MODULE_COUNTS:
        raise ValueError("INVALID_FROZEN_MODULE_SIZE")
    labelings=[]
    for presentation in itertools.combinations(names,1):
        remaining=[x for x in names if x not in presentation]
        for colour in itertools.combinations(remaining,3):
            labelings.append({name:("presentation" if name in presentation else
                                   "colour" if name in colour else "outline")
                              for name in names})
    if len(labelings)!=280:
        raise ValueError("EXPECTED_280_EXACT_ALLOCATIONS")
    out={}
    for scale in ("within_rv","among_rv"):
        def median_gap(labels):
            within=[r[scale] for r in pairs if labels[r["left"]]==labels[r["right"]]]
            cross=[r[scale] for r in pairs if labels[r["left"]]!=labels[r["right"]]]
            if len(within)!=9 or len(cross)!=19:
                raise ValueError("MODULE_PAIR_NUMBER_CHANGED")
            return float(np.median(within)-np.median(cross))
        observed=median_gap(real)
        null=np.array([median_gap(z) for z in labelings])
        out[scale]={"observed_module_minus_cross_median_RV":observed,
                    "one_sided_exact_label_p":float((null>=observed-1e-14).mean()),
                    "null_95pct_range":[float(v) for v in np.quantile(null,[.025,.975])],
                    "n_unique_allocations":280,
                    "conditional_on_observed_RV":True}
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--aggregate",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    v=a.parse_args()
    raw=v.aggregate.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=INPUT_SHA:
        raise ValueError("UNVERIFIED_OR_CHANGED_V4_SOURCE_AGGREGATE")
    data=json.loads(raw)
    s=data["S1"]
    if (data["analysis_id"]!="azami_v4_S1_S2_extreme_climate_and_phenotypic_LEGO_20261024"
        or s["common_all_eight_trait_n"]!=156683
        or s["common_all_eight_trait_taxa"]!=174):
        raise ValueError("V4_SOURCE_COHORT_DRIFT")
    results={"analysis_id":"v4_s1_s2_exact_phenotype_module_label_null_20261024",
             "frozen_aggregate_sha":INPUT_SHA,
             "base_original":exact(s["raw_28_pairs"]),
             "climate_scenarios":{},
             "interpretation":"Conditional post-hoc label null for visual measurement modules, not independently inherited genetic LEGO or adaptation."}
    for name,scenario in s["scenario_results"].items():
        if scenario["status"]=="ESTIMATED_POSTHOC_EXPLORATORY":
            results["climate_scenarios"][name]=exact(scenario["pairs"])
        else:
            results["climate_scenarios"][name]={"status":scenario["status"]}
    v.output.parent.mkdir(parents=True,exist_ok=True)
    v.output.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"raw":results["base_original"],
                      "n_scenarios":len(results["climate_scenarios"])},indent=2))


if __name__=="__main__":
    main()
