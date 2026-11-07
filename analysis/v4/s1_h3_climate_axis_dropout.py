"""V4 3.3: source-authenticated climate analogue leave-one-axis-out sensitivity.

Existing positive H3c hypothesis and 24x199 permutation source run stay frozen.
This new run uses original 3 *climate* variables, not phenotype-selected
climate axes. Variable sets pre-fixed to full 3 and three 2-axis omissions.
Dimension adjustment sqrt(d/3) is outcome-blind; comparison is exploratory.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from analysis.v4 import s1_h3_climate_analogue_alignment as h3
from analysis.v4.s1_climate_trait_hypervolume import (
    load_source,geo_aggregate,standardise_cells
)

OPTIONS={
    "BIO12_BIO01_RSDS":(0,1,2),
    "BIO12_RSDS":(0,2),
    "BIO01_RSDS":(1,2),
    "BIO12_BIO01":(0,1),
}
DRAWS=12
PERMUTATIONS=99

def compare_cells(z,degree):
    results={}
    base_threshold=h3.CLIMATE_ANALOG_Z
    original_min=h3.MIN_SUCCESSFUL_DRAWS
    try:
        h3.MIN_SUCCESSFUL_DRAWS=6
        for label,indices in OPTIONS.items():
            altered=z.copy()
            original=z[["_z0","_z1","_z2"]].to_numpy(float)
            components=np.zeros_like(original)
            components[:,:len(indices)]=original[:,list(indices)]
            for i in range(3):
                altered["_z"+str(i)]=components[:,i]
            h3.CLIMATE_ANALOG_Z=float(base_threshold*np.sqrt(len(indices)/3))
            r=h3.analogue_analysis(altered,draws=DRAWS,
                                   permutations=PERMUTATIONS,
                                   seed=20261107+degree)
            results[label]={"climate_axis_count":len(indices),
                            "matched_climate_axis_indices":[int(i) for i in indices],
                            "threshold_scaled_for_dimensionality":
                                float(h3.CLIMATE_ANALOG_Z),
                            "analysis":r}
    finally:
        h3.CLIMATE_ANALOG_Z=base_threshold
        h3.MIN_SUCCESSFUL_DRAWS=original_min
    return results

def run(source,receipt,output):
    frame,prov=load_source(source,receipt)
    out={}
    for deg in (1,5):
        cells,support=geo_aggregate(frame,deg)
        z,scaling=standardise_cells(cells)
        out[str(deg)+"degree"]={
            "source_support":support,
            "source_equal_taxon_standardisation":scaling["equal_taxon_total_weight"],
            "climate_candidates":compare_cells(z,deg)}
    doc={"schema":"azami_v4_h3c_climate_axis_dropout_sensitivity_v1",
         "status":"source_original_H3c_exploratory_sensitivity",
         "source_provenance":prov,
         "n_original_s1_observations":313350,
         "original_H3c_24x199_Holm2_positive_untouched":True,
         "candidate_selection_only_source_climate_not_photo_outcome":True,
         "new_exploratory_draws_per_case":DRAWS,
         "new_exploratory_permutations_per_case":PERMUTATIONS,
         "threshold_rule":"0.9 * sqrt(number_of_climate_axes/3) on z scores",
         "n_evaluated_candidates_per_resolution":len(OPTIONS),
         "results":out}
    banned={"obs_id","user_id","accepted_key","photo_id",
            "taxon_name","latitude","longitude","photographer_id"}
    def audit(v):
        if isinstance(v,dict):
            for k,item in v.items():
                if k in banned:raise ValueError("PRIVATE_H3C_IDENTIFIER")
                audit(item)
        elif isinstance(v,list):
            for item in v:audit(item)
    audit(doc)
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    Path(output).write_text(json.dumps(doc,indent=2,ensure_ascii=False,
                                      allow_nan=False)+"\n",encoding="utf-8")
    return doc

def main():
    ap=argparse.ArgumentParser()
    for flag in ("source","receipt","output"):
        ap.add_argument("--"+flag,required=True)
    a=ap.parse_args()
    out=run(a.source,a.receipt,a.output)
    print(json.dumps({"status":out["status"],"n_original":out["n_original_s1_observations"],
                      "case_status":{k:{z:v["analysis"]["status"] for z,v in x["climate_candidates"].items()}
                        for k,x in out["results"].items()}}),flush=True)

if __name__=="__main__":main()
