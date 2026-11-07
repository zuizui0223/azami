"""Outcome-blind S2 source climate selection for native/introduced shared adjustment.

Requires only original 297622 observation IDs, climate, region, observer and
geographic/date fields. NO floral outcome, prior significance, RF rank or prior
regression effects can enter the selector. Elevation is diagnosed separately.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v4.s1_s2_extreme_lego import ALL_ENV
from analysis.v4.s1_s2_extreme_range_models import _get_dated_balanced
from analysis.v4.xonly_climate_selection import select_climate_variables

CLIMATE=tuple(x for x in ALL_ENV if x!="topo_elevation")
ELEVATION="topo_elevation"
SPECIES=("Cirsium arvense","Cirsium vulgare")
assert len(CLIMATE)==11 and ELEVATION in ALL_ENV

def common_90pct_interval(a,b):
    a5,a95=np.quantile(a,[.05,.95]);b5,b95=np.quantile(b,[.05,.95])
    lo=max(a5,b5);hi=min(a95,b95)
    if lo>hi:return 0.,0.
    return float(np.mean((a>=lo)&(a<=hi))),float(np.mean((b>=lo)&(b<=hi)))

def select_from_origins(focal,additional):
    if (len(focal)!=297622 or len(additional)!=297622 or
        not focal.obs_id.is_unique or not additional.obs_id.is_unique or
        set(focal.obs_id)!=set(additional.obs_id)):
        raise ValueError("S2_ORIGINAL_ID_OR_COUNT_MISMATCH")
    # Never read/merge a photographic measurement table.
    needed=["obs_id",*[k for k in ALL_ENV if k not in focal]]
    if not set(needed)<=set(additional) or (set(needed)-{"obs_id"}) & set(focal):
        raise ValueError("S2_ENVIRONMENT_COLUMN_DRIFT")
    joined=focal.merge(additional[needed],on="obs_id",validate="one_to_one",how="left",sort=False)
    if not set(ALL_ENV)<=set(joined):raise ValueError("S2_INCOMPLETE_12_ENV_COLUMNS")
    if set(joined.taxon_name)!=set(SPECIES) or set(joined.range_state)!={"native","introduced"}:
        raise ValueError("S2_SOURCE_SPECIES_OR_RANGE_CHANGED")
    summaries=[]
    for species in SPECIES:
        cells={};support={}
        for status in ("native","introduced"):
            rows=joined.loc[(joined.taxon_name==species)&(joined.range_state==status)]
            n_original=len(rows)
            rows=rows.loc[np.isfinite(rows[list(ALL_ENV)].to_numpy(float)).all(axis=1)].copy()
            balanced=_get_dated_balanced(rows)
            centres=balanced.groupby("_cell",sort=True)[list(ALL_ENV)].median()
            if len(centres)<15:raise ValueError("S2_BLIND_GEO_CELLS_INSUFFICIENT")
            cells[status]=centres
            support[status]={"original_observation_n":int(n_original),
                             "all12_finite_n":int(len(rows)),
                             "observer_cell_unique_n":int(len(balanced)),
                             "geographic_cells_n":int(len(centres))}
        a,b=cells["native"],cells["introduced"]
        overlap={}
        available=[]
        for p in CLIMATE:
            pa,pb=common_90pct_interval(a[p].to_numpy(float),b[p].to_numpy(float))
            keep=bool(min(pa,pb)>=.10)
            overlap[p]={"native_cell_common_interval_fraction":pa,
                        "introduced_cell_common_interval_fraction":pb,
                        "eligible_on_X_support_only":keep}
            if keep:available.append(p)
        if len(available)<3:raise ValueError("S2_BLIND_NO_COMMON_CLIMATE_AXES")
        pooled=np.vstack([a[available].to_numpy(float),b[available].to_numpy(float)])
        w=np.r_[np.full(len(a),.5/len(a)),np.full(len(b),.5/len(b))]
        chosen=select_climate_variables(pooled,available,weights=w,
                                        target=.90,max_k=min(8,len(available)),
                                        bootstraps=24,
                                        seed=20261107+SPECIES.index(species))
        # Elevation is evaluated as an X-only residual geographic coordinate.
        alt=np.r_[a[ELEVATION].to_numpy(float),b[ELEVATION].to_numpy(float)]
        subset=[available.index(p) for p in chosen["selected_variables"]]
        xx=np.column_stack([np.ones(len(pooled)),pooled[:,subset]])
        beta=np.linalg.lstsq(xx*np.sqrt(w)[:,None],alt*np.sqrt(w),rcond=1e-10)[0]
        mu=np.average(alt,weights=w)
        ss=np.dot(w,(alt-mu)**2)
        residual_var=float(np.dot(w,(alt-xx@beta)**2)/ss) if ss>1e-12 else None
        summaries.append({"species":species,"native_and_introduced_source_support":support,
                          "candidate_climate_overlap_X_only":overlap,
                          "selected_shared_covariate_set_X_only":chosen,
                          "elevation_residual_X_variance_fraction":residual_var,
                          "uses_same_selected_climate_columns_in_both_ranges":True})
    return {"schema":"azami_s2_outcome_blind_common_climate_selection_v1",
            "status":"ORIGINAL_S2_ENVIRONMENT_ONLY",
            "original_source_observation_n":297622,
            "candidate_climate_names":list(CLIMATE),
            "geography_topographic_elevation_separate":True,
            "trait_values_not_loaded":True,
            "past_results_P_values_and_RF_rankings_not_loaded":True,
            "source_note":"original WCVP two species; no extra proxy species",
            "results":summaries}

def run(focal_path,extra_path,output):
    original=pd.read_csv(focal_path,dtype={"obs_id":str},low_memory=False)
    extra=pd.read_csv(extra_path,dtype={"obs_id":str},low_memory=False)
    expected={"obs_id","taxon_name","range_state","user_id","latitude",
              "longitude","observed_on"}
    if not expected<=set(original):raise ValueError("S2_ORIGINAL_SOURCE_COORDINATES_MISSING")
    result=select_from_origins(original,extra)
    banned={"obs_id","accepted_key","taxon_name","latitude","longitude",
            "user_id","photo_id","photographer_id","individual_prediction"}
    def audit(node):
        if isinstance(node,dict):
            for key,value in node.items():
                if key in banned:raise ValueError("S2_BLIND_EXPORT_PRIVATE_KEY")
                audit(value)
        elif isinstance(node,list):
            for value in node:audit(value)
    audit(result)
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    return result

def main():
    p=argparse.ArgumentParser()
    for key in ("focal","extra","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    result=run(a.focal,a.extra,a.output)
    print(json.dumps({"status":result["status"],"chosen":{
        x["species"]:x["selected_shared_covariate_set_X_only"]["selected_variables"]
        for x in result["results"]}}),flush=True)
if __name__=="__main__":main()
