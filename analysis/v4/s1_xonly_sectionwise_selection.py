"""Strictly outcome-blind V4 S1 covariate selection for sections 3.1 and 3.3.

Reads authenticated environmental predictors, observation taxon and coordinates.
No trait values, trait eligibility, previous analysis P values, RF importances,
or PLS weights enter selection. Returns anonymous summaries only.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v4.species_environment_existing_inputs import load_existing_inputs
from analysis.v4.species_environment_scan import PREDICTORS
from analysis.v4.xonly_climate_selection import select_climate_variables

CLIMATE=tuple(s for s in PREDICTORS if s!="topo_elevation")
assert len(CLIMATE)==22

def cells_by_taxon(frame,degree,*,minimum_photos=5):
    data=frame.loc[frame["accepted_key"].notna(),
                   ["accepted_key","latitude","longitude",*CLIMATE]].copy()
    data=data.dropna(subset=CLIMATE)
    valid=(data.latitude.between(-90,90)&data.longitude.between(-180,180))
    data=data.loc[valid].copy()
    data["a"]=np.floor((data.latitude+90)/degree).astype(int)
    data["b"]=np.floor((data.longitude+180)/degree).astype(int)
    g=data.groupby(["accepted_key","a","b"],observed=True,sort=True)
    unit=g[list(CLIMATE)].median()
    unit["source_n"]=g.size()
    unit=unit.reset_index()
    summed=unit.groupby("accepted_key",sort=True,observed=True)["source_n"].sum()
    eligible=summed.index[summed>=minimum_photos]
    return unit.loc[unit.accepted_key.isin(eligible)].copy(),{
        "source_complete_climate_photos":int(len(data)),
        "source_taxa_with_minimum_photo_support":int(len(eligible)),
        "taxon_cell_units":int(len(unit))}

def run_source(frame,receipt_sha):
    required={"obs_id","accepted_key","latitude","longitude",*CLIMATE}
    if not required <= set(frame) or len(frame)!=313350:
        raise ValueError("S1_X_ONLY_SOURCE_GRAIN_INVALID")
    # Section 3.1: among taxa and within taxa distinguish different X structures.
    unit,coverage=cells_by_taxon(frame,1)
    group=unit.groupby("accepted_key",sort=True,observed=True)
    between=group[list(CLIMATE)].mean()
    if len(between)<50:raise ValueError("S1_X_ONLY_FEW_TAXA")
    between_select=select_climate_variables(
        between.to_numpy(float),CLIMATE,bootstraps=24,seed=20261107)
    counts=group.size()
    eligible=counts[counts>=3].index
    within=unit.loc[unit.accepted_key.isin(eligible)].copy()
    centered=within[list(CLIMATE)]-within.groupby("accepted_key")[list(CLIMATE)].transform("mean")
    within_weights=1/within.accepted_key.map(within.groupby("accepted_key").size()).to_numpy(float)
    within_select=select_climate_variables(
        centered.to_numpy(float),CLIMATE,weights=within_weights,bootstraps=24,seed=20261108)
    # Section 3.3: select X only in geography-repeated source climate cells.
    # This is intentionally NOT conditioned on measured flower-trait values.
    h3={}
    for degree in (1,5):
        sample,meta=cells_by_taxon(frame,degree,minimum_photos=40)
        g=sample.groupby("accepted_key",observed=True,sort=True)
        retained=g.size();retained=retained[retained>=8].index
        sample=sample.loc[sample.accepted_key.isin(retained)].copy()
        if sample.accepted_key.nunique()<25:raise ValueError("S1_H3_X_ONLY_CELLS_UNSUPPORTED")
        weights=1/sample.accepted_key.map(sample.groupby("accepted_key").size()).to_numpy(float)
        selection=select_climate_variables(
            sample[list(CLIMATE)].to_numpy(float),CLIMATE,weights=weights,
            bootstraps=24,seed=20261107+degree)
        h3[str(degree)+"degree"]={"n_source_taxa":int(sample.accepted_key.nunique()),
                                 "n_source_cells":int(len(sample)),
                                 "selected_climate_X":selection}
    output={"schema":"azami_v4_sectionwise_xonly_s1_v1",
            "status":"original_S1_authenticated_X_only_blind_selection",
            "original_native_observations":313350,
            "climate_source_receipt_sha256":str(receipt_sha),
            "candidate_predictors":list(CLIMATE),
            "uses_original_or_previous_flower_trait_results":False,
            "uses_photographic_trait_columns":False,
            "elevation_separate_geographic_covariate":True,
            "section31":{
                "between_taxa_n":int(len(between)),
                "within_taxa_n":int(len(eligible)),
                "within_taxon_cell_units":int(len(within)),
                "climate_complete_source_support":coverage,
                "between_taxa_X_only":between_select,
                "within_taxa_X_only":within_select},
            "section33":{"role":"source_geographic_climate_distance_representation",
                         "comparison_set_has_not_been_chosen_by_H3c_p_values":True,
                         "by_resolution":h3},
            "fixed_GEB_S1_S2_H3c_untouched":True}
    banned={"obs_id","accepted_key","taxon","taxon_name","latitude","longitude",
            "user_id","photographer_id","photo_id","individual_prediction"}
    def audit(obj):
        if isinstance(obj,dict):
            for key,val in obj.items():
                if key in banned:raise ValueError("S1_XONLY_PRIVATE_IDENTIFIER_EXPORTED")
                audit(val)
        elif isinstance(obj,list):
            for val in obj:audit(val)
    audit(output)
    return output

def run_existing(native,decoded_receipt,seven,nine,nine_manifest,source_receipt,
                 source_sha,output):
    table,receipt=load_existing_inputs(native,decoded_receipt,seven,nine,
                                       source_receipt,nine_manifest,
                                       source_receipt_sha256=source_sha,expected_n=313350)
    geo=pd.read_csv(native,usecols=["obs_id","latitude","longitude"],
                    dtype={"obs_id":str},low_memory=False)
    if len(geo)!=313350 or not geo.obs_id.is_unique:raise ValueError("S1_XONLY_COORD_ID_INVALID")
    ids=table.obs_id.copy()
    table=table.merge(geo,on="obs_id",validate="one_to_one",sort=False,how="left")
    if not table.obs_id.equals(ids):raise ValueError("S1_XONLY_ORIGINAL_JOIN_DRIFT")
    # Drop all response columns BEFORE calling any selection algorithm.
    source=table[["obs_id","accepted_key","latitude","longitude",*CLIMATE]].copy()
    result=run_source(source,source_sha)
    result["original_native_table_sha256"]=receipt["input_digests"]["native"]
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    return result

def main():
    p=argparse.ArgumentParser()
    for f in ("native","decoded-receipt","seven","nine","nine-manifest",
              "source-receipt","source-sha","output"):
        p.add_argument("--"+f,required=True)
    a=p.parse_args()
    result=run_existing(a.native,a.decoded_receipt,a.seven,a.nine,
                        a.nine_manifest,a.source_receipt,a.source_sha,a.output)
    print(json.dumps({"status":result["status"],
                      "between":result["section31"]["between_taxa_X_only"]["selected_variables"],
                      "within":result["section31"]["within_taxa_X_only"]["selected_variables"],
                      "h3_1deg":result["section33"]["by_resolution"]["1degree"]["selected_climate_X"]["selected_variables"]}))
if __name__=="__main__":main()
