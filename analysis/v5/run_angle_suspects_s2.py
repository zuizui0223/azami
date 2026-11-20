"""Run a bounded V5 S2 image-angle audit without rerunning downstream niche models.

Invoked only inside the private, source-authenticated GitHub Actions runner.
Exports aggregated diagnostic counts (no IDs, coordinates, taxon-level
records or measurements) for exactly 79,421 focal V5 source records.
"""
from pathlib import Path
import os
import pandas as pd
from analysis.v4.native_introduced_opportunity import authority_source,angular_cells
from analysis.v5.native_pca_env import background_mesh
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v5.audit_angle_suspects import export_s2

PRIVATE=Path("work/submission/private")
OUT=Path("work/submission/public")

def main():
    auth=authority_source((PRIVATE/"authority350.zip").read_bytes())
    wild=pd.read_csv(PRIVATE/"focal_wild.csv.gz",
                     dtype={"obs_id":str,"user_id":str},low_memory=False)
    if len(wild)!=297622 or wild.obs_id.duplicated().any():
        raise ValueError("ORIGINAL_WCVP_SOURCE_ID_CHANGED")
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(
        ["taxon_name","range_state","_cell","user_id"]
    ).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(
        wild,land,
        source_meta=Path("analysis/v4/environment_source_metadata.json"),
        work=PRIVATE/"angle_qc_rasters"
    )
    trait_cols=[
        "orientation_angle_degrees_median",
        "corolla_lab_chroma_median",
        "shape_aspect_ratio_median"
    ]
    measured=pd.read_csv(PRIVATE/"original-measured.csv.gz",
        usecols=["obs_id",*trait_cols],dtype={"obs_id":str})
    if not measured.obs_id.is_unique:
        raise ValueError("ORIGINAL_MEASUREMENT_DUPLICATE")
    wild=wild.drop(columns=trait_cols,errors="ignore").merge(
        measured,on="obs_id",how="left",validate="one_to_one")
    export_s2(wild,OUT/"v5_angle_suspect_s2.json")
    # Independently profile the original usable-head support and within-photo
    # angular dispersion on exactly the same five-climate, three-trait cohort.
    import json
    import numpy as np
    from analysis.v4.native_cached_climate import VARIABLES
    from analysis.v5.audit_angle_repeatability import summarise
    numerical=[*trait_cols,*VARIABLES]
    for col in numerical:
        wild[col]=pd.to_numeric(wild[col],errors="coerce")
    selected=wild.loc[
        np.isfinite(wild[numerical].to_numpy(float)).all(axis=1)
        &wild[trait_cols[0]].between(0,180)
    ]
    if len(selected)!=79421 or not selected.obs_id.is_unique:
        raise ValueError("S2_DIFFERENT_COHORT_IN_REPEATABILITY_GATE")
    original=pd.read_csv(PRIVATE/"original-measured.csv.gz",
                         dtype={"obs_id":str},low_memory=False)
    audit=summarise(selected.obs_id.to_numpy(),original)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"v5_angle_repeatability_s2.json").write_text(
        json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("V5_ANGLE_REPEATABILITY_S2_AGGREGATED_START")
    print(json.dumps(audit,separators=(",",":"),ensure_ascii=False))
    print("V5_ANGLE_REPEATABILITY_S2_AGGREGATED_END")

if __name__=="__main__":
    main()
