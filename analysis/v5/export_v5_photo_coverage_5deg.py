"""Export only privacy-filtered V5 5-degree photo measurement cells.

Dataframe is accepted only after the original V5 source/WCVP/CHELSA-5
restoration and three-trait measurements have been joined.  The raw
photo observations and their coordinates remain private in CI.
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from analysis.v4.native_cached_climate import VARIABLES

TRAITS = ["orientation_angle_degrees_median","corolla_lab_chroma_median","shape_aspect_ratio_median"]
EXPECTED = {
 ("Cirsium arvense","native"):18505,
 ("Cirsium arvense","introduced"):22269,
 ("Cirsium vulgare","native"):13115,
 ("Cirsium vulgare","introduced"):25532,
}
def export(data:pd.DataFrame, target:Path)->dict:
    needed=TRAITS+list(VARIABLES)+["latitude","longitude","taxon_name","range_state"]
    if not set(needed).issubset(data):
        raise ValueError("V5_SOURCE_FIELDS_NOT_COMPLETE")
    d=data.loc[data.taxon_name.isin({s for s,_ in EXPECTED}) &
               data.range_state.isin({"native","introduced"}),needed].copy()
    for c in TRAITS+list(VARIABLES)+["latitude","longitude"]:
        d[c]=pd.to_numeric(d[c], errors="coerce")
    d=d.replace([np.inf,-np.inf], np.nan).dropna(subset=needed[:-2])
    d=d.loc[d[TRAITS[0]].between(0,180) &
            d.latitude.between(-90,90) & d.longitude.between(-180,180)]
    counts=d.groupby(["taxon_name","range_state"]).size()
    got={key:int(counts.get(key,0)) for key in EXPECTED}
    if got!=EXPECTED:
        raise ValueError("V5_THREE_TRAIT_COHORT_COUNT_MISMATCH: "+repr(got))
    d=d.copy()
    d["lat_bin"]=np.minimum(35,np.floor((d.latitude+90)/5)).astype(int)
    d["lon_bin"]=np.minimum(71,np.floor((d.longitude+180)/5)).astype(int)
    cells=d.groupby(["taxon_name","range_state","lat_bin","lon_bin"],as_index=False).size()
    cells=cells.rename(columns={"taxon_name":"species","size":"record_count"})
    public=cells.loc[cells.record_count>=5].sort_values(["species","range_state","lat_bin","lon_bin"])
    target=Path(target);target.parent.mkdir(parents=True,exist_ok=True)
    public.to_csv(target,index=False)
    receipt={
      "source":"V5 exact three-trait and five-annual-climate photo cohort",
      "total_records":sum(EXPECTED.values()),
      "stratum_totals":{s+"_"+r:int(n) for (s,r),n in EXPECTED.items()},
      "shown_cells":int(len(public)),
      "shown_records":int(public.record_count.sum()),
      "suppressed_records":int(cells.record_count.sum()-public.record_count.sum()),
      "privacy":"only 5-degree bins with at least five records; no individual ids or coordinates",
      "csv_sha256":hashlib.sha256(target.read_bytes()).hexdigest(),
    }
    target.with_suffix(".receipt.json").write_text(
      json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    return receipt
