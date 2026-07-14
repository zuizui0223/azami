#!/usr/bin/env python3
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import pandas as pd


def load_base():
    path = Path(__file__).with_name("76_run_exhaustive_within_species_environment_rf.py")
    spec = importlib.util.spec_from_file_location("rf76", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--observation", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-batch-size", type=int, default=1024)
    args = p.parse_args()
    mod = load_base()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.observation, low_memory=False)
    lat = pd.to_numeric(df["latitude"], errors="coerce")
    lon = pd.to_numeric(df["longitude"], errors="coerce")
    if not (lat.notna() & lon.notna()).all():
        raise SystemExit("Invalid coordinates")
    df["spatial_block_10deg"] = [mod.spatial_block_10deg(a,b) for a,b in zip(lat,lon)]
    df, meta = mod.extract_environment(df, args.sample_batch_size)
    env_path = out / "strict_spatial_thinned_with_climate_topography_soil.csv"
    df.to_csv(env_path, index=False, encoding="utf-8-sig")
    climate=list(mod.CLIMATE); topo=list(mod.TOPOGRAPHY); soil=[f"soil_{p}_0_30cm" for p in mod.SOIL_PROPERTIES]
    cov=[]
    for name in climate+topo+soil:
        s=pd.to_numeric(df[name], errors="coerce")
        cov.append({"predictor":name,"n_nonmissing":int(s.notna().sum()),"complete_fraction":float(s.notna().mean())})
    pd.DataFrame(cov).to_csv(out/"environment_predictor_coverage.csv",index=False,encoding="utf-8-sig")
    (out/"environment_metadata.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")
    print({"observations":len(df),"species":df["taxon_name"].nunique(),"min_coverage":min(x["complete_fraction"] for x in cov)})

if __name__ == "__main__": main()
