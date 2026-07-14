#!/usr/bin/env python3
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import numpy as np, pandas as pd


def load_base():
    path = Path(__file__).with_name("76_run_exhaustive_within_species_environment_rf.py")
    spec = importlib.util.spec_from_file_location("rf76", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def main():
    p=argparse.ArgumentParser(description="RF nonlinear prediction baseline for the primary hierarchical SPDE-INLA analysis.")
    p.add_argument("--environment", required=True); p.add_argument("--trait", required=True); p.add_argument("--out-dir", required=True)
    p.add_argument("--n-splits", type=int, default=5); p.add_argument("--n-estimators", type=int, default=300); p.add_argument("--min-samples-leaf", type=int, default=5)
    p.add_argument("--random-state", type=int, default=20260713); p.add_argument("--permutation-max-n", type=int, default=2000); p.add_argument("--permutation-repeats", type=int, default=1)
    args=p.parse_args(); mod=load_base()
    if args.trait not in mod.TRAITS: raise SystemExit(f"Unknown trait: {args.trait}")
    out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True); df=pd.read_csv(args.environment,low_memory=False)
    climate=list(mod.CLIMATE); topo=list(mod.TOPOGRAPHY); soil=[f"soil_{p}_0_30cm" for p in mod.SOIL_PROPERTIES]
    groups={"climate":climate,"climate_topography":climate+topo,"climate_soil":climate+soil,"climate_topography_soil":climate+topo+soil}
    fold_rows=[]; imp_rows=[]
    for name,preds in groups.items():
        rows,imps=mod.fit_group_cv(df,args.trait,name,preds,args); fold_rows.extend(rows); imp_rows.extend(imps)
    folds=pd.DataFrame(fold_rows); folds.to_csv(out/"rf_spatial_cv_folds.csv",index=False,encoding="utf-8-sig")
    summary=mod.summarize_cv(folds); summary.to_csv(out/"rf_spatial_cv_summary.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame(imp_rows).to_csv(out/"rf_full_model_permutation_importance.csv",index=False,encoding="utf-8-sig")
    if len(summary)!=4: raise SystemExit(f"Expected 4 model-group summaries, got {len(summary)}")
    wide=summary.set_index("model_group")["mean_weighted_r2"]; base=wide.get("climate",np.nan)
    ab=pd.DataFrame([{"trait":args.trait,"analysis_role":"nonlinear_prediction_baseline","r2_climate":base,"r2_climate_topography":wide.get("climate_topography",np.nan),"r2_climate_soil":wide.get("climate_soil",np.nan),"r2_full":wide.get("climate_topography_soil",np.nan),"delta_topography_over_climate":wide.get("climate_topography",np.nan)-base,"delta_soil_over_climate":wide.get("climate_soil",np.nan)-base,"delta_full_over_climate":wide.get("climate_topography_soil",np.nan)-base}])
    ab.to_csv(out/"rf_environment_group_ablation.csv",index=False,encoding="utf-8-sig")
    (out/"trait_rf_report.json").write_text(json.dumps({"trait":args.trait,"analysis_role":"nonlinear_prediction_baseline","primary_analysis":"hierarchical_spde_inla","rows":len(summary)},indent=2),encoding="utf-8")
    print(summary.to_string(index=False)); print(ab.to_string(index=False))

if __name__=="__main__": main()
