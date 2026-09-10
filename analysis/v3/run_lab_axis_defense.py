#!/usr/bin/env python3
"""Post-hoc L*a*b* interpretation defense for the frozen v2 chroma-radiation anchor.

This keeps the frozen v2 result unchanged. It reconstructs taxon-level a* and b*
from the same 143-taxon L*/C*/hue cohort used by run_colour_space_defense.py,
then asks whether radiation is associated with darkness (L*), red-magenta
expression (a*) or blue-yellow direction (b*).
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v3 import run_biological_axis_reanalysis as axis

MEMBERS=["corolla_lab_lightness","corolla_lab_chroma","corolla_hue_sin","corolla_hue_cos"]
PRED="chelsa_rsds_mean"

def z(x):
    a=np.asarray(x,float); s=a.std(ddof=0)
    if not np.isfinite(a).all() or s<=0: raise ValueError("no finite variation")
    return (a-a.mean())/s

def scalar_beta(y,x):
    yz=z(y); xz=z(x); return float(xz@yz/(xz@xz))

def perm_scalar(y,x,perms,rng):
    yz=z(y); xz=z(x); den=float(xz@xz); obs=float(xz@yz/den); exc=0
    for _ in range(perms):
        xp=rng.permutation(xz); sim=float(xp@yz/den); exc += abs(sim)>=abs(obs)-1e-15
    return obs,(exc+1)/(perms+1)

def perm_joint(frame,cols,x,perms,rng):
    Y=np.column_stack([z(frame[c]) for c in cols]); xz=z(x); den=float(xz@xz)
    beta=xz@Y/den; mag=float(np.linalg.norm(beta)); exc=0
    for _ in range(perms):
        xp=rng.permutation(xz); sim=xp@Y/den; exc += np.linalg.norm(sim)>=mag-1e-15
    return beta,mag,(exc+1)/(perms+1)

def bootstrap(frame,reps,seed):
    rng=np.random.default_rng(seed); n=len(frame); rows=[]
    for r in range(reps):
        b=frame.iloc[rng.integers(0,n,n)]
        x=b[PRED].to_numpy(float)
        try:
            vals={c:scalar_beta(b[c],x) for c in ["corolla_lab_lightness","corolla_lab_chroma","lab_a_reconstructed","lab_b_reconstructed"]}
        except ValueError:
            continue
        rows.append({"replicate":r,"beta_L":vals["corolla_lab_lightness"],"beta_C":vals["corolla_lab_chroma"],"beta_a":vals["lab_a_reconstructed"],"beta_b":vals["lab_b_reconstructed"]})
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--traits",type=Path,required=True); ap.add_argument("--environment",type=Path,required=True); ap.add_argument("--out-dir",type=Path,required=True); ap.add_argument("--permutations",type=int,default=9999); ap.add_argument("--bootstrap",type=int,default=2000); ap.add_argument("--seed",type=int,default=20260910)
    a=ap.parse_args(); a.out_dir.mkdir(parents=True,exist_ok=True)
    t,e=axis.load_data(a.traits,a.environment)
    med,_=axis.taxon_endpoint_table(t,MEMBERS,5)
    env=e.groupby("taxon_name")[axis.PREDICTORS].median()
    f=med.join(env,how="inner").dropna(subset=[PRED]).copy()
    if len(f)!=143: raise SystemExit(f"expected exact 143-taxon colour cohort, got {len(f)}")
    h=np.arctan2(f["corolla_hue_sin"].to_numpy(float),f["corolla_hue_cos"].to_numpy(float))
    c=f["corolla_lab_chroma"].to_numpy(float)
    f["lab_a_reconstructed"]=c*np.cos(h)
    f["lab_b_reconstructed"]=c*np.sin(h)
    x=f[PRED].to_numpy(float); rng=np.random.default_rng(a.seed)
    scalars={}
    for col,label in [("corolla_lab_lightness","L_star"),("corolla_lab_chroma","C_star"),("lab_a_reconstructed","a_star"),("lab_b_reconstructed","b_star")]:
        b,p=perm_scalar(f[col],x,a.permutations,rng); scalars[label]={"beta_std":b,"p_perm":p}
    beta3,mag3,p3=perm_joint(f,["corolla_lab_lightness","lab_a_reconstructed","lab_b_reconstructed"],x,a.permutations,rng)
    boot=bootstrap(f,a.bootstrap,a.seed+1); boot.to_csv(a.out_dir/"lab_axis_taxon_bootstrap.csv",index=False)
    f.reset_index()[["taxon_name",*MEMBERS,"lab_a_reconstructed","lab_b_reconstructed",PRED]].to_csv(a.out_dir/"lab_axis_joint_cohort.csv",index=False)
    def bs(col):
        s=boot[col]; return {"median":float(s.median()),"low95":float(s.quantile(.025)),"high95":float(s.quantile(.975)),"probability_negative":float((s<0).mean()),"probability_positive":float((s>0).mean())}
    report={"analysis_id":"ch1_v3_lab_axis_defense_20260910","claim_boundary":"post-hoc interpretation defense only; a* and b* are reconstructed from the frozen taxon-level C* and hue direction and are not new measured endpoints or anthocyanin concentrations","n_taxa":int(len(f)),"radiation_slopes":scalars,"joint_Lab_displacement":{"component_betas":{"L_star":float(beta3[0]),"a_star":float(beta3[1]),"b_star":float(beta3[2])},"effect_magnitude":mag3,"p_perm":p3},"bootstrap":{"replicates":int(len(boot)),"L_star":bs("beta_L"),"C_star":bs("beta_C"),"a_star":bs("beta_a"),"b_star":bs("beta_b")},"interpretation":{"L_star":"higher values are lighter; negative radiation slope would support darkening","a_star":"positive direction is red-magenta; not equivalent to anthocyanin concentration","b_star":"positive direction is yellow and negative direction is blue; shifts help distinguish purple/blue versus warmer colour direction","C_star":"distance from the neutral axis; lower C* means desaturation, not darkness"}}
    (a.out_dir/"lab_axis_defense_report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2))
if __name__=="__main__": main()
