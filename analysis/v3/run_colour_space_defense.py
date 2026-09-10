#!/usr/bin/env python3
"""Post-hoc colour-space defense for the frozen v2 chroma-radiation anchor.

This does not redefine the v2 multiplicity family or headline result. It asks
whether the negative chroma-radiation association is accompanied by lower L*
(darkening) or instead represents desaturation/hue reorganization.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analysis.v3 import run_biological_axis_reanalysis as axis

MEMBERS=["corolla_lab_lightness","corolla_lab_chroma","corolla_hue_sin","corolla_hue_cos"]
PRED="chelsa_rsds_mean"

def z(x):
    a=np.asarray(x,float); s=a.std(ddof=0)
    if not np.isfinite(a).all() or s<=0: raise ValueError("no finite variation")
    return (a-a.mean())/s

def perm_scalar(y,x,perms,rng):
    yz=z(y); xz=z(x); den=float(xz@xz); obs=float(xz@yz/den); exc=0
    for _ in range(perms):
        xp=rng.permutation(xz); sim=float(xp@yz/den); exc += abs(sim)>=abs(obs)-1e-15
    return obs,(exc+1)/(perms+1)

def perm_joint(Y,x,perms,rng):
    Yz=np.column_stack([z(Y[:,j]) for j in range(Y.shape[1])]); xz=z(x); den=float(xz@xz)
    beta=xz@Yz/den; mag=float(np.linalg.norm(beta)); exc=0
    for _ in range(perms):
        xp=rng.permutation(xz); sim=xp@Yz/den; exc += np.linalg.norm(sim)>=mag-1e-15
    return beta,mag,(exc+1)/(perms+1)

def bootstrap(frame,reps,seed):
    rng=np.random.default_rng(seed); n=len(frame); rows=[]
    for r in range(reps):
        b=frame.iloc[rng.integers(0,n,n)]
        x=b[PRED].to_numpy(float)
        try:
            bl,_=perm_scalar(b["corolla_lab_lightness"],x,0,rng)
            bc,_=perm_scalar(b["corolla_lab_chroma"],x,0,rng)
        except Exception:
            continue
        rows.append({"replicate":r,"beta_L":bl,"beta_C":bc,"both_negative":bool(bl<0 and bc<0)})
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--traits",type=Path,required=True); ap.add_argument("--environment",type=Path,required=True); ap.add_argument("--out-dir",type=Path,required=True); ap.add_argument("--permutations",type=int,default=9999); ap.add_argument("--bootstrap",type=int,default=2000); ap.add_argument("--seed",type=int,default=20260910)
    a=ap.parse_args(); a.out_dir.mkdir(parents=True,exist_ok=True)
    t,e=axis.load_data(a.traits,a.environment)
    med,_=axis.taxon_endpoint_table(t,MEMBERS,5)
    env=e.groupby("taxon_name")[axis.PREDICTORS].median()
    f=med.join(env,how="inner").dropna(subset=[PRED]).copy()
    if len(f)<50: raise SystemExit("joint colour cohort too small")
    x=f[PRED].to_numpy(float)
    rng=np.random.default_rng(a.seed)
    scalars={}
    for m in ["corolla_lab_lightness","corolla_lab_chroma"]:
        b,p=perm_scalar(f[m].to_numpy(float),x,a.permutations,rng); scalars[m]={"beta_std":b,"p_perm":p}
    beta_h,mag_h,p_h=perm_joint(f[["corolla_hue_sin","corolla_hue_cos"]].to_numpy(float),x,a.permutations,rng)
    beta4,mag4,p4=perm_joint(f[MEMBERS].to_numpy(float),x,a.permutations,rng)
    boot=bootstrap(f,a.bootstrap,a.seed+1)
    boot.to_csv(a.out_dir/"colour_space_taxon_bootstrap.csv",index=False)
    f.reset_index()[["taxon_name",*MEMBERS,PRED]].to_csv(a.out_dir/"colour_space_joint_cohort.csv",index=False)
    summary={
      "analysis_id":"ch1_v3_colour_space_defense_20260910",
      "claim_boundary":"post-hoc interpretation defense only; frozen v2 chroma-radiation result and multiplicity family unchanged",
      "n_taxa_joint_colour_cohort":int(len(f)),
      "radiation_component_slopes":scalars,
      "hue_joint":{"component_betas":{"sin":float(beta_h[0]),"cos":float(beta_h[1])},"effect_magnitude":mag_h,"p_perm":p_h},
      "joint_L_C_hue":{"component_betas":{m:float(v) for m,v in zip(MEMBERS,beta4)},"effect_magnitude":mag4,"p_perm":p4},
      "bootstrap":{"replicates":int(len(boot)),"beta_L_median":float(boot.beta_L.median()),"beta_L_low95":float(boot.beta_L.quantile(.025)),"beta_L_high95":float(boot.beta_L.quantile(.975)),"beta_C_median":float(boot.beta_C.median()),"beta_C_low95":float(boot.beta_C.quantile(.025)),"beta_C_high95":float(boot.beta_C.quantile(.975)),"probability_L_negative":float((boot.beta_L<0).mean()),"probability_C_negative":float((boot.beta_C<0).mean()),"probability_both_negative":float(boot.both_negative.mean())},
      "interpretation_rule":"C*<0 alone means lower chroma, not darkness. Darkening requires an accompanying negative L* displacement; hue components describe colour-direction reorganization but are not biochemical anthocyanin measurements."
    }
    (a.out_dir/"colour_space_defense_report.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2))
if __name__=="__main__": main()
