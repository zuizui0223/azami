"""V5 annual five-axis sensitivity to photographed-angle measurement error.

Reconstruct the frozen V5 native analysis using original source-verified
157,822 photographic measurement records. S1 source stays private. Only
per-scenario aggregates may leave GitHub Actions.

Scenarios are not empirical error estimates. The mirror and bbox-shift
p95 calibration values come from an earlier independent technical audit,
not per-observation V5 error labels. A rainfall-linked offset is a
synthetic differential-bias challenge, never an estimated camera-roll trend.

Primary estimands:
- within-taxon and cross-taxon analogue excess distance under the same
  24 rarefactions, 499 intact-vector phenotype permutations, 1° and 5°;
- within-taxon-pair/all-eligible-cell signed BIO12->image-angle coefficient,
  with original taxon bootstrap (pointwise interval).

Angle deletion gives a different 2D distance estimand, and must not be
read as an ordinary change of the original 3D effect magnitude.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v4.s1_climate_trait_hypervolume import (
    geo_aggregate, standardise_cells
)
from analysis.v4 import s1_h3_climate_analogue_alignment as h

ANNUAL = ("chelsa_bio04", "chelsa_bio15")
ANGLE = "orientation_image_vertical_angle"
RAIN = "chelsa_bio12"
EXPECTED_N = 157822
GRAIN = {1:(63,7117),5:(31,814)}
THRESHOLD = h.CLIMATE_ANALOG_Z * np.sqrt(5 / 3)
N_PERM = 499
N_DRAWS = 24
FROZEN_WITHIN = {1:0.08421227833291306,5:0.23836924606855578}
FROZEN_CROSS = {1:0.030937948818972237,5:0.06532956260499098}
FROZEN_BETA = {1:0.1171,5:0.1829}

def analyse(frame:pd.DataFrame, degrees:int, *, ignore_angle:bool=False)->dict:
    cells,_=geo_aggregate(frame,degrees,ANNUAL)
    z,_=standardise_cells(cells,ANNUAL)
    taxa=h.eligible_taxa(z)
    assert (len(taxa),len(cells))==GRAIN[degrees]
    for t,(_,g) in zip(taxa,z.groupby("accepted_key",sort=True)):
        assert len(t["x"])==len(g)
        t["x"]=np.column_stack([
            t["x"],g[["_z6","_z7"]].to_numpy(float)])
        if ignore_angle:
            t["y"]=t["y"][:,1:]
    rng=np.random.default_rng(20261109)
    draw_rng=np.random.default_rng(20261109)
    perm_rng=np.random.default_rng(20261111)
    streams=np.random.default_rng(20261110)
    maps=[
        np.array([perm_rng.permutation(len(t["y"])) for _ in range(N_PERM)])
        for t in taxa
    ]
    mult=np.random.default_rng(20261112).multinomial(
        len(taxa),np.full(len(taxa),1/len(taxa)),size=N_PERM)
    mult=np.vstack([np.ones(len(taxa)),mult])
    u,v=np.triu_indices(8,1)
    cross_scores=[]
    within=[]
    grams=[]
    for draw in range(N_DRAWS):
        # Original paired climate/phenotype picks and original rarefaction.
        x,y,sector=h._sample_taxa(taxa,rng)
        picks=[draw_rng.choice(len(t["x"]),8,replace=False) for t in taxa]
        if not np.allclose(y,np.stack([t["y"][ix] for t,ix in zip(taxa,picks)])):
            raise ValueError("V5_PHENO_SAMPLING_DRIFT")
        # h.matched_pairs expects 3D phenotype vectors. For angle-deleted
        # scenarios, supply original 3D coordinates *only for the climate mask*,
        # then replace d with actual 2D Euclidean phenotype distances.
        if ignore_angle:
            y3=np.stack([
                z.loc[z.accepted_key.eq(t["taxon"]),["_z3","_z4","_z5"]].to_numpy(float)[ix]
                for t,ix in zip(taxa,picks)
            ])
            pair,support=h.matched_pairs(x,y3,sector,THRESHOLD)
            if pair is not None:
                pair=dict(pair,d=np.linalg.norm(
                    y[pair["i"],:,None,:]-y[pair["j"],None,:,:],
                    axis=3))
        else:
            pair,support=h.matched_pairs(x,y,sector,THRESHOLD)
        if support["status"]!="qualified":
            continue
        seed=int(streams.integers(2**31))
        score,perms=h.observed_and_permutations(
            pair,len(taxa),np.random.default_rng(seed),N_PERM)
        cross_scores.append([score,*perms])
        # Within-taxon sampling and permutations reproduce signed_pairs.py
        # exactly, including its original condition on qualified cross pairs.
        near=np.linalg.norm(x[:,u]-x[:,v],axis=2)<=THRESHOLD
        ok=(near.sum(axis=1)>=2)&(near.sum(axis=1)<len(u))
        if not ok.any():
            continue
        y_ok=y[ok]; mask=near[ok]
        dist=np.linalg.norm(y_ok[:,u]-y_ok[:,v],axis=2)
        obs=float(np.mean(dist.mean(axis=1)-
            (dist*mask).sum(axis=1)/mask.sum(axis=1)))
        null=np.zeros(N_PERM,float)
        for tid in np.flatnonzero(ok):
            indexes=maps[tid][:,picks[tid]]
            shuffled=taxa[tid]["y"][indexes]
            d=np.linalg.norm(shuffled[:,u]-shuffled[:,v],axis=2)
            null+=d.mean(axis=1)-d[:,near[tid]].mean(axis=1)
        within.append([int(ok.sum()),obs,*(null/ok.sum())])
        if ignore_angle:
            continue
        # Source V5 all-cells/within-pair 5-axis linear contrast.
        i,j=pair["i"],pair["j"]
        dx=(x[i,:,None,:]-x[j,None,:,:]).reshape(len(i),64,5)
        dy=(y[i,:,None,:]-y[j,None,:,:]).reshape(len(i),64,3)
        weights=np.full((len(i),64),1/64.0)
        mx=(dx*weights[:,:,None]).sum(axis=1,keepdims=True)
        my=(dy*weights[:,:,None]).sum(axis=1,keepdims=True)
        rx=dx-mx;ry=dy-my
        xx=np.einsum("nka,nkb,nk->nab",rx,rx,weights)
        xy=np.einsum("nka,nkb,nk->nab",rx,ry,weights)
        bw=mult[:,i]*mult[:,j]
        bw=bw/bw.sum(axis=1,keepdims=True)
        grams.append((
            bw@xx.reshape(-1,25),
            bw@xy.reshape(-1,15)
        ))
    if len(within)!=N_DRAWS or len(cross_scores)!=N_DRAWS:
        raise ValueError(f"V5_DRAW_FAMILY_DRIFT_{degrees}: "+
                         f"{len(within)}/{len(cross_scores)}")
    obs,*perms=np.mean(within,axis=0)
    null=np.asarray(perms)
    base=np.mean(null)
    within_excess=float(obs-base)
    within_p=(1+np.count_nonzero(
        np.abs(null-base)>=abs(within_excess)))/(N_PERM+1)
    obs,*perms=np.mean(cross_scores,axis=0)
    null=np.asarray(perms)
    base=np.mean(null)
    cross_excess=float(obs-base)
    cross_p=(1+np.count_nonzero(
        np.abs(null-base)>=abs(cross_excess)-1e-12))/(N_PERM+1)
    result={
        "degrees":degrees,"taxa":len(taxa),
        "taxon_cells":len(cells),"draws":len(within),
        "dimension":2 if ignore_angle else 3,
        "within_excess":within_excess,"within_p":float(within_p),
        "cross_excess":cross_excess,"cross_p":float(cross_p),
    }
    if grams:
        xx=np.mean([g[0] for g in grams],axis=0).reshape(-1,5,5)
        xy=np.mean([g[1] for g in grams],axis=0).reshape(-1,5,3)
        coeff=np.linalg.solve(xx,xy)
        lo,hi=np.quantile(coeff[1:],[.025,.975],axis=0)
        result.update({
            "BIO12_to_image_angle_beta":float(coeff[0,0,0]),
            "BIO12_to_image_angle_pointwise_lo":float(lo[0,0]),
            "BIO12_to_image_angle_pointwise_hi":float(hi[0,0]),
            "BIO15_to_chroma_beta":float(coeff[0,4,1]),
        })
    return result

def run(selected:pd.DataFrame, out:Path)->dict:
    assert len(selected)==EXPECTED_N and selected.obs_id.is_unique
    original=selected.copy()
    orig=selected[ANGLE].to_numpy(float)
    rain=selected[RAIN].to_numpy(float)
    if not (np.isfinite(orig).all() and (orig>=0).all() and
            (orig<=180).all() and np.isfinite(rain).all()):
        raise ValueError("ANGLE_OR_RAIN_INVALID")
    # Within-source-taxon rainfall standardized offsets.
    groups=selected.groupby("accepted_key",sort=False)[RAIN]
    means=groups.transform("mean").to_numpy(float)
    sds=groups.transform("std").to_numpy(float)
    z=np.clip((rain-means)/np.maximum(sds,1e-8),-2,2)
    scenario_specs=[
        ("baseline_3trait", None,None,False),
        ("without_angle_2trait",None,None,True),
        ("random_mirror_p95_4p67",4.67/1.95996398454,0,False),
        ("random_bbox_5pct_p95_54p1",54.1/1.95996398454,0,False),
        ("rain_correlated_plus10deg_per_within_taxon_SD",None,+10,False),
        ("rain_correlated_minus10deg_per_within_taxon_SD",None,-10,False),
    ]
    scenarios=[]
    for name,sigma,bias,drop in scenario_specs:
        reps=5 if sigma is not None else 1
        for replicate in range(reps):
            if sigma is not None:
                rand=np.random.default_rng(
                    2026112000+replicate+(0 if sigma<10 else 100))
                angle=np.clip(orig+rand.normal(0,sigma,len(orig)),0,180)
                frame=original.assign(**{ANGLE:angle})
            elif bias is not None:
                frame=original.assign(**{ANGLE:np.clip(orig+bias*z,0,180)})
            else:
                frame=original
            for degree in (1,5):
                result=analyse(frame,degree,ignore_angle=drop)
                result["scenario"]=name
                result["replicate"]=replicate
                scenarios.append(result)
                if name=="baseline_3trait":
                    if not np.isclose(result["within_excess"],
                                     FROZEN_WITHIN[degree],atol=1e-7):
                        raise ValueError(f"ORIGINAL_V5_WITHIN_NOT_REPRODUCED_{degree}: "+
                                         str(result["within_excess"]))
                    if not np.isclose(result["cross_excess"],
                                     FROZEN_CROSS[degree],atol=1e-7):
                        raise ValueError(f"ORIGINAL_V5_CROSS_NOT_REPRODUCED_{degree}: "+
                                         str(result["cross_excess"]))
                    if not np.isclose(result["BIO12_to_image_angle_beta"],
                                     FROZEN_BETA[degree],atol=.0001):
                        raise ValueError(f"ORIGINAL_V5_BIO12_NOT_REPRODUCED_{degree}: "+
                                         str(result["BIO12_to_image_angle_beta"]))
    df=pd.DataFrame(scenarios)
    # Original correction family is the two spatial grains, separately
    # for within-taxon and between-taxon tests. Never count perturbation
    # replicates as distinct original biological data.
    for (scenario,rep),idx in df.groupby(["scenario","replicate"]).groups.items():
        for name in ("within","cross"):
            df.loc[idx,name+"_holm"]=multipletests(
                df.loc[idx,name+"_p"],method="holm")[1]
    exported=df.to_dict("records")
    summarised=[]
    for (name,degree),p in df.groupby(["scenario","degrees"],sort=False):
        row={
            "scenario":name,"degrees":int(degree),
            "simulations":len(p),"dimensionality":int(p.dimension.iloc[0]),
        }
        for key in ("within_excess","within_holm","cross_excess","cross_holm",
                    "BIO12_to_image_angle_beta",
                    "BIO12_to_image_angle_pointwise_lo",
                    "BIO12_to_image_angle_pointwise_hi"):
            z=pd.to_numeric(p[key],errors="coerce").dropna()
            row[key+"_median"]=float(z.median()) if len(z) else None
            row[key+"_min"]=float(z.min()) if len(z) else None
            row[key+"_max"]=float(z.max()) if len(z) else None
        summarised.append(row)
    receipt={
      "schema":"v5_annual_5axis_angle_error_perturbation_v1",
      "source_cohort_records":int(len(selected)),
      "source_definition":"exact original 157822 V5 source photo-derived observations",
      "baseline_frozen_reproduced":True,
      "errors_are_synthetic":True,
      "mirror_reference_p95_degrees":4.67,
      "five_percent_bbox_shift_reference_p95_degrees":54.1,
      "random_p95_convention":"sigma = reference absolute-error p95 / 1.95996, iid zero-mean Gaussian record-level perturbation",
      "bias_convention":"±10° per within-taxon BIO12 SD, z clipped to ±2; hypothetical, not estimated",
      "angle_removed":"two-dimensional trait estimand, not same-scale 3D effect; unchanged source record selection",
      "inference_warning":"Neither random scenarios nor marginally significant coefficients validate gravity, photography roll, detector bias, trait evolution or causality",
      "results":summarised,
      "per_replicate":exported,
    }
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",
                   encoding="utf-8")
    print("V5_ANGLE_ERROR_SENSITIVITY_JSON_START")
    print(json.dumps(summarised,ensure_ascii=False,separators=(",",":")))
    print("V5_ANGLE_ERROR_SENSITIVITY_JSON_END")
    return receipt
