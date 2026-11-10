"""V5 Part I: annual vs warm-quarter climate + both seasonality axes.

Same source photos, taxa, cells, rarefaction, null draws and signed-pair
taxonomic bootstrap. Does NOT use climate rarity. Not phylogenetic or causal.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v4.s1_climate_trait_hypervolume import (
    AXES, PHENOTYPE, load_source, geo_aggregate, standardise_cells
)
from analysis.v4 import s1_h3_climate_analogue_alignment as h

EXTRA=("chelsa_bio04","chelsa_bio15","chelsa_bio18","chelsa_bio05")
CLIMATE={
    "annual_plus_seasonality":("_z0","_z1","_z2","_z6","_z7"),
    "warm_plus_seasonality":("_z8","_z9","_z2","_z6","_z7")
}
CLIMATE_NAMES={
    "annual_plus_seasonality":("BIO12","BIO01","RSDS","BIO04","BIO15"),
    "warm_plus_seasonality":("BIO18","BIO05","RSDS","BIO04","BIO15")
}
TRAIT_NAMES=("angle","chroma","outline")
N_DRAWS=24
N_NULL=199
N_BOOT=199
THRESH=h.CLIMATE_ANALOG_Z*np.sqrt(5/3)


def verified_photo_source(native, decoded, receipt):
    """Exactly reproduce frozen 157,822 S1 photo selection *before* warm QC."""
    ref,proof=load_source(decoded,receipt,EXTRA[:2])
    require=[*AXES,*EXTRA[:2],"latitude","longitude"]
    mask=np.isfinite(native[require].to_numpy(float)).all(axis=1)
    for field in PHENOTYPE:
        mask &= native[field+"__eligible"].fillna(False).eq(True).to_numpy(bool)
    mask &= native[PHENOTYPE[0]].between(0,180).to_numpy(bool)
    sub=native.loc[mask,["accepted_key","latitude","longitude",*AXES,*EXTRA]].copy()
    if len(ref)!=len(sub):
        raise ValueError("S1_ANNUAL_SOURCE_PHOTO_MASK_DRIFT")
    np.testing.assert_array_equal(ref.accepted_key.to_numpy(),sub.accepted_key.to_numpy())
    np.testing.assert_allclose(ref[["latitude","longitude",*AXES,*EXTRA[:2]]].to_numpy(float),
                               sub[["latitude","longitude",*AXES,*EXTRA[:2]]].to_numpy(float),rtol=0,atol=1e-10)
    warm_mask=np.isfinite(sub[list(EXTRA[2:])].to_numpy(float)).all(axis=1)
    included=sub.loc[warm_mask].copy()
    if len(included)<50000:
        raise ValueError("S1_WARM_COMMON_COHORT_TOO_SMALL")
    return included,{
        "verified_annual_source_photo_n":len(ref),
        "matched_annual_and_warm_source_photo_n":len(included),
        "warm_source_excluded_original_eligible":int((~warm_mask).sum()),
        "photo_trait_eligible_same_mask_both_climates":True,
        "original_decoded_source_sha":proof["decoded_sha256"],
    }


def _fixed_taxa(z):
    taxa=[]
    for key,g in z.groupby("accepted_key",sort=True):
        if len(g)<8:continue
        sec=h.geographic_sector(g.longitude.to_numpy(float))
        sector=int(np.bincount(sec,minlength=3).argmax())
        taxa.append({
            "sector":sector,"y":g[["_z3","_z4","_z5"]].to_numpy(float),
            "x":{name:g[list(cols)].to_numpy(float) for name,cols in CLIMATE.items()}
        })
    if len(taxa)<25:
        raise ValueError("S1_INSUFFICIENT_TAXA")
    return taxa


def _bootstrap_moments(x,y,i,j,mult):
    """One vote per cross-sector taxon pair with taxon-block bootstrap."""
    dx=(x[i,:,None,:]-x[j,None,:,:]).reshape(len(i),64,5)
    dy=(y[i,:,None,:]-y[j,None,:,:]).reshape(len(i),64,3)
    mx=dx.mean(axis=1)
    my=dy.mean(axis=1)
    dwx=dx-mx[:,None,:]
    dwy=dy-my[:,None,:]
    w=mult[:,i]*mult[:,j]
    if (w.sum(axis=1)<=0).any():
        raise ValueError("TAXON_BLOCK_BOOTSTRAP_SECTOR_MASS_ZERO")
    w=w/w.sum(axis=1,keepdims=True)
    mats={}
    for typ,vx,vy in (("between_pair",mx,my),("within_pair",dwx,dwy)):
        if typ=="between_pair":
            xx=np.einsum("pa,pb->pab",vx,vx)
            xy=np.einsum("pa,pb->pab",vx,vy)
        else:
            xx=np.einsum("pna,pnb->pab",vx,vx)/64
            xy=np.einsum("pna,pnb->pab",vx,vy)/64
        mats[typ]=(w@xx.reshape(len(i),25),w@xy.reshape(len(i),15))
    return mats


def run(native,decoded,receipt):
    d,proof=verified_photo_source(native,decoded,receipt)
    results={"status":"SOURCE_MATCHED_V5_S1_ANNUAL_VS_WARM_5AXIS",
             "cohort":proof,"models":list(CLIMATE_NAMES),
             "model_axis_names":{k:list(v) for k,v in CLIMATE_NAMES.items()},
             "climate_standardised_with_equal_total_weight_per_taxon":True,
             "same_photos_cells_taxa_rarefactions":True,
             "s1_cross_sector_permuted_null":True,
             "model_comparison_exploratory_after_previous_V5_results":True,
             "not_adaptation_or_phylogenetic_convergence":True,
             "resolution":[],"cross":[],"within":[],"directional":[],"paired_model_differences":[]}
    cross_all=[];within_all=[]
    for deg in (1,5):
        cells,coverage=geo_aggregate(d,deg,EXTRA)
        z,_=standardise_cells(cells,EXTRA)
        taxa=_fixed_taxa(z)
        n=len(taxa)
        sectors=np.array([t["sector"] for t in taxa])
        i,j=np.where(np.triu(sectors[:,None]!=sectors[None,:],1))
        if len(i)<40:
            raise ValueError("S1_INSUFFICIENT_FIXED_CROSS_SECTOR_TAXON_PAIRS")
        boot_rng=np.random.default_rng(20261112+deg)
        multinomial=boot_rng.multinomial(n,np.full(n,1/n),size=N_BOOT)
        mult=np.vstack([np.ones(n),multinomial])
        rng=np.random.default_rng(20261109+deg)
        summaries={model:{"cross_o":[],"cross_n":[],"cross_pairs":[],"within_o":[],"within_n":[],"within_n_taxa":[],
                          "moments":{key:[] for key in ("between_pair","within_pair")}}
                   for model in CLIMATE}
        u,v=np.triu_indices(8,1)
        for rep in range(N_DRAWS):
            picks=[rng.choice(len(t["y"]),8,replace=False) for t in taxa]
            y=np.stack([taxa[k]["y"][pick] for k,pick in enumerate(picks)])
            # Identical phenotype permutations for both models in this draw.
            permrng=np.random.default_rng(20261111+deg*1000+rep)
            permutations=[np.array([permrng.permutation(8) for _ in range(N_NULL)])
                          for _ in taxa]
            for name in CLIMATE:
                x=np.stack([taxa[k]["x"][name][pick] for k,pick in enumerate(picks)])
                summary=summaries[name]
                pairs,support=h.matched_pairs(x,y,sectors,THRESH)
                if support["status"]=="qualified":
                    obs,null=h.observed_and_permutations(
                        pairs,n,np.random.default_rng(20261119+deg*1000+rep),N_NULL)
                    summary["cross_o"].append(obs)
                    summary["cross_n"].append(null)
                    summary["cross_pairs"].append(support["n_pairs"])
                near=np.linalg.norm(x[:,u]-x[:,v],axis=2)<=THRESH
                usable=(near.sum(1)>=2)&(near.sum(1)<len(u))
                if usable.any():
                    dist=np.linalg.norm(y[:,u]-y[:,v],axis=2)
                    obs=np.mean((dist.mean(1)-(dist*near).sum(1)/np.maximum(near.sum(1),1))[usable])
                    null=np.zeros(N_NULL)
                    for t in np.flatnonzero(usable):
                        yp=y[t][permutations[t]]
                        shuffled=np.linalg.norm(yp[:,u]-yp[:,v],axis=2)
                        null+=shuffled.mean(1)-(shuffled*near[t]).sum(1)/near[t].sum()
                    summary["within_o"].append(obs)
                    summary["within_n"].append(null/usable.sum())
                    summary["within_n_taxa"].append(int(usable.sum()))
                moments=_bootstrap_moments(x,y,i,j,mult)
                for k in moments:
                    summary["moments"][k].append(moments[k])
        results["resolution"].append({
            "degrees":deg,"taxa":n,"taxon_cells":len(cells),
            "eligible_source_photos":len(d),"cross_sector_taxon_pairs_available":len(i)})
        for model,values in summaries.items():
            if len(values["cross_o"])<12 or len(values["within_o"])<12:
                raise ValueError("S1_INSUFFICIENT_CLIMATE_ANALOG_DRAW_"+str(deg)+"_"+model)
            for typ,prefix in (("cross","cross"),("within","within")):
                obs=float(np.mean(values[prefix+"_o"]))
                null=np.mean(values[prefix+"_n"],axis=0)
                center=float(np.mean(null))
                excess=obs-center
                p=float((1+(np.abs(null-center)>=abs(excess)-1e-12).sum())/(N_NULL+1))
                rec={
                    "degrees":deg,"model":model,"n_taxa":n,"draws":len(values[prefix+"_o"]),
                    "observed":obs,"permutation_null_mean":center,
                    "excess_trait_similarity":float(excess),"p":p,
                    "photo_mask_identical_both_models":True,
                    "climate_threshold_RMS_5axis_z":THRESH}
                if typ=="cross":rec["mean_taxon_pairs"]=float(np.mean(values["cross_pairs"]))
                else:rec["mean_supported_taxa"]=float(np.mean(values["within_n_taxa"]))
                results[typ].append(rec)
            for scale,repetitions in values["moments"].items():
                xx=np.mean([q[0] for q in repetitions],axis=0).reshape(N_BOOT+1,5,5)
                xy=np.mean([q[1] for q in repetitions],axis=0).reshape(N_BOOT+1,5,3)
                if np.linalg.matrix_rank(xx[0])!=5:
                    raise ValueError("S1_SINGULAR_SIGNED_PAIR_"+str(deg)+"_"+model)
                betas=np.linalg.solve(xx,xy)
                for idx,axis in enumerate(CLIMATE_NAMES[model]):
                    for tid,trait in enumerate(TRAIT_NAMES):
                        lo,hi=np.quantile(betas[1:,idx,tid],[.025,.975])
                        results["directional"].append({
                            "degrees":deg,"model":model,"scale":scale,
                            "axis":axis,"trait":trait,"standardised_beta":float(betas[0,idx,tid]),
                            "bootstrap_ci_pointwise_95":[float(lo),float(hi)],
                            "taxon_bootstraps":N_BOOT,"draws":len(repetitions),
                            "cross_sector_pairs_fixed_across_models":True,
                            "climate_matching_not_used_to_select_signed_pairs":True})
                summaries[model]["moments"][scale]=betas
        for scale in ("within_pair","between_pair"):
            annual=summaries["annual_plus_seasonality"]["moments"][scale]
            warm=summaries["warm_plus_seasonality"]["moments"][scale]
            # Corresponding positions are rain, temperature, RSDS, temp
            # seasonality, rain seasonality. Both are in each cohort's 1-SD units.
            delta=warm-annual
            for idx,name in enumerate(("rain","temperature","radiation","temperature_seasonality","precipitation_seasonality")):
                for jtrait,trait in enumerate(TRAIT_NAMES):
                    lo,hi=np.quantile(delta[1:,idx,jtrait],[.025,.975])
                    results["paired_model_differences"].append({
                        "degrees":deg,"scale":scale,"axis_family":name,"trait":trait,
                        "warm_minus_annual_std_beta":float(delta[0,idx,jtrait]),
                        "bootstrap_ci_pointwise_95":[float(lo),float(hi)],
                        "shared_source_taxa_and_cells":True})
    for kind in ("cross","within"):
        tab=results[kind]
        p=multipletests([r["p"] for r in tab],method="holm")[1]
        for row,v in zip(tab,p):
            row["holm_4_models_x_grains"]=float(v)
    return results
