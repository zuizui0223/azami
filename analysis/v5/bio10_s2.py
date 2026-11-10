"""V5 Part II: annual and two alternative warm-quarter climate niches, both with seasonality.

Both models: exactly identical wild source, photographic phenotypes,
regions, random draws, smoothing parameters, and region-label permutation.
No climate rarity scalar and no phenotype-determined climate choices.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_introduced_opportunity import angular_cells
from analysis.v5.native_pca_env import (
    SPECIES_KEYS, _native_and_introduced_codes, _probability_content,
    _weights_per_5degree_cell
)
from analysis.v5.novel_head_overlap import analysis as head_novel_analysis
from analysis.v5.background_paired_overlap import (
    TRAITS, SEED, _head_cells, _region_background,
    block_key, regional_density, correct_occupancy, overlap_matrix,
    project as photo_project, grid as photo_grid, infer
)
from analysis.v5.bio10_s1 import CLIMATE_NAMES
# S1 adds BIO04/BIO15 as extras because its old source had only 3 axes.
# S2 already includes both in VARIABLES: list just newly sampled layers,
# otherwise DataFrame column selection duplicates BIO04/BIO15.
EXTRA=("chelsa_bio18","chelsa_bio05","chelsa_bio10")

MODELS={
    "annual_plus_seasonality":tuple(VARIABLES),
    "warm_plus_seasonality":(
        "chelsa_bio05","chelsa_bio04","chelsa_bio18",
        "chelsa_bio15","chelsa_rsds_mean"),
    "warm_bio10_plus_seasonality":(
        "chelsa_bio10","chelsa_bio04","chelsa_bio18",
        "chelsa_bio15","chelsa_rsds_mean")
}
REFERENCE_ANNUAL_D={"Cirsium arvense":.474,"Cirsium vulgare":.442}
REFERENCE_OLD_RHO={"Cirsium arvense":.07109894,"Cirsium vulgare":-.107866}
HEADER=("range_state","block","i","j")
MINREGIONS=3
DRAWS=32
BSIZE=8


def full_climate_pca(background, axes, normalise_pc=False):
    """Fixed common weighted PCA, half area per domain, like V5 PCA-env."""
    data=background[list(axes)].to_numpy(float)
    wt=background._area_weight.to_numpy(float).copy()
    states=background.range_state.to_numpy(str)
    for label in ("native","introduced"):
        chosen=states==label
        if chosen.sum()<200 or (wt[chosen]<=0).any():
            raise ValueError("INVALID_WCVP_BACKGROUND_AREA")
        wt[chosen]/=(2*wt[chosen].sum())
    wt/=wt.sum()
    center=np.average(data,axis=0,weights=wt)
    scale=np.sqrt(np.average((data-center)**2,axis=0,weights=wt))
    if np.any(scale<=1e-10):
        raise ValueError("ZERO_CLIMATE_BACKGROUND_PCA_AXIS")
    z=(data-center)/scale
    covariance=(z*wt[:,None]).T@z
    eig,vec=np.linalg.eigh(covariance)
    order=np.argsort(eig)[::-1]
    eig=eig[order]
    basis=vec[:,order[:2]].copy()
    for k in range(2):
        j=int(np.argmax(np.abs(basis[:,k])))
        if basis[j,k]<0:basis[:,k]*=-1
    if eig[1]<=1e-10:raise ValueError("CLIMATE_BACKGROUND_PCA_RANK_LT_2")
    def transform(dataframe):
        x=dataframe[list(axes)].to_numpy(float) if isinstance(dataframe,pd.DataFrame) else np.asarray(dataframe,float)
        result=((x-center)/scale)@basis
        if normalise_pc:result=result/np.sqrt(eig[:2])
        return result
    return transform,float(eig[:2].sum()/eig.sum())


def _niche(species, wild, land, authority, model):
    axes=MODELS[model]
    status=_native_and_introduced_codes(authority,SPECIES_KEYS[species])
    env=land.loc[land.TDWG_code.isin(status)].copy()
    env["range_state"]=env.TDWG_code.map(status)
    if len(env)<5000:raise ValueError("INSUFFICIENT_NICHE_LAND")
    complete=env[np.isfinite(env[list(axes)].to_numpy(float)).all(axis=1)].copy()
    if complete._area_weight.sum()/env._area_weight.sum()<.8:
        raise ValueError("FULL_LAND_BACKGROUND_COVERAGE_LT_80PC")
    ni=complete.loc[complete.range_state=="native"]
    ii=complete.loc[complete.range_state=="introduced"]
    occ=wild.loc[wild.taxon_name==species].copy()
    occ["_cell"]=angular_cells(occ.latitude,occ.longitude)
    occ=occ.sort_values("obs_id").drop_duplicates(
        ["taxon_name","range_state","_cell","user_id"]).copy()
    occ=occ.loc[np.isfinite(occ[list(axes)].to_numpy(float)).all(axis=1)].copy()
    no=occ.loc[occ.range_state=="native"]
    io=occ.loc[occ.range_state=="introduced"]
    if min(len(no),len(io))<100:raise ValueError("INSUFFICIENT_NICHE_SOURCE_OBSERVATIONS")
    transform,pc_fraction=full_climate_pca(complete,axes)
    xbg=transform(complete)
    low=np.quantile(xbg,.005,axis=0);high=np.quantile(xbg,.995,axis=0)
    margin=np.maximum((high-low)*.1,.1)
    edges=[np.linspace(low[k]-margin[k],high[k]+margin[k],81) for k in range(2)]
    def distribution(avail,observed):
        b=transform(avail);o=transform(observed)
        aw=avail._area_weight.to_numpy(float)
        ow=_weights_per_5degree_cell(observed)
        hbg=np.histogram2d(b[:,0],b[:,1],bins=edges,weights=aw/aw.sum())[0]
        hoc=np.histogram2d(o[:,0],o[:,1],bins=edges,weights=ow/ow.sum())[0]
        if hbg.sum()<.95 or hoc.sum()<.80:
            raise ValueError("PCA_NICHE_GRID_TRUNCATES_BACKGROUND_OR_OCCUPIED")
        fbg=gaussian_filter(hbg,1.5,mode="constant")
        focc=gaussian_filter(hoc,1.5,mode="constant")
        floor=max(float(fbg.max())*.0001,1e-12)
        d=focc/np.maximum(fbg,floor)
        d[fbg<floor]=0
        if d.sum()<=0:raise ValueError("NICHE_NO_COMMON_PCA_SUPPORT")
        return d/d.sum()
    a=distribution(ni,no)
    b=distribution(ii,io)
    d=float(np.minimum(a,b).sum())
    native95=_probability_content(a,.95)
    novelty=float(b[~native95].sum())
    # The verified V5 annual baseline used the *same* full wild WCVP cohort
    # and exactly these background weights, PC bins and bandwidth.
    if model=="annual_plus_seasonality" and abs(d-REFERENCE_ANNUAL_D[species])>.0015:
        raise ValueError("ANNUAL_CLIMATE_NICHE_BASELINE_NOT_REPRODUCED_"+species)
    row={
        "species":species,"model":model,"climate_axes":list(axes),
        "n_full_wild_native":len(no),"n_full_wild_introduced":len(io),
        "n_background_native":len(ni),"n_background_introduced":len(ii),
        "full_range_climate_overlap_D":d,
        "introduced_corrected_mass_outside_native95":novelty,
        "pc12_background_variance_fraction":pc_fraction,
        "same_wild_source_both_models":True,
        "availability_corrected":True,"pca_env_style_not_ecospat_exact":True,
    }
    return row,(transform,edges,native95),complete


def _photo_cells(photo_wild, measured):
    """Reproduce the exact old V5 photo and 1-degree-cell cohort (no imputation)."""
    frame=photo_wild.drop(columns=TRAITS,errors="ignore").merge(
        measured[["obs_id",*TRAITS]],on="obs_id",how="left",validate="one_to_one")
    mask=frame[TRAITS[0]].between(0,180) & np.isfinite(frame[[*TRAITS,*VARIABLES,*EXTRA]].to_numpy(float)).all(axis=1)
    frame=frame.loc[mask].copy()
    frame["coarse"]=angular_cells(frame.latitude,frame.longitude)
    frame=frame.sort_values("obs_id").drop_duplicates(
        ["taxon_name","range_state","coarse","user_id"]).copy()
    frame["i"]=(frame.latitude+90).clip(upper=179.99999).astype(int)
    frame["j"]=((frame.longitude+180)%360).astype(int)
    frame["block"]=(frame.i//20)*18+frame.j//20
    agg=frame.groupby(["taxon_name",*HEADER],sort=True)
    c=agg[[*TRAITS,*VARIABLES,*EXTRA]].median().reset_index()
    c["photos"]=agg.size().to_numpy()
    c["latitude"]=agg.latitude.median().to_numpy()
    c["longitude"]=agg.longitude.median().to_numpy()
    old=_head_cells(photo_wild)
    for species,b in c.groupby("taxon_name",sort=True):
        original=old[species]
        b=b.reset_index(drop=True)
        for var in [*HEADER,*TRAITS,*VARIABLES]:
            if var in HEADER:
                np.testing.assert_array_equal(b[var].to_numpy(),original[var].to_numpy())
            else:
                np.testing.assert_allclose(b[var].to_numpy(float),original[var].to_numpy(float),atol=1e-12,rtol=0)
        if int(b.photos.sum())!=int(original.photos.sum()):
            raise ValueError("S2_ORIGINAL_PHOTO_COHORT_CHANGED")
    if set(c.taxon_name)!=set(SPECIES_KEYS):
        raise ValueError("S2_TWO_SPECIES_COHORT_MISSING")
    return frame,c


def _regional_pair(c,env,model,species,niche_info):
    """Replicate successful region-pair D with arbitrary matched climate axes."""
    axes=MODELS[model]
    work=c.reset_index(drop=True)
    tr,var=full_climate_pca(env,axes,normalise_pc=True)
    xb=tr(env)
    xp=tr(work)
    coords=np.vstack([xb,xp])
    edges=[np.linspace(coords[:,j].min()-.5,coords[:,j].max()+.5,65) for j in range(2)]
    xold,_=photo_project(work,list(axes))
    y3,_=photo_project(work,TRAITS)
    y2,_=photo_project(work,TRAITS[:2])
    head_edges={3:photo_grid(y3),2:photo_grid(y2)}
    photo_edges=photo_grid(xold)
    work=work.copy();work["_idx"]=np.arange(len(work))
    reg=_region_background(env,work)
    keys={}
    for side in ("native","introduced"):
        keys[side]=sorted((side,int(block)) for block in work[work.range_state==side].block.unique()
                          if (side,int(block)) in reg)
    if min(map(len,keys.values()))<MINREGIONS:
        raise ValueError("REGIONAL_BACKGROUND_INADEQUATE_"+species+"_"+model)
    index={k:work.loc[(work.range_state==k[0])&(work.block==k[1]),"_idx"].to_numpy(int)
           for k in keys["native"]+keys["introduced"]}
    result=[]
    for h in (.35,.5,.75):
        backgrounds={k:regional_density(tr(reg[k]["available"]),edges,h,
                                        weights=reg[k]["available"]._area_weight.to_numpy(float))
                     for k in index}
        for dims,traitpoints in ((3,y3),(2,y2)):
            rng=np.random.default_rng(SEED)
            corrected=np.zeros((len(keys["native"]),len(keys["introduced"])))
            uncorrected=np.zeros_like(corrected)
            heads=np.zeros_like(corrected)
            supported_fraction=1.
            for draw in range(DRAWS):
                drawn=[(k,rng.choice(index[k],BSIZE,replace=False))
                       for k in keys["native"]+keys["introduced"]]
                cdist=[];pdens=[];hdens=[]
                for key,pick in drawn:
                    oc=regional_density(xp[pick],edges,h)
                    adjusted,support=correct_occupancy(oc,backgrounds[key])
                    supported_fraction=min(supported_fraction,support)
                    cdist.append(adjusted.ravel())
                    pdens.append(regional_density(xold[pick],photo_edges,h).ravel())
                    hdens.append(regional_density(traitpoints[pick],head_edges[dims],h).ravel())
                n=len(keys["native"])
                corrected+=overlap_matrix(cdist[:n],cdist[n:])/DRAWS
                uncorrected+=overlap_matrix(pdens[:n],pdens[n:])/DRAWS
                heads+=overlap_matrix(hdens[:n],hdens[n:])/DRAWS
            rho,p=infer(corrected,heads)
            rawrho,rawp=infer(uncorrected,heads)
            loo=[]
            for axis in (0,1):
                for j in range(corrected.shape[axis]):
                    a=np.delete(corrected,j,axis).ravel()
                    b=np.delete(heads,j,axis).ravel()
                    loo.append(float(np.corrcoef(rankdata(a),rankdata(b))[0,1]))
            result.append({
                "species":species,"model":model,"head_dimensions":dims,"bandwidth_pc_sd":h,
                "background_corrected_rho":float(rho),
                "background_corrected_permutation_p":float(p),
                "photo_climate_rho_same_cells":float(rawrho),
                "photo_climate_permutation_p":float(rawp),
                "rho_change_due_to_background":float(rho-rawrho),
                "n_native_regions":len(keys["native"]),
                "n_introduced_regions":len(keys["introduced"]),
                "n_region_pairs":int(corrected.size),
                "n_eligible_photos":int(c.photos.sum()),
                "min_occupancy_on_background":supported_fraction,
                "leave_one_region_min":min(loo),"leave_one_region_max":max(loo),
                "background_pc12_variance_fraction":var,
                "within_species_native_introduced_region_matrix_permutations":4999,
            })
    return result


def _novel_head(species,frame,novelty_spec,model):
    """Matched regions and full heads in native-analogous vs native-outside95 climate."""
    transform,edges,native95=novelty_spec
    d=frame.copy()
    coords=transform(d)
    ix=np.column_stack([np.searchsorted(e,coords[:,j],side="right")-1
                        for j,e in enumerate(edges)])
    inside=((ix>=0)&(ix<80)).all(axis=1)
    d=d.loc[inside].copy()
    labels=np.array([not native95[a,b] for a,b in ix[inside]],bool)
    d["novel"]=labels.astype(int)
    nat=d.loc[d.range_state=="native"].copy()
    intro=d.loc[d.range_state=="introduced"].copy()
    # Same original observations and photos in both models; only their
    # environment-defined climate-class labels may change.
    ans=head_novel_analysis(species,nat,intro,.95)
    primary=[v for v in ans if v["design"]=="geographic_pairs" and v["bandwidth"]==.5
             and v["traits"] in ("angle+chroma+outline","angle+chroma")]
    if len(primary)!=2:
        raise ValueError("HEAD_NOVEL_PRIMARY_PAIR_TEST_MISSING")
    rows=[]
    for r in primary:
        rows.append({
            "species":species,"model":model,"head_traits":r["traits"],
            "pair_design":r["design"],"n_native_cells":r["n_native"],
            "n_analogous_cells":r["n_analogous"],"n_novel_cells":r["n_novel"],
            "n_matched_cell_pairs":r["pairs"],"n_geographic_blocks":r["blocks"],
            "D_native_analogous":r["D_native_analogous"],
            "D_native_novel":r["D_native_novel"],
            "D_drop_analogous_minus_novel":r["drop_in_overlap"],
            "regional_label_swap_p":r["p_swap"],
            "same_original_head_photo_cohort":True,
        })
    return rows


def run(wild, land, authority, photo_wild, measured):
    photo,cells=_photo_cells(photo_wild,measured)
    result={
        "status":"SOURCE_MATCHED_V5_S2_ANNUAL_VS_BIO05_VS_BIO10_5AXIS",
        "models":{k:list(v) for k,v in MODELS.items()},
        "same_original_photo_records_in_all_models":True,
        "no_climate_rarity_continuous_score_used":True,
        "results_selected_after_prior_V5_exploratory":True,
        "full_niche":[],"regional_overlap":[],"novel_head_overlap":[]}
    for sp in SPECIES_KEYS:
        spcells=cells.loc[cells.taxon_name==sp].copy().drop(columns=["taxon_name"])
        spframe=photo.loc[photo.taxon_name==sp].copy()
        for model in MODELS:
            niche,info,env=_niche(sp,wild,land,authority,model)
            result["full_niche"].append(niche)
            result["regional_overlap"].extend(_regional_pair(spcells,env,model,sp,info))
            result["novel_head_overlap"].extend(_novel_head(sp,spframe,info,model))
    # Multiplicity is a single 2 species x 2 climate definitions family;
    # different traits / bandwidths are sensitivity, not another first shot.
    reg=result["regional_overlap"]
    primary=[r for r in reg if r["head_dimensions"]==3 and r["bandwidth_pc_sd"]==.5]
    if len(primary)!=6:raise ValueError("S2_BIO10_REGIONAL_PRIMARY_FAMILY_DRIFT")
    vals=multipletests([r["background_corrected_permutation_p"] for r in primary],method="holm")[1]
    for v,r in zip(vals,primary):r["holm_6_spec_species_model"]=float(v)
    novelty=result["novel_head_overlap"]
    p=[r for r in novelty if r["head_traits"]=="angle+chroma+outline"]
    if len(p)!=6:raise ValueError("S2_BIO10_NOVEL_HEAD_PRIMARY_FAMILY_DRIFT")
    vals=multipletests([r["regional_label_swap_p"] for r in p],method="holm")[1]
    for v,r in zip(vals,p):r["holm_6_spec_species_model"]=float(v)
    result["cohort"]={
        "source_wild_n":int(len(wild)),"photo_eligible_n":int(len(photo)),
        "photo_by_species":{sp:int(cells.loc[cells.taxon_name==sp,"photos"].sum()) for sp in SPECIES_KEYS},
        "photo_membership_fixed_before_climate_novelty_classification":True,
        "no_unprotected_individual_ids_public":True,
    }
    return result
