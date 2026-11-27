"""Direct spatial audit of the existing conditional enrichment curves.

Retain all three traits, two climate definitions, fixed cells and 125 bins.
No new significance tests. Reproduce original conditional score tests;
describe regional contributions and an equal-region-weight refit. Count
whether two tail contrasts occur in the same cell pair or regional aggregate.
A cell pair or block is not a biological population or a selection episode.
"""
from pathlib import Path
import json
from itertools import product
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from analysis.v5.phenotype_curve_mechanism import (
    ROOT,SPECIES,DEFINITIONS,VARIABLES,NAMES,EXPECTED_BETA,
    authority_source,angular_cells,background_mesh,apply_original_five_CHELSA,
    load_three_axis,_background_species_audit,_pca_env_2d,run_species,
    X,K,fit,eta_profile,marginals,
)
from analysis.v5.phenotype_within_species_enrichment import sign_patterns,group_sums
OUT=ROOT/"phenotype_curve_regional_structure"

# The original within-species enrichment CSV (artifact 15100013741,
# SHA256 fafbf9f1aa28e54c4e0b803c59d41d93629f7df34f8930253f73fc114189fe7b).
# Order: chroma, angle, outline, then their quadratic terms.
EXPECTED_POS={
 (SPECIES[0],"native_core95"):[6,2,6,2,5,5],
 (SPECIES[0],"density_surplus"):[8,9,6,7,7,6],
 (SPECIES[1],"native_core95"):[10,7,3,3,11,9],
 (SPECIES[1],"density_surplus"):[9,5,5,7,9,7],
}
EXPECTED_P={
 (SPECIES[0],"native_core95"):[.578125,.05859375,.36328125,.125,.62890625,.8046875],
 (SPECIES[0],"density_surplus"):[.767333984375,.693359375,.4088134765625,.9366455078125,.482177734375,.9527587890625],
 (SPECIES[1],"native_core95"):[.018310546875,.43310546875,.4376220703125,.0093994140625,.0101318359375,.896240234375],
 (SPECIES[1],"density_surplus"):[.31890869140625,.2010498046875,.404296875,.72412109375,.2508544921875,.634765625],
}

def contribution(s):
    s=np.asarray(s,float)
    pos=float(np.maximum(s,0).sum())
    neg=float(np.minimum(s,0).sum())
    mass=pos-neg
    return dict(net=float(s.sum()),positive_mass=pos,negative_mass=neg,
                cancellation=float(1-abs(s.sum())/mass) if mass>1e-12 else np.nan,
                largest_absolute_share=float(np.max(np.abs(s))/mass) if mass>1e-12 else np.nan,
                positive_units=int((s>0).sum()),negative_units=int((s<0).sum()),
                zero_units=int((s==0).sum()))

def score_diagnostic(a,b,blocks,sp,definition):
    rows=[]
    flips=sign_patterns(blocks)
    for j in range(6):
        nuisance=np.delete(np.arange(6),j)
        _,ss,hh=fit(a,b,X,nuisance)
        q=ss[:,j]-ss[:,nuisance]@np.linalg.solve(hh[np.ix_(nuisance,nuisance)],hh[nuisance,j])
        qg=group_sums(q[:,None],blocks).ravel()
        stat=flips@qg
        p=float(np.mean(np.abs(stat)>=abs(stat[0])-1e-12))
        assert (qg>0).sum()==EXPECTED_POS[(sp,definition)][j]
        np.testing.assert_allclose(p,EXPECTED_P[(sp,definition)][j],rtol=0,atol=1e-12)
        np.testing.assert_allclose(q.sum(),qg.sum(),rtol=0,atol=1e-10)
        for level,v in (("cell_pair",q),("region",qg)):
            rows.append(dict(trait=NAMES[j%3],term="directional" if j<3 else "quadratic",
                    level=level,units=len(v),original_score_p=p,**contribution(v)))
    return rows

def regional_weights(blocks):
    _,ix,counts=np.unique(blocks,return_inverse=True,return_counts=True)
    w=len(blocks)/(len(counts)*counts[ix])
    np.testing.assert_allclose(w.sum(),len(blocks),rtol=0,atol=1e-10)
    np.testing.assert_allclose(np.bincount(ix,weights=w),len(blocks)/len(counts))
    return w

def weighted_objective(beta,a,b,w):
    # Both groups retain equal total mass per pair; only the pair's
    # contribution to the profiled objective is reweighted.
    z=np.zeros((len(a),5))
    loss,_,p=eta_profile(a,b,beta,z,())
    grad=((a+b)*p-b)@X
    return float(w@loss),w@grad

def weighted_fit(a,b,w,start):
    w=np.asarray(w,float)
    assert np.all(w>0) and w.shape==(len(a),)
    result=minimize(weighted_objective,start,args=(a,b,w),jac=True,method="L-BFGS-B",
                    options=dict(maxiter=300,ftol=1e-14,gtol=1e-8,maxls=40))
    val,grad=weighted_objective(result.x,a,b,w)
    if (not np.isfinite(val) or not np.isfinite(result.x).all()
            or np.max(np.abs(grad))>1e-5 or np.max(np.abs(result.x))>12):
        raise ValueError("unqualified equal-region fit")
    return result.x

def weight_diagnostic(a,b,blocks,beta):
    w=regional_weights(blocks)
    rows=[];curves=[]
    # Both fits use the exact same cells and fixed ranks.
    try:
        weighted=weighted_fit(a,b,w,beta)
    except (ValueError,np.linalg.LinAlgError):
        weighted=None
    for name,bx in (("equal_pair",beta),("equal_region",weighted)):
        for j,trait in enumerate(NAMES):
            for term,k in (("directional",j),("quadratic",j+3)):
                row=dict(weighting=name,trait=trait,term=term,status="estimated" if bx is not None else "fit_failed")
                if bx is not None:
                    row.update(coefficient=float(bx[k]),original_coefficient=float(beta[k]),
                               difference=float(bx[k]-beta[k]))
                rows.append(row)
            for rank in range(-2,3):
                row=dict(weighting=name,trait=trait,native_quintile=rank+3,status="estimated" if bx is not None else "fit_failed")
                if bx is not None:
                    row["enrichment_vs_middle"]=float(np.exp(bx[j]*rank+bx[j+3]*rank*rank))
                curves.append(row)
    return rows,curves

def tail_sign(x):
    # Prevent floating point roundoff from creating spurious sign changes.
    return np.where(x>1e-12,1,np.where(x< -1e-12,-1,0))

def sign_counts(delta,w,level,weighting,trait):
    s=tail_sign(delta[:,[0,4]])
    w=np.asarray(w,float);w=w/w.sum()
    out=[]
    for l,h in product((-1,0,1),repeat=2):
        sel=(s[:,0]==l)&(s[:,1]==h)
        out.append(dict(trait=trait,level=level,weighting=weighting,
            lower_sign=l,upper_sign=h,units=len(s),count=int(sel.sum()),weight_fraction=float(w[sel].sum())))
    assert sum(r["count"] for r in out)==len(s)
    np.testing.assert_allclose(sum(r["weight_fraction"] for r in out),1)
    return out

def tail_diagnostic(a,b,blocks):
    aa,bb=marginals(a),marginals(b)
    np.testing.assert_allclose(aa.sum(2),1,atol=1e-10)
    np.testing.assert_allclose(bb.sum(2),1,atol=1e-10)
    ids=np.unique(blocks)
    sizes=np.array([(blocks==g).sum() for g in ids])
    ga=np.stack([aa[blocks==g].mean(0) for g in ids])
    gb=np.stack([bb[blocks==g].mean(0) for g in ids])
    w=regional_weights(blocks)
    rows=[];profiles=[]
    for j,trait in enumerate(NAMES):
        d=bb[:,j]-aa[:,j];dg=gb[:,j]-ga[:,j]
        np.testing.assert_allclose(d.mean(0),np.average(dg,axis=0,weights=sizes),atol=1e-10)
        for delta,weights,level,weighting in (
          (d,np.ones(len(d)),"cell_pair","equal_pair"),
          (d,w,"cell_pair","equal_region"),
          (dg,sizes,"region","equal_pair"),
          (dg,np.ones(len(ids)),"region","equal_region"),
        ):
            rows.extend(sign_counts(delta,weights,level,weighting,trait))
        for name,weights in (("equal_pair",np.ones(len(a))),("equal_region",w)):
            pa=np.average(aa[:,j],axis=0,weights=weights)
            pb=np.average(bb[:,j],axis=0,weights=weights)
            for k in range(5):
                profiles.append(dict(trait=trait,weighting=name,native_quintile=k+1,
                     analogous=float(pa[k]),novel=float(pb[k]),difference=float(pb[k]-pa[k])))
    return rows,profiles

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=load_three_axis(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    scores=[];weights=[];curves=[];tails=[];marginal=[]
    for definition in DEFINITIONS:
        for sp in SPECIES:
            d=out[(sp,definition)]
            a,b,g=d["hist_left"],d["hist_right"],d["blocks"]
            meta=dict(species=sp,definition=definition,pairs=len(a),blocks=len(np.unique(g)))
            beta,_,_=fit(a,b,X,np.arange(6))
            np.testing.assert_allclose(beta,EXPECTED_BETA[(sp,definition)],atol=1e-6,rtol=0)
            scores.extend(dict(**meta,**r) for r in score_diagnostic(a,b,g,sp,definition))
            wrows,crows=weight_diagnostic(a,b,g,beta)
            weights.extend(dict(**meta,**r) for r in wrows)
            curves.extend(dict(**meta,**r) for r in crows)
            trows,mrows=tail_diagnostic(a,b,g)
            tails.extend(dict(**meta,**r) for r in trows)
            marginal.extend(dict(**meta,**r) for r in mrows)
    pd.DataFrame(scores).to_csv(OUT/"adjusted_score_contributions.csv",index=False)
    pd.DataFrame(weights).to_csv(OUT/"regional_weight_coefficients.csv",index=False)
    pd.DataFrame(curves).to_csv(OUT/"regional_weight_curves.csv",index=False)
    pd.DataFrame(tails).to_csv(OUT/"paired_tail_signs.csv",index=False)
    pd.DataFrame(marginal).to_csv(OUT/"marginal_profiles.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
