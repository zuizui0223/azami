"""Fixed matched-cell water-response test for the two focal species.
Species-specific climate novelty, species-relative native chroma tails;
shared physical climate scales. Cross-sectional phenotype association, not
selection, causation, adaptation rate or a common-garden treatment.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
from scipy.linalg import null_space
from scipy.spatial import Delaunay, QhullError
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_species_contrast import (
    ROOT,NAMES,COLS,EXPECTED,EXPECTED_BETA,joint_records,T,eligible,
    load_three_axis,make_cells,match_cells,head_cells,background_mesh,
    _background_species_audit,_pca_env_2d,_probability_content,
    apply_original_five_CHELSA,authority_source,angular_cells,VARIABLES,
    X,K,native_cuts,rank_codes,cell_hist,fit,
)
from analysis.v5.phenotype_redistribution import cellkeys
OUT=ROOT/"phenotype_water_interaction"
SPECIES=("Cirsium arvense","Cirsium vulgare")
DEFINITIONS=("native_core95","density_surplus")
ENDPOINTS=("existing_high20","existing_low20")
RAIN=list(VARIABLES).index("chelsa_bio12")
SEASON=list(VARIABLES).index("chelsa_bio15")
WATER_TERMS=(("rain",RAIN+1),("seasonality",SEASON+1),("rain_x_seasonality",6))
P=7
SEED=20261127

def climate_cells(frame,cells):
    # The earlier climate-context model used all complete raw three-trait
    # records to form cell climate medians. Preserve that exposure definition.
    f=frame.groupby("cellkey")[list(VARIABLES)].median()
    x=f.reindex(cellkeys(cells)).to_numpy(float)
    assert np.isfinite(x).all()
    return x

def existing_tail_matrix(v,ix,n,cuts,lo,hi):
    j=NAMES.index("chroma")
    z=v[:,j]
    flags=np.column_stack(((z>=cuts[j,-1])&(z<=hi[j]),
                           (z>=lo[j])&(z<cuts[j,0])))
    counts=np.bincount(ix,minlength=n)
    out=np.stack([np.bincount(ix,weights=flags[:,k],minlength=n)
                  for k in range(2)],axis=1)
    return np.divide(out,counts[:,None],out=np.full_like(out,np.nan),
                     where=counts[:,None]>0)

def common_scale(a,b):
    # Equal species contributions; pooled native cell climate means and
    # variance. For water use fixed physical steps, 500 mm and 20 BIO15 units.
    center=(a.mean(0)+b.mean(0))/2
    scale=np.sqrt(((a-center)**2).mean(0)/2+((b-center)**2).mean(0)/2)
    scale[RAIN]=500.
    scale[SEASON]=20.
    assert np.isfinite(center).all() and (scale>0).all()
    return center,scale

def difference_features(left,right,center,scale):
    a=(left-center)/scale
    b=(right-center)/scale
    delta=b-a
    # The difference in level-products is the correct pair difference of
    # the nonlinear surface. Do NOT use delta_rain * delta_seasonality.
    cross=b[:,RAIN]*b[:,SEASON]-a[:,RAIN]*a[:,SEASON]
    return np.column_stack((np.ones(len(a)),delta,cross))

def prepare(d,center,scale):
    return dict(**d,design=difference_features(d["left"],d["right"],center,scale))

def summarize_exposure(sp,definition,d):
    rows=[]
    for j,name in enumerate(VARIABLES):
        a,b=d["left"][:,j],d["right"][:,j]
        rows.append(dict(species=sp,definition=definition,variable=name,
            pairs=len(a),blocks=len(np.unique(d["blocks"])),
            analogous_median=float(np.median(a)),novel_median=float(np.median(b)),
            mean_delta=float(np.mean(b-a)),median_delta=float(np.median(b-a)),
            delta_q05=float(np.quantile(b-a,.05)),delta_q95=float(np.quantile(b-a,.95))))
    for j,name in enumerate(ENDPOINTS):
        a,b=d["tail_left"][:,j],d["tail_right"][:,j]
        rows.append(dict(species=sp,definition=definition,variable=name,
            pairs=len(a),blocks=len(np.unique(d["blocks"])),
            analogous_median=float(np.median(a)),novel_median=float(np.median(b)),
            analogous_mean=float(a.mean()),novel_mean=float(b.mean()),
            mean_delta=float(np.mean(b-a)),median_delta=float(np.median(b-a)),
            delta_q05=float(np.quantile(b-a,.05)),delta_q95=float(np.quantile(b-a,.95))))
    return rows

def in_hull(ref,target):
    # Diagnose joint support without fitting an outcome or extrapolating.
    # Standardize by reference variation for numerical conditioning only.
    mu=ref.mean(0);sd=ref.std(0)
    if np.any(sd<1e-10) or len(ref)<=ref.shape[1]:
        return None
    a=(ref-mu)/sd;b=(target-mu)/sd
    if np.linalg.matrix_rank(a)<a.shape[1]:
        return None
    try:
        return Delaunay(np.unique(a,axis=0)).find_simplex(b,tol=1e-9)>=0
    except QhullError:
        return None

def overlap_diagnostic(definition,prepared):
    rows=[]
    for i,d in enumerate(prepared):
        other=prepared[1-i]
        # Joint support of water contrasts, including their level interaction.
        a=d["design"][:,[RAIN+1,SEASON+1,6]]
        b=other["design"][:,[RAIN+1,SEASON+1,6]]
        good=in_hull(b,a)
        # Also require both endpoint water levels to be within the other
        # species' corresponding endpoint hull. A descriptive diagnostic,
        # not a claim of full five-dimensional or geographic matching.
        gl=in_hull(other["left"][:,[RAIN,SEASON]],d["left"][:,[RAIN,SEASON]])
        gr=in_hull(other["right"][:,[RAIN,SEASON]],d["right"][:,[RAIN,SEASON]])
        for kind,g in (("water_difference_hull",good),("both_endpoint_hulls",None if gl is None or gr is None else gl&gr),
                       ("intersection",None if good is None or gl is None or gr is None else good&gl&gr)):
            row=dict(species=SPECIES[i],definition=definition,diagnostic=kind,
                     all_pairs=len(a),all_blocks=len(np.unique(d["blocks"])),
                     status="estimated" if g is not None else "degenerate_hull")
            if g is not None:
                row.update(supported_pairs=int(g.sum()),supported_fraction=float(g.mean()),
                           supported_blocks=len(np.unique(d["blocks"][g])))
            rows.append(row)
    return rows

def supported_design(d):
    z=d["design"]
    return (len(z)>=max(25,3*P) and len(np.unique(d["blocks"]))>=8
            and np.isfinite(z).all() and np.linalg.matrix_rank(z)==P
            and np.linalg.cond(z)<1e6)

def contrast_vectors():
    rows=[]
    for name,j in WATER_TERMS:
        a=np.zeros(2*P);b=a.copy();a[j]=1;b[P+j]=1
        rows.extend([(name,SPECIES[0],a),(name,SPECIES[1],b),
                     (name,"C_vulgare_minus_C_arvense",b-a)])
    return rows

def stack_model(prepared,endpoint):
    k=ENDPOINTS.index(endpoint)
    zz=[];yy=[];bb=[]
    for i,d in enumerate(prepared):
        n=len(d["design"])
        z=np.zeros((n,2*P));z[:,i*P:(i+1)*P]=d["design"]
        zz.append(z);yy.append(d["tail_right"][:,k]-d["tail_left"][:,k]);bb.append(d["blocks"])
    return np.vstack(zz),np.concatenate(yy),np.concatenate(bb)

def score_sign_test(z,y,blocks,c):
    # Null-space reparameterization: impose c' beta = 0, allow ALL nuisance
    # coefficients to vary. Residualize the one contrast direction. Joint
    # block flips preserve shared geography across species. With fitted
    # nuisance this is an approximate score-symmetry test, not randomization.
    nuisance=z@null_space(c[None,:])
    x=z@c/(c@c)
    xr=x-nuisance@np.linalg.lstsq(nuisance,x,rcond=None)[0]
    yr=y-nuisance@np.linalg.lstsq(nuisance,y,rcond=None)[0]
    ids=np.unique(blocks)
    s=np.array([xr[blocks==g]@yr[blocks==g] for g in ids])
    if len(ids)<=15:
        flips=np.array(list(product([-1.,1.],repeat=len(ids))))
        add=0
    else:
        flips=np.random.default_rng(SEED).choice([-1.,1.],size=(19999,len(ids)))
        add=1
    v=flips@s
    p=(np.count_nonzero(np.abs(v)>=abs(s.sum())-1e-12)+add)/(len(v)+add)
    return float(p),len(v),int((s>0).sum())

def fit_water(prepared,endpoint):
    meta=dict(ca_pairs=len(prepared[0]["design"]),cv_pairs=len(prepared[1]["design"]),
              ca_blocks=len(np.unique(prepared[0]["blocks"])),cv_blocks=len(np.unique(prepared[1]["blocks"])))
    ok=all(supported_design(d) for d in prepared)
    if not ok:
        return [dict(**meta,term=name,comparison=sp,status="insufficient_support",p_t=1.,p_score=1.)
                for name,sp,c in contrast_vectors()]
    z,y,blocks=stack_model(prepared,endpoint)
    beta=np.linalg.lstsq(z,y,rcond=None)[0]
    resid=y-z@beta
    ids=np.unique(blocks)
    # Species-specific small-sample correction on influence contributions.
    # Combine aligned regional influences before forming the covariance,
    # so shared-region cross-species dependence is not discarded.
    influences=np.zeros((len(ids),2*P))
    for i,d in enumerate(prepared):
        ii=slice(i*P,(i+1)*P)
        x=d["design"]
        k=ENDPOINTS.index(endpoint)
        r=d["tail_right"][:,k]-d["tail_left"][:,k]-x@beta[ii]
        g=d["blocks"]
        b=len(np.unique(g));n=len(g)
        bread=np.linalg.inv(x.T@x)
        correction=np.sqrt(b/(b-1)*(n-1)/(n-P))
        for j,group in enumerate(ids):
            influences[j,ii]=(x[g==group].T@r[g==group])@bread*correction
    cov=influences.T@influences
    df=min(meta["ca_blocks"],meta["cv_blocks"])-1
    crit=student_t.ppf(.975,df)
    rows=[]
    for name,sp,c in contrast_vectors():
        est=float(c@beta)
        se=float(np.sqrt(max(c@cov@c,0)))
        if se<1e-12:
            rows.append(dict(**meta,term=name,comparison=sp,status="degenerate_se",p_t=1.,p_score=1.))
            continue
        p,patterns,npos=score_sign_test(z,y,blocks,c)
        rows.append(dict(**meta,term=name,comparison=sp,status="estimated",
            coefficient=est,cluster_se=se,lo=est-crit*se,hi=est+crit*se,df=df,
            p_t=float(2*student_t.sf(abs(est/se),df)),p_score=p,
            union_blocks=len(ids),score_patterns=patterns,positive_score_blocks=npos))
    return rows

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,valid=rank_codes(nv,cuts,lo,hi)
    assert valid.all()
    support=np.bincount(rc,minlength=K)>0
    assert support.sum()==125
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    blocks=c.iloc[left].block.to_numpy()
    np.testing.assert_array_equal(blocks,c.iloc[right].block.to_numpy())
    assert np.all(kept[np.r_[left,right]]>0)
    beta,_,_=fit(h[left],h[right],X,np.arange(6))
    np.testing.assert_allclose(beta,EXPECTED_BETA[(sp,definition)],atol=1e-6,rtol=0)
    tails=existing_tail_matrix(iv,ii,len(c),cuts,lo,hi)
    climate=climate_cells(fi,c)
    assert np.isfinite(tails[np.r_[left,right]]).all()
    out[(sp,definition)]=dict(native_climate=climate_cells(fn,n),
        left=climate[left],right=climate[right],blocks=blocks,
        tail_left=tails[left],tail_right=tails[right])

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    ex=EXPECTED[sp]
    assert (len(n),len(full))==ex[:2]
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    labeled["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[2]
    analyze(sp,"native_core95",n,c,pairs,fn,fi,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[3]
    analyze(sp,"density_surplus",n,c,pairs,fn,fi,out)

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
    coef=[];profiles=[];support=[];scale_rows=[]
    for definition in DEFINITIONS:
        ca=out[(SPECIES[0],definition)]
        cv=out[(SPECIES[1],definition)]
        center,scale=common_scale(ca["native_climate"],cv["native_climate"])
        for j,axis in enumerate(VARIABLES):
            scale_rows.append(dict(definition=definition,axis=axis,center=center[j],scale=scale[j]))
        prepared=[prepare(d,center,scale) for d in (ca,cv)]
        for sp,d in zip(SPECIES,prepared):
            profiles.extend(summarize_exposure(sp,definition,d))
        support.extend(overlap_diagnostic(definition,prepared))
        for endpoint in ENDPOINTS:
            rows=fit_water(prepared,endpoint)
            coef.extend(dict(definition=definition,endpoint=endpoint,**r) for r in rows)
    d=pd.DataFrame(coef)
    assert len(d)==36
    for p in ("p_t","p_score"):
        d["holm_all36_"+p]=multipletests(d[p].fillna(1),method="holm")[1]
        d["holm_definition18_"+p]=d.groupby("definition")[p].transform(lambda p:multipletests(p.fillna(1),method="holm")[1])
    d.to_csv(OUT/"water_interactions.csv",index=False)
    pd.DataFrame(profiles).to_csv(OUT/"paired_profiles.csv",index=False)
    pd.DataFrame(support).to_csv(OUT/"water_common_support.csv",index=False)
    pd.DataFrame(scale_rows).to_csv(OUT/"common_climate_scales.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
