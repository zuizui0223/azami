"""Exploratory observed environmental composition and within-climate phenotype distributions.

Freeze the continuous B. Same original eligible cells and photographic
records. Two diagnostics, no selection/causal inference, no new P search:
(1) climate-caliper rematching within existing regional blocks, partition
    marginal density differences into matched and unmatched components;
(2) continuous local climate profiles with equal regional mass in the two
    groups, retaining the entire phenotype distribution (not quintiles).
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from analysis.v5.phenotype_continuous_enrichment import (
    ROOT,SPECIES,DEFINITIONS,VARIABLES,RAIN,SEASON,NAMES,COLS,EXPECTED,
    EXPECTED_BETA,joint_records,eligible,load_three_axis,make_cells,match_cells,
    head_cells,background_mesh,_background_species_audit,_pca_env_2d,
    _probability_content,apply_original_five_CHELSA,authority_source,
    angular_cells,X,K,native_cuts,rank_codes,cell_hist,fit,climate_cells,
    native_reference,raw_pairs,fit_raw,kde_mass,curve_grid,BANDWIDTH,
)
OUT=ROOT/"phenotype_environment_explanation"
# Same per-axis RMS standardized analogue tolerance as Part 1
# (0.9 Euclidean distance for 3 climatic dimensions). Never tuned to outcomes.
CALIPER=.9/np.sqrt(3.)
MAX_AXIS=1.
MIN_PAIRS=25
MIN_REGIONS=8
ENV_BW=.25
MIN_LOCAL_ESS=20.
MIN_LOCAL_REGIONS=3.
GRID=41
# Existing continuous analysis: order must match NAMES, linear then quadratic.
CONT={
 (SPECIES[0],"native_core95"):[.305947,-.330475,.656272,-.174406,.010560,-.311833],
 (SPECIES[1],"native_core95"):[.456362,.003420,-.072591,-.165495,.175251,.017357],
 (SPECIES[0],"density_surplus"):[.053459,-.046998,-.087837,-.033377,.026150,.026261],
 (SPECIES[1],"density_surplus"):[.120866,-.139475,.019680,.039265,.057016,-.018520],
}
SCHEMES={axis:(j,) for j,axis in enumerate(VARIABLES)}
SCHEMES.update(water_joint=(RAIN,SEASON),all_five=tuple(range(len(VARIABLES))))

def match_one(a,b,axes):
    """Maximum cardinality then minimum RMS climate distance.

    Dummies allow every left cell to be unmatched. A dummy's cost is
    greater than the sum of any admissible real distances, ensuring that
    one extra matched pair always dominates a cost improvement.
    """
    a,b=np.asarray(a,float),np.asarray(b,float)
    n,m=len(a),len(b)
    if not n or not m:
        return np.empty(0,int),np.empty(0,int)
    delta=a[:,None,axes]-b[None,:,axes]
    rms=np.sqrt(np.mean(delta**2,axis=2))
    valid=(rms<=CALIPER+1e-12)&(np.max(np.abs(delta),axis=2)<=MAX_AXIS+1e-12)
    dummy=(n+1)*(1+CALIPER)
    cost=np.full((n,m+n),dummy)
    cost[:,:m]=np.where(valid,rms,dummy*(n+2))
    l,r=linear_sum_assignment(cost)
    good=r<m;l,r=l[good],r[good]
    assert valid[l,r].all()
    return l,r

def rematch(left,right,blocks,z,axes):
    la=[];ra=[];ga=[]
    for g in np.unique(blocks):
        ids=np.flatnonzero(blocks==g)
        l,r=match_one(z[left[ids]],z[right[ids]],axes)
        la.extend(left[ids[l]]);ra.extend(right[ids[r]]);ga.extend([g]*len(l))
    a,b=np.asarray(la,int),np.asarray(ra,int)
    assert len(np.unique(a))==len(a)==len(b)==len(np.unique(b))
    return a,b,np.asarray(ga,dtype=blocks.dtype)

def group_balance(x,z,a,b,axes,meta):
    if not len(a):
        return [dict(**meta,axis=axis,matched_axis=j in axes,status="empty")
                for j,axis in enumerate(VARIABLES)]
    rows=[]
    for j,axis in enumerate(VARIABLES):
        delta=z[b,j]-z[a,j]
        rows.append(dict(**meta,axis=axis,matched_axis=j in axes,status="descriptive",
                  analogous_mean=float(x[a,j].mean()),novel_mean=float(x[b,j].mean()),
                  signed_difference_native_sd=float(delta.mean()),
                  mean_abs_paired_difference_native_sd=float(abs(delta).mean()),
                  p90_abs_paired_difference_native_sd=float(np.quantile(abs(delta),.9)),
                  max_abs_paired_difference_native_sd=float(abs(delta).max())))
    return rows

def partition_density(a,b,selected_a,selected_b):
    """Exact finite-sample additive partition of the marginal density gap.

    a,b: same numbers of cell-weighted phenotype densities; selected indices
    give an equal number of cells from each group. No imputed outcome.
    """
    n=len(a);m=len(selected_a)
    assert n==len(b) and m==len(selected_b)
    full=b.mean(0)-a.mean(0)
    if m:
        shared=b[selected_b].mean(0)-a[selected_a].mean(0)
        contrib=m/n*shared
    else:
        shared=np.full_like(full,np.nan);contrib=np.zeros_like(full)
    ua=np.setdiff1d(np.arange(n),selected_a)
    ub=np.setdiff1d(np.arange(n),selected_b)
    assert len(ua)==len(ub)==n-m
    non=(b[ub].sum(0)-a[ua].sum(0))/n
    np.testing.assert_allclose(full,contrib+non,atol=1e-10,rtol=0)
    return full,shared,contrib,non

def local_weights(xa,xb,blocks,at):
    """Continuous compactly supported climate kernels, common regional mass.

    Local region weights are proportional to the minimum kernel mass
    in the two groups; no borrowing across regions or extrapolation.
    xa and xb are on the same native-SD climate scale.
    """
    u=(xa-at)/ENV_BW;v=(xb-at)/ENV_BW
    ka=np.exp(-.5*u*u)*(abs(u)<=2.)
    kb=np.exp(-.5*v*v)*(abs(v)<=2.)
    wa=np.zeros(len(xa));wb=np.zeros(len(xb))
    masses=[]
    for g in np.unique(blocks):
        use=blocks==g;sa=ka[use].sum();sb=kb[use].sum()
        common=min(sa,sb)
        if common>0:
            wa[use]=ka[use]*common/sa
            wb[use]=kb[use]*common/sb
            masses.append(common)
    total=wa.sum()
    if total<=0:
        return wa,wb,0.,0.,0.,0
    np.testing.assert_allclose(total,wb.sum(),atol=1e-10,rtol=0)
    for g in np.unique(blocks):
        np.testing.assert_allclose(wa[blocks==g].sum(),wb[blocks==g].sum(),atol=1e-10,rtol=0)
    wa/=total;wb/=total
    mass=np.asarray(masses)/total
    return wa,wb,float(1/(wa@wa)),float(1/(wb@wb)),float(1/(mass@mass)),len(masses)

def weighted_moments(v,ix,valid,kept,cells):
    vv=v[valid];ii=ix[valid]
    w=1/kept[ii]
    first=np.zeros((len(kept),3));second=np.zeros_like(first)
    np.add.at(first,ii,w[:,None]*vv)
    np.add.at(second,ii,w[:,None]*vv**2)
    return first[cells],second[cells]

def coef_rows(out,meta,beta,center,scale,qq,status):
    for j,trait in enumerate(NAMES):
        r=dict(**meta,trait=trait,status=status,native_center=float(center[j]),native_sd=float(scale[j]))
        if beta is not None:
            b,c=beta[j],beta[j+3]
            r.update(linear=float(b),quadratic=float(c),
               turning_point=float(center[j]-b*scale[j]/(2*c)) if abs(c)>1e-9 else np.nan,
               turning_point_inside_native_q05_q95=bool(abs(c)>1e-9 and
                     qq[j,0]<=center[j]-b*scale[j]/(2*c)<=qq[j,-1]))
            for value in curve_grid(qq,j):
                z=(value-center[j])/scale[j]
                out["curves"].append(dict(**meta,trait=trait,value=float(value),native_z=float(z),
                     log_relative_enrichment=float(b*z+c*z*z),
                     relative_enrichment=float(np.exp(b*z+c*z*z))))
        out["coefficients"].append(r)

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,good=rank_codes(nv,cuts,lo,hi);assert good.all()
    support=np.bincount(rc,minlength=K)>0;assert support.sum()==125
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    code,valid=rank_codes(iv,cuts,lo,hi);valid&=support[code]
    np.testing.assert_array_equal(np.bincount(ii[valid],minlength=len(c)),kept)
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    blocks=c.iloc[left].block.to_numpy()
    np.testing.assert_array_equal(blocks,c.iloc[right].block.to_numpy())
    old=fit(h[left],h[right],X,np.arange(6))[0]
    np.testing.assert_allclose(old,EXPECTED_BETA[(sp,definition)],atol=1e-6,rtol=0)
    center,scale,qq=native_reference(nv,ni,nw,len(n))
    original=raw_pairs(iv,ii,valid,left,right,kept,center,scale,blocks)
    beta=fit_raw(original,np.arange(6))[0]
    # Frozen continuous fit (published to six decimal places).
    np.testing.assert_allclose(beta,CONT[(sp,definition)],atol=2e-6,rtol=0)
    x=climate_cells(fi,c);xn=climate_cells(fn,n)
    mu=xn.mean(0);sd=xn.std(0)
    assert np.isfinite(sd).all() and (sd>0).all()
    z=(x-mu)/sd
    base=dict(species=sp,definition=definition,original_pairs=len(left),original_regions=len(np.unique(blocks)))
    # Save both curves and unconditional smooth densities. Full-record
    # density units are inverse physical trait units.
    densities={}
    for j in range(3):
        grid=curve_grid(qq,j);bw=BANDWIDTH*scale[j]
        norm=1/(np.sqrt(2*np.pi)*bw)
        aa=kde_mass(iv[valid,j],ii[valid],kept,left,grid,bw)*norm
        bb=kde_mass(iv[valid,j],ii[valid],kept,right,grid,bw)*norm
        densities[j]=(grid,aa,bb)
    left_index={v:i for i,v in enumerate(left)}
    right_index={v:i for i,v in enumerate(right)}
    comparisons=[("original",(),left,right,blocks)]
    for scheme,axes in SCHEMES.items():
        a,b,g=rematch(left,right,blocks,z,axes)
        comparisons.append((scheme,axes,a,b,g))
    for scheme,axes,a,b,g in comparisons:
        m=len(a);ng=len(np.unique(g))
        status="eligible" if m>=MIN_PAIRS and ng>=MIN_REGIONS else "limited_support"
        meta=dict(**base,scheme=scheme,axes=";".join(VARIABLES[j] for j in axes),
                  pairs=m,regions=ng,retained_fraction=m/len(left))
        out["cohort"].append(dict(**meta,status=status,
             analogous_records=int(kept[a].sum()),novel_records=int(kept[b].sum()),
             mean_rms_distance=float(np.sqrt(np.mean((z[b][:,axes]-z[a][:,axes])**2,axis=1)).mean())
                               if m and len(axes) else np.nan))
        out["balance"].extend(group_balance(x,z,a,b,axes,meta))
        if scheme=="original":
            coef_rows(out,meta,beta,center,scale,qq,"original_reproduced")
            continue
        if status=="eligible":
            try:
                d=raw_pairs(iv,ii,valid,a,b,kept,center,scale,g)
                fit_beta=fit_raw(d,np.arange(6))[0]
                coef_rows(out,meta,fit_beta,center,scale,qq,"descriptive_refit")
            except (ValueError,np.linalg.LinAlgError,AssertionError):
                coef_rows(out,meta,None,center,scale,qq,"fit_failed")
        else:
            coef_rows(out,meta,None,center,scale,qq,"limited_support")
        ia=np.array([left_index[v] for v in a],int)
        ib=np.array([right_index[v] for v in b],int)
        for j,(grid,aa,bb) in densities.items():
            full,shared,contrib,non=partition_density(aa,bb,ia,ib)
            for k,value in enumerate(grid):
                ma=float(aa[ia,k].mean()) if m else np.nan
                mb=float(bb[ib,k].mean()) if m else np.nan
                out["density_partition"].append(dict(**meta,trait=NAMES[j],value=float(value),
                     original_analogous_density=float(aa[:,k].mean()),
                     original_novel_density=float(bb[:,k].mean()),
                     matched_analogous_density=ma,matched_novel_density=mb,
                     original_difference=float(full[k]),matched_difference=float(shared[k]),
                     matched_contribution=float(contrib[k]),unmatched_contribution=float(non[k]),
                     status=status))
    # Direct P(phenotype | measured climate) with continuous local kernels,
    # equalizing the regional mixture at every evaluated climate value.
    # No new regression or climate-to-fitness conversion.
    f_a,q_a=weighted_moments(iv,ii,valid,kept,left)
    f_b,q_b=weighted_moments(iv,ii,valid,kept,right)
    for axis in range(len(VARIABLES)):
        xa,xb=z[left,axis],z[right,axis]
        low,high=np.quantile(np.r_[xa,xb],[.05,.95])
        for at in np.linspace(low,high,GRID):
            wa,wb,esa,esb,ges,nreg=local_weights(xa,xb,blocks,at)
            meta=dict(**base,axis=VARIABLES[axis],climate_value=float(mu[axis]+sd[axis]*at),
                  climate_z=float(at),climate_bandwidth=float(ENV_BW*sd[axis]),
                  effective_analogous_cells=esa,effective_novel_cells=esb,
                  effective_regions=ges,regions_with_overlap=nreg,
                  status="descriptive" if min(esa,esb)>=MIN_LOCAL_ESS and ges>=MIN_LOCAL_REGIONS else "sparse")
            if not nreg:
                out["climate_profiles"].append(dict(**meta,trait="all",status_empty=True))
                continue
            for j,trait in enumerate(NAMES):
                va=float(wa@f_a[:,j]);vb=float(wb@f_b[:,j])
                sa=float(np.sqrt(max(wa@q_a[:,j]-va*va,0)))
                sb=float(np.sqrt(max(wb@q_b[:,j]-vb*vb,0)))
                out["climate_profiles"].append(dict(**meta,trait=trait,
                    analogous_mean=va,novel_mean=vb,mean_difference=vb-va,
                    analogous_sd=sa,novel_sd=sb,sd_difference=sb-sa))
                # For plotting, use every second point of the existing
                # 101-point trait evaluation grid; input traits remain raw.
                grid,aa,bb=densities[j]
                da=wa@aa;db=wb@bb
                for k in range(0,len(grid),2):
                    out["conditional_density"].append(dict(**meta,trait=trait,
                         value=float(grid[k]),analogous_density=float(da[k]),
                         novel_density=float(db[k]),density_difference=float(db[k]-da[k])))

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
    out={k:[] for k in ("cohort","balance","coefficients","curves","density_partition","climate_profiles","conditional_density")}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    for key,rows in out.items():
        pd.DataFrame(rows).to_csv(OUT/(key+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
