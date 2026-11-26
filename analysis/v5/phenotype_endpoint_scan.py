"""Exploratory endpoint-wide native-reference distribution scan.

Keep the two V5 climate assignments and original 3-trait complete-case geographic
pairs fixed before inspecting additional archived image endpoints. A source
column/usable-head gate and pair attrition report accompany every endpoint.
No individual observations, locations or block IDs are exported.
"""
import json
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_redistribution import ROOT, SEED, eligible, cellkeys
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.invasion_partition import head_cells
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d, _probability_content
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.s1_s2_extreme_lego import S2_RENAMING
from analysis.v4.native_integration import FEATURES

OUT = ROOT/"phenotype_endpoint_scan"
CORE = dict(zip(["orientation_image_vertical_angle","corolla_lab_chroma",
                 "capitulum_outline_aspect_ratio"],T))
PROBS = np.arange(1,20)/20
MODES = ["record_equalcell","cell_summary"]
DEFS = ["native_core95","density_surplus"]

def load_measures(path):
    header=set(pd.read_csv(path,nrows=0).columns)
    if not set(["obs_id",*T])<=header: raise ValueError("Original V5 trait columns missing")
    mapping={}
    for endpoint,components in FEATURES.items():
        for key in components:
            candidates=list(dict.fromkeys([CORE.get(key,key),key,*S2_RENAMING[key]]))
            col=next((c for c in candidates if c in header),None)
            flag=None
            if col:
                options=[col.removesuffix("_median")+"_n_usable_heads",key+"_n_usable_heads"]
                flag=next((c for c in options if c in header),None)
            mapping[key]=(col,flag)
    use=list(dict.fromkeys(["obs_id",*T,*[c for p in mapping.values() for c in p if c]]))
    raw=pd.read_csv(path,usecols=use,dtype={"obs_id":str},low_memory=False)
    if raw.obs_id.duplicated().any():raise ValueError("Original measurement ID duplicated")
    m=raw[["obs_id",*T]].copy()
    for key,(col,flag) in mapping.items():
        x=pd.to_numeric(raw[col],errors="coerce") if col else pd.Series(np.nan,index=raw.index)
        x=x.where(np.isfinite(x))
        if flag:x=x.where(pd.to_numeric(raw[flag],errors="coerce")>0)
        if key=="orientation_image_vertical_angle":x=x.where(x.between(0,180))
        m["scan_"+key]=x
    return m,mapping

def cell_average_matrix(v,ids,n):
    v=np.asarray(v,float)
    if v.ndim==1:v=v[:,None]
    out=np.zeros((n,v.shape[1]))
    np.add.at(out,ids,v)
    count=np.bincount(ids,minlength=n)
    result=np.full_like(out,np.nan)
    np.divide(out,count[:,None],where=count[:,None]>0,out=result)
    return result

def values(frame,cells,columns,mode,circular):
    idx=pd.Index(cellkeys(cells)).get_indexer(frame.cellkey)
    x=frame[columns].to_numpy(float)
    ok=(idx>=0)&np.isfinite(x).all(1)
    if circular:ok &= np.linalg.norm(x,axis=1)>1e-8
    ix,x=idx[ok],x[ok]
    if circular:
        # Archived sine and cosine are a pair; no linear averaging across the hue seam.
        theta=np.arctan2(x[:,0],x[:,1])
        if mode=="cell_summary":
            m=cell_average_matrix(np.column_stack([np.sin(theta),np.cos(theta)]),ix,len(cells))
            good=np.isfinite(m).all(1)&(np.linalg.norm(m,axis=1)>1e-8)
            theta=np.arctan2(m[good,0],m[good,1]);ix=np.flatnonzero(good)
        x=theta
    else:
        x=x[:,0]
        if mode=="cell_summary":
            g=pd.Series(x).groupby(ix).median()
            x=g.to_numpy();ix=g.index.to_numpy(int)
    counts=np.bincount(ix,minlength=len(cells))
    return x,ix,1/counts[ix],int(ok.sum())

def empirical_rank(x,ref,w):
    """Weighted midrank against frozen native observations, including ties."""
    idx=np.argsort(ref,kind="stable")
    ref=ref[idx];w=w[idx]
    cm=np.r_[0,np.cumsum(w)/w.sum()]
    a=np.searchsorted(ref,x,side="left")
    b=np.searchsorted(ref,x,side="right")
    return (cm[a]+cm[b])/2

def design_matrices(x,ix,n,ref,rw,circular):
    if circular:
        # Circular Fourier discrepancy: weighted first three harmonics.
        cols=np.column_stack([f(k*x)/k for k in [1,2,3] for f in [np.sin,np.cos]])
        return cell_average_matrix(cols,ix,n),cell_average_matrix(np.column_stack([np.sin(x),np.cos(x)]),ix,n),None
    u=empirical_rank(x,ref,rw)
    feat=cell_average_matrix(u[:,None]<=PROBS[None,:],ix,n)
    lo,hi=ref.min(),ref.max()
    summary=cell_average_matrix(np.column_stack([x,u,(u<=.2)&(x>=lo),(u>=.8)&(x<=hi),(x<lo)|(x>hi)]),ix,n)
    return feat,summary,u

def permute(blocks,counts):
    # All pairs in a geographic block change label together; native reference fixed.
    k=len(counts);exact=k<=15
    rng=np.random.default_rng(SEED)
    flips=np.array(list(product([-1.,1.],repeat=k))) if exact else rng.choice([-1.,1.],(9999,k))
    observed=blocks.sum(0)/counts.sum()
    null=flips@blocks/counts.sum()
    score=float(np.sqrt(np.mean(observed**2)))
    stat=np.sqrt(np.mean(null**2,axis=1))
    p=(np.count_nonzero(stat>=score-1e-12)+int(not exact))/(len(flips)+int(not exact))
    return score,p,exact,len(flips)

def scalar_summary(x):
    names=["mean_raw","mean_native_rank","existing_low20","existing_high20","outside_native_observed_range"]
    return dict(zip(names,x))

def circular_summary(x):
    strength=np.linalg.norm(x)
    return dict(mean_hue_degrees=(np.degrees(np.arctan2(x[0],x[1]))%360 if strength>1e-8 else np.nan),
                circular_resultant=strength)

def run_scan(sp,definition,c,pairs,fn,fi,n,mapping,output):
    total=sum(len(a) for a,b in pairs)
    for trait,components in FEATURES.items():
        circular=(trait=="hue")
        cols=["scan_"+key for key in components]
        for mode in MODES:
            row=dict(species=sp,definition=definition,representation=mode,trait=trait,
                     fixed_pairs=total,fixed_blocks=len(pairs),
                     source_columns=";".join(str(mapping[key][0] or "unavailable") for key in components),
                     usable_head_gate=";".join(str(mapping[key][1] or "finite_only") for key in components))
            ref,rid,rw,nobs=values(fn,n,cols,mode,circular)
            x,ix,w,iobs=values(fi,c,cols,mode,circular)
            native_cells=len(np.unique(rid))
            row.update(native_cells=native_cells,native_observations=nobs,introduced_observations=iobs)
            if len(ref)==0 or len(x)==0:
                output["tests"].append(dict(**row,status="source_unavailable",pairs=0,blocks=0,p=1.))
                continue
            a,sa,_=design_matrices(ref,rid,len(n),ref,rw,circular)
            b,sb,_=design_matrices(x,ix,len(c),ref,rw,circular)
            finite=np.isfinite(b).all(1)
            matched=[(l[finite[l]&finite[r]],r[finite[l]&finite[r]]) for l,r in pairs]
            matched=[(l,r) for l,r in matched if len(l)]
            nn=sum(len(l) for l,r in matched)
            nleft=sum(finite[l].sum() for l,r in pairs);nright=sum(finite[r].sum() for l,r in pairs)
            row.update(pairs=nn,blocks=len(matched),pair_retention=nn/total if total else 0,
                       left_available=int(nleft),right_available=int(nright))
            groups=[("native",sa[np.isfinite(sa).all(1)])]
            if nn:
                li=np.concatenate([l for l,r in matched]);ri=np.concatenate([r for l,r in matched])
                groups.extend([("analogue",sb[li]),("novel",sb[ri])])
            for name,m in groups:
                if len(m):
                    s=circular_summary(m.mean(0)) if circular else scalar_summary(m.mean(0))
                    output["profiles"].append(dict(**row,group=name,cells=len(m),**s))
            # Keep nonqualified endpoints in the Holm family with conservative p=1.
            if native_cells<50 or len(matched)<8 or nn<20 or nn<.8*total:
                output["tests"].append(dict(**row,status="insufficient_coverage",p=1.))
                continue
            counts=np.array([len(l) for l,r in matched])
            blocksum=np.stack([(b[r]-b[l]).sum(0) for l,r in matched])
            effect,p,exact,nperm=permute(blocksum,counts)
            output["tests"].append(dict(**row,status="estimated",distribution_distance=effect,p=p,
                                        exact=exact,permutations=nperm))
            # Signed rank shift and raw means/tail proportions remain descriptive; no cherry-picked tests.
            delta=sb[ri].mean(0)-sb[li].mean(0)
            if not circular:
                ds=np.stack([(sb[r]-sb[l]).sum(0) for l,r in matched])
                rng=np.random.default_rng(SEED)
                bs=rng.multinomial(len(counts),np.full(len(counts),1/len(counts)),size=3999)
                bd=(bs@ds)/(bs@counts)[:,None]
                for j,(label,val) in enumerate(scalar_summary(delta).items()):
                    lo,hi=np.quantile(bd[:,j],[.025,.975])
                    output["contrasts"].append(dict(**row,endpoint=label,difference=val,lo=lo,hi=hi))

def run_species(sp,sn,si,transform,edges,nat,inv,mapping,output):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860,30,183),"Cirsium vulgare":(1113,1139,146,316)}[sp]
    assert (len(n),len(full))==expected[:2]
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    on=_probability_content(nat,.95)
    labeled["novel"]=(~on[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c);assert sum(len(a) for a,b in pairs)==expected[2]
    run_scan(sp,"native_core95",c,pairs,fn,fi,n,mapping,output)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c);assert sum(len(a) for a,b in pairs)==expected[3]
    run_scan(sp,"density_surplus",c,pairs,fn,fi,n,mapping,output)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,
        source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured,mapping=load_measures(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    output={k:[] for k in ["tests","profiles","contrasts"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        target=frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"]
        assert abs(1-abs(nat-inv).sum()/2-target)<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,mapping,output)
    for name,rows in output.items():
        d=pd.DataFrame(rows)
        if name=="tests":
            assert len(d)==2*len(FEATURES)*len(DEFS)*len(MODES)
            assert (d.groupby(["definition","representation"]).size()==2*len(FEATURES)).all()
            d["holm"]=d.groupby(["definition","representation"]).p.transform(lambda p:multipletests(p.fillna(1),method="holm")[1])
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
