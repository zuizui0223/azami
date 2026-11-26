"""Exploratory native-reference phenotype frequency redistribution.

No individual observations, coordinates, observer IDs or block IDs are exported.
Climate-only group definitions precede all phenotype comparisons. The native
reference is contemporary, not a sampled ancestral population.
"""
from pathlib import Path
import json
from itertools import product
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.native_pca_env import (
    background_mesh, _background_species_audit, _pca_env_2d, _probability_content
)
from analysis.v5.invasion_partition import background, head_cells
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES

ROOT = Path("work/submission")
OUT = ROOT / "redistribution"
Q = np.array([.025,.20,.40,.60,.80,.975])
TRAITS = dict(zip(("angle","chroma","outline"), T))
LABELS = ("below_q025","q025_q20","q20_q40","q40_q60","q60_q80","q80_q975","above_q975")
SEED = 20261126


def eligible(frame):
    d = frame.loc[frame[T[0]].between(0,180)].dropna(subset=T).copy()
    d["i"] = np.minimum(179, np.floor(d.latitude+90)).astype(int)
    d["j"] = np.floor((d.longitude+180)%360).astype(int)
    d["cellkey"] = d.i * 360 + d.j
    return d


def cellkeys(c):
    return (c.i*360+c.j).to_numpy()


def weighted_quantile(x, w, qs):
    order = np.argsort(x, kind="stable")
    x,w = np.asarray(x)[order],np.asarray(w)[order]
    cumulative = np.cumsum(w)/w.sum()
    return x[np.minimum(np.searchsorted(cumulative,qs,side="left"),len(x)-1)]


def rows_for(frame,cells,rep):
    if rep == "median":
        return cells[T].to_numpy(),np.arange(len(cells)),np.ones(len(cells))
    ids = pd.Index(cellkeys(cells)).get_indexer(frame.cellkey)
    ok=ids>=0
    ix=ids[ok]
    assert len(ix)>0 and np.all(np.bincount(ix,minlength=len(cells))>0)
    return frame.loc[ok,T].to_numpy(),ix,1/np.bincount(ix,minlength=len(cells))[ix]


def cell_hist(codes,ids,n,k):
    h=np.zeros((n,k),float)
    np.add.at(h,(ids,codes),1)
    den=h.sum(1)
    assert (den>0).all()
    return h/den[:,None]


def cell_average(x,ids,n):
    return np.bincount(ids,weights=x,minlength=n)/np.bincount(ids,minlength=n)


def encode(values,ids,n,cuts,cols):
    codes=np.zeros(len(values),int)
    for j,c in zip(cols,cuts):
        codes=codes*(len(c)+1)+np.searchsorted(c,values[:,j],side="right")
    return cell_hist(codes,ids,n,int(np.prod([len(c)+1 for c in cuts])))


def tv(a,b):
    return .5*np.abs(a-b).sum(axis=-1)


def randomization(a,b,ref,counts,rng):
    """a,b are block sums, not independent cells. Preserve paired block size."""
    k=len(counts); denom=counts.sum()
    left,right=a.sum(0)/denom,b.sum(0)/denom
    diff=right-left
    exact=k<=15
    flips=np.array(list(product([0,1],repeat=k)),dtype=float) if exact else rng.integers(0,2,(9999,k))
    shift=flips@(b-a)/denom
    ls=left+shift; rs=right-shift
    diffs=rs-ls
    p=(np.sum(np.abs(diffs)>=np.abs(diff)-1e-12,axis=0)+int(not exact))/(len(flips)+int(not exact))
    # Resample paired blocks. The reference is fixed; this is conditional uncertainty.
    bs=rng.multinomial(k,np.full(k,1/k),size=3999)
    bd=(bs@(b-a))/(bs@counts)[:,None]
    lo,hi=np.quantile(bd,[.025,.975],axis=0)
    return left,right,diff,p,lo,hi,ls,rs,exact,len(flips)


def paired_features(a,b,ref,pairs,seed):
    counts=np.array([len(x) for x,y in pairs])
    aa=np.stack([a[x].sum(0) for x,y in pairs])
    bb=np.stack([b[y].sum(0) for x,y in pairs])
    return randomization(aa,bb,ref,counts,np.random.default_rng(seed))


def probability_p(null,observed,exact):
    return (np.count_nonzero(np.abs(null)>=abs(observed)-1e-12)+int(not exact))/(len(null)+int(not exact))


def distributions(sp,cover,rep,frame_n,frame_i,n,c,pairs,avail,all_i,
                  summary,bins,tests,thresholds,profiles):
    vn,ni,nw=rows_for(frame_n,n,rep)
    vi,ii,iw=rows_for(frame_i,c,rep)
    va,ai,aw=rows_for(frame_i,all_i,rep)
    matches_a=np.concatenate([a for a,b in pairs]) if pairs else np.array([],int)
    matches_b=np.concatenate([b for a,b in pairs]) if pairs else np.array([],int)
    info=dict(species=sp,coverage=cover,representation=rep,n_native=len(n),
              n_sim=int((c.novel==0).sum()),n_new=int((c.novel==1).sum()),
              pairs=len(matches_a),blocks=len(pairs))
    valid=len(pairs)>=3 and len(matches_a)>=8
    def emit(group,kind,native,inter,whole,labels,kind_tests):
        ref=native.mean(0)
        masks=[("native",native,np.ones(len(n),bool)),
               ("introduced_all",whole,np.ones(len(all_i),bool)),
               ("analogue_all",inter,c.novel.to_numpy()==0),
               ("novel_all",inter,c.novel.to_numpy()==1)]
        if valid:
            masks += [("analogue_matched",inter,matches_a),("novel_matched",inter,matches_b)]
        if cover==.95:
            masks += [(g,whole,m) for g,m in avail.items()]
        for name,mat,idx in masks:
            x=mat[idx]
            if len(x)<8: continue
            for label,v,r in zip(labels,x.mean(0),ref):
                summary.append(dict(**info,feature=group,kind=kind,group=name,
                                    category=label,cells=len(x),value=v,native_reference=r))
        if not valid: return
        l,r,d,p,lo,hi,ls,rs,exact,nperm=paired_features(inter,inter,ref,pairs,SEED)
        for label,v,z,a,b,pv,low,high in zip(labels,l,r,ref,d,p,lo,hi):
            bins.append(dict(**info,feature=group,kind=kind,category=label,native_reference=a,
                             matched_analogue=v,matched_novel=z,difference=b,lo=low,hi=high,
                             p=pv,exact=exact,permutations=nperm,family=kind_tests))
        if kind in ("marginal","joint"):
            observed=tv(l,r)
            departure=tv(r,ref)-tv(l,ref)
            for test,effect,null in [("distribution_TV",observed,tv(ls,rs)),
                                    ("departure_from_native",departure,tv(rs,ref)-tv(ls,ref))]:
                tests.append(dict(**info,feature=group,kind=kind,test=test,effect=effect,
                                  p=probability_p(null,effect,exact),exact=exact,permutations=nperm))
    for j,(name,col) in enumerate(TRAITS.items()):
        qs=np.quantile(vn[:,j],Q) if rep=="median" else weighted_quantile(vn[:,j],nw,Q)
        minimum,maximum=vn[:,j].min(),vn[:,j].max()
        thresholds.append(dict(**info,feature=name,minimum=minimum,maximum=maximum,
                               **{f"q{q:g}":v for q,v in zip(Q,qs)}))
        vals=[encode(v,idx,k,[qs],[j]) for v,idx,k in [(vn,ni,len(n)),(vi,ii,len(c)),(va,ai,len(all_i))]]
        emit(name,"marginal",*vals,LABELS,"marginal_bins")
        # Explicitly distinguish observed native support from statistical tails.
        def special(v,idx,k):
            x=v[:,j]
            return np.column_stack([cell_average(z,idx,k) for z in
                [(x>=minimum)&(x<=maximum)&(x<qs[1]),
                 (x>=minimum)&(x<=maximum)&(x>=qs[-2]),
                 x<minimum,x>maximum,x]])
        vals=[special(v,idx,k) for v,idx,k in [(vn,ni,len(n)),(vi,ii,len(c)),(va,ai,len(all_i))]]
        emit(name,"support_and_mean",*vals,("existing_low20","existing_high20","below_observed_min","above_observed_max","mean_raw"),"support_and_mean")
    for cols,label in [([0,1],"angle+chroma"),([0,1,2],"angle+chroma+outline")]:
        cuts=[np.quantile(vn[:,j],[1/3,2/3]) if rep=="median" else
              weighted_quantile(vn[:,j],nw,[1/3,2/3]) for j in cols]
        vals=[encode(v,idx,k,cuts,cols) for v,idx,k in [(vn,ni,len(n)),(vi,ii,len(c)),(va,ai,len(all_i))]]
        labs=["/".join(p) for p in product(("low","mid","high"),repeat=len(cols))]
        emit(label,"joint",*vals,labs,"joint_bins")
        ref=vals[0].mean(0)
        classes=np.array([0 if p==0 else 1 if p<=.02 else 2 for p in ref])
        states=[np.column_stack([v[:,classes==k].sum(1) for k in range(3)]) for v in vals]
        emit(label,"joint_support",*states,("unobserved_native_combination","rare_existing_combination","common_existing_combination"),"joint_support")


def climate_profile(sp,coverage,n,c,pairs,frame_n,frame_i,out):
    def climate(frame,cells):
        q=frame.groupby("cellkey")[list(VARIABLES)].median()
        return q.reindex(cellkeys(cells)).to_numpy()
    a,b=climate(frame_n,n),climate(frame_i,c)
    groups=[("native",a),("analogous",b[c.novel==0]),("novel",b[c.novel==1])]
    if pairs:
        groups += [("analogous_matched",b[np.concatenate([x for x,y in pairs])]),
                   ("novel_matched",b[np.concatenate([y for x,y in pairs])])]
    for name,v in groups:
        if len(v)<8:continue
        for j,axis in enumerate(VARIABLES):
            lo,mid,hi=np.nanquantile(v[:,j],[.25,.5,.75])
            out.append(dict(species=sp,coverage=coverage,group=name,cells=len(v),
                            axis=axis,q25=lo,median=mid,q75=hi))


def run_species(sp,sn,si,env,transform,edges,nat,inv,outputs):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True)
    full_i=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860,30),"Cirsium vulgare":(1113,1139,146)}[sp]
    assert (len(n),len(full_i))==expected[:2],(sp,len(n),len(full_i))
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    bn,bi=[background(env[env._range==k],transform,edges) for k in ["native","introduced"]]
    an,ai=[_probability_content(v,.99) for v in (bn,bi)]
    hc=head_cells(si,transform,edges)
    ix=pd.Index(cellkeys(hc)).get_indexer(cellkeys(full_i))
    keep=ix>=0
    hu,hv=hc.u.to_numpy()[ix[keep]],hc.v.to_numpy()[ix[keep]]
    for cover in [.9,.95,.99]:
        on=_probability_content(nat,cover)
        labeled=si.loc[inside].copy()
        labeled["novel"]=(~on[uv[inside,0],uv[inside,1]]).astype(int)
        c=make_cells(labeled)
        c=c.loc[c.novel!=.5].reset_index(drop=True)
        pairs=match_cells(c)
        if cover==.95:
            assert sum(len(a) for a,b in pairs)==expected[2],(sp,[(len(a),len(b)) for a,b in pairs])
        avail={}
        for name,mask in [
            ("availability_core",ai&an&on),
            ("availability_expansion",ai&an&~on),
            ("availability_outside_native_background",ai&~an),
            ("availability_outside_intro_background",~ai)]:
            a=np.zeros(len(full_i),bool);a[keep]=mask[hu,hv];avail[name]=a
        climate_profile(sp,cover,n,c,pairs,fn,fi,outputs["climate_profiles"])
        for rep in ["median","record_equalcell"]:
            distributions(sp,cover,rep,fn,fi,n,c,pairs,avail,full_i,
                          outputs["frequencies"],outputs["paired_categories"],
                          outputs["omnibus"],outputs["thresholds"],outputs["climate_profiles"])
        outputs["cohort"].append(dict(species=sp,coverage=cover,native_cells=len(n),
            introduced_cells=len(full_i),analogous_cells=int((c.novel==0).sum()),
            novel_cells=int((c.novel==1).sum()),pairs=sum(len(a) for a,b in pairs),
            blocks=len(pairs),climate_D=1-abs(nat-inv).sum()/2,
            introduced_density_outside_native=inv[~on].sum()))


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
    measured=pd.read_csv(ROOT/"private/original-measured.csv.gz",usecols=["obs_id",*T],dtype={"obs_id":str})
    wild=wild.drop(columns=T,errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    outputs={k:[] for k in ["frequencies","paired_categories","omnibus","thresholds","climate_profiles","cohort"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        target=frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"]
        assert abs(1-abs(nat-inv).sum()/2-target)<1e-5
        run_species(sp,sn,si,env,transform,edges,nat,inv,outputs)
    for name,rows in outputs.items():
        d=pd.DataFrame(rows)
        if name in ("paired_categories","omnibus") and len(d):
            keys=["coverage","representation", "family" if name=="paired_categories" else "test"]
            d["holm"]=d.groupby(keys).p.transform(lambda p:multipletests(p,method="holm")[1])
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)


if __name__=="__main__":
    main()
