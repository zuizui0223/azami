"""Can the contemporary native climate-phenotype relationship predict the introduced shift?

Exploratory native-only transport diagnostic. No extrapolation beyond the
geographically held-out native climate distance envelope, no causal or genetic
claims. All private observation and spatial identifiers stay in the runner.
"""
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_redistribution import ROOT, SEED, eligible, cellkeys
from analysis.v5.phenotype_endpoint_scan import load_measures, values, design_matrices, permute
from analysis.v5.phenotype_native_calibration import native_splits
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.invasion_partition import head_cells
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d, _probability_content
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_integration import FEATURES

OUT=ROOT/"phenotype_native_transport"
K=20
MIN_DRAW_FRACTION=.8
MIN_PAIR_FRACTION=.5
NATIVE_MODE="record_equalcell"


def climate(frame,cells):
    return frame.groupby("cellkey")[list(VARIABLES)].median().reindex(cellkeys(cells)).to_numpy(float)


def neighbors(train,query,k=K):
    mu=train.mean(0);sd=train.std(0)
    if not np.isfinite(sd).all() or (sd<=0).any() or len(train)<k:
        raise ValueError("Degenerate climate reference")
    d,i=cKDTree((train-mu)/sd).query((query-mu)/sd,k=list(range(1,k+1)))
    return d.mean(1),i


def predict_native(train_x,train_y,cal_x,intro_x,k=K):
    d,ix=neighbors(train_x,np.vstack([cal_x,intro_x]),k)
    p=train_y[ix].mean(1)
    nc=len(cal_x)
    threshold=np.quantile(d[:nc],.95)
    return p[:nc],p[nc:],d[nc:]<=threshold,threshold


def paired_scan(sp,definition,n,c,pairs,fn,fi,splits,mapping,out):
    left=np.concatenate([a for a,b in pairs]);right=np.concatenate([b for a,b in pairs])
    block=np.concatenate([np.full(len(a),j,int) for j,(a,b) in enumerate(pairs)])
    assert len(left)==len(right) and len(left)>0
    xn,xi=climate(fn,n),climate(fi,c)
    assert np.isfinite(xn).all() and np.isfinite(xi).all()
    for name,components in FEATURES.items():
        circ=name=="hue"
        cols=["scan_"+x for x in components]
        nv,ni,nw,_=values(fn,n,cols,NATIVE_MODE,circ)
        iv,ii,iw,_=values(fi,c,cols,NATIVE_MODE,circ)
        base=dict(species=sp,definition=definition,trait=name,representation=NATIVE_MODE,
                  fixed_pairs=len(left),fixed_blocks=len(pairs),
                  source_columns=";".join(str(mapping[key][0] or "unavailable") for key in components),
                  usable_head_gate=";".join(str(mapping[key][1] or "finite_only") for key in components))
        if not len(nv) or not len(iv):
            out["tests"].append(dict(**base,status="missing_trait",p=1.))
            continue
        counts=np.zeros(len(left),int)
        obs_sum=pred_sum=None
        osum=psum=None
        errors=[];baseline=[];positive=0;draws=0;thresholds=[]
        for tr,ca in splits:
            rmask=np.isin(ni,tr)
            if np.unique(ni[rmask]).size<100:continue
            ref,rw=nv[rmask],nw[rmask]
            f_n,s_n,_=design_matrices(nv,ni,len(n),ref,rw,circ)
            f_i,s_i,_=design_matrices(iv,ii,len(c),ref,rw,circ)
            train=tr[np.isfinite(f_n[tr]).all(1)]
            cal=ca[np.isfinite(f_n[ca]).all(1)]
            if len(train)<100 or len(cal)<50:continue
            # Additional descriptive signed statistics; raw units are not pooled
            # with native-percentile summaries. Hue uses only circular moments.
            if circ:
                sn=s_n;si=s_i
            else:
                sn=s_n[:,1:4];si=s_i[:,1:4]
            train_y=np.column_stack([f_n[train],sn[train]])
            pc,pi,support,threshold=predict_native(xn[train],train_y,xn[cal],xi)
            p=f_n.shape[1]
            e=float(np.square(f_n[cal]-pc[:,:p]).sum())
            b=float(np.square(f_n[cal]-f_n[train].mean(0)).sum())
            if b<=1e-12:continue
            errors.append(e);baseline.append(b);positive+=int(e<b);draws+=1;thresholds.append(threshold)
            good=support[left]&support[right]&np.isfinite(f_i[left]).all(1)&np.isfinite(f_i[right]).all(1)
            if obs_sum is None:
                obs_sum=np.zeros((len(left),p));pred_sum=np.zeros_like(obs_sum)
                osum=np.zeros((len(left),si.shape[1]));psum=np.zeros_like(osum)
            obs_sum[good]+=f_i[right[good]]-f_i[left[good]]
            pred_sum[good]+=pi[right[good],:p]-pi[left[good],:p]
            osum[good]+=si[right[good]]-si[left[good]]
            psum[good]+=pi[right[good],p:]-pi[left[good],p:]
            counts[good]+=1
        if not draws:
            out["tests"].append(dict(**base,status="no_native_validation",p=1.))
            continue
        q2=1-sum(errors)/sum(baseline)
        good=counts>=MIN_DRAW_FRACTION*len(splits)
        groups=np.unique(block[good])
        nn=int(good.sum())
        base.update(native_draws=draws,native_cv_q2=q2,native_positive_cv_fraction=positive/draws,
                    mean_native_support_distance95=float(np.mean(thresholds)),
                    supported_pairs=nn,supported_blocks=len(groups),pair_retention=nn/len(left))
        if nn==0:
            out["tests"].append(dict(**base,status="no_supported_pairs",p=1.))
            continue
        observed=obs_sum[good]/counts[good,None]
        predicted=pred_sum[good]/counts[good,None]
        residual=observed-predicted
        o,p,r=observed.mean(0),predicted.mean(0),residual.mean(0)
        obs_rms=float(np.sqrt(np.mean(o**2)));pred_rms=float(np.sqrt(np.mean(p**2)))
        alignment=float(o@p/(np.linalg.norm(o)*np.linalg.norm(p))) if np.linalg.norm(o)*np.linalg.norm(p)>1e-12 else np.nan
        base.update(observed_rms=obs_rms,predicted_rms=pred_rms,residual_rms=float(np.sqrt(np.mean(r**2))),alignment=alignment)
        labels=["hue_sin","hue_cos"] if circ else ["mean_native_rank","existing_low20","existing_high20"]
        oo=(osum[good]/counts[good,None]).mean(0);pp=(psum[good]/counts[good,None]).mean(0)
        for j,label in enumerate(labels):
            out["contrasts"].append(dict(**base,endpoint=label,observed_difference=oo[j],native_prediction=pp[j],residual=oo[j]-pp[j]))
        if q2<=0 or positive/draws<.5:
            out["tests"].append(dict(**base,status="native_prediction_not_validated",p=1.))
            continue
        if len(groups)<8 or nn<20 or nn<MIN_PAIR_FRACTION*len(left):
            out["tests"].append(dict(**base,status="insufficient_common_climate_support",p=1.))
            continue
        # Conditional approximate whole-block sign test on residual contrast,
        # not a randomized assignment or a test of natural selection.
        bsum=np.stack([residual[block[good]==j].sum(0) for j in groups])
        n_per=np.array([np.sum(block[good]==j) for j in groups])
        effect,pv,exact,nperm=permute(bsum,n_per)
        out["tests"].append(dict(**base,status="estimated",p=pv,exhaustive_sign_patterns=exact,permutations=nperm))


def run_species(sp,sn,si,transform,edges,nat,inv,mapping,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860,30,183),"Cirsium vulgare":(1113,1139,146,316)}[sp]
    assert (len(n),len(full))==expected[:2]
    splits=native_splits(n.block.to_numpy())
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy();on=_probability_content(nat,.95)
    labeled["novel"]=(~on[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c);assert sum(len(a) for a,b in pairs)==expected[2]
    paired_scan(sp,"native_core95",n,c,pairs,fn,fi,splits,mapping,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c);assert sum(len(a) for a,b in pairs)==expected[3]
    paired_scan(sp,"density_surplus",n,c,pairs,fn,fi,splits,mapping,out)


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
    out={k:[] for k in ["tests","contrasts"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,mapping,out)
    for name,rows in out.items():
        d=pd.DataFrame(rows)
        if name=="tests":
            assert len(d)==32
            d["holm"]=d.groupby("definition").p.transform(lambda p:multipletests(p.fillna(1),method="holm")[1])
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
