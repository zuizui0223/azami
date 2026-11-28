"""Prepare the frozen cell-median cohorts for the R hypervolume package.
Private inputs stay inside the runner; only modeled random points/aggregates
are exported by the R stage. Do not export this input file.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from analysis.v5.part2_hypervolume import (
    ROOT, EXPECTED, VARIABLES, T, authority_source, angular_cells, background_mesh,
    apply_original_five_CHELSA, make_cells, match_cells, head_cells,
    _background_species_audit, _pca_env_2d, _probability_content)
from analysis.v5.part2_novelty_cloud import make_joint_cells

OUT=ROOT/"private/r_hypervolume_inputs.csv"
AXES_TRAIT=("angle","chroma","outline")
AXES_CLIMATE=tuple(VARIABLES)

def add(out,scales,sp,space,design,group,frame,cols,axes,mu,sd):
    a=frame[list(cols)].to_numpy(float)
    assert len(a)>0 and np.isfinite(a).all() and np.all(sd>0)
    z=(a-mu)/sd
    d=pd.DataFrame(z,columns=[f"x{i+1}" for i in range(z.shape[1])])
    d["species"]=sp;d["space"]=space;d["design"]=design;d["group"]=group
    out.append(d)
    scales.append(dict(species=sp,space=space,design=design,group=group,
                       cells=len(a),blocks=int(frame.block.nunique())))
    
def main():
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=pd.read_csv(ROOT/"private/original-measured.csv.gz",usecols=["obs_id",*T],dtype={"obs_id":str})
    wild=wild.drop(columns=T,errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out=[];cohorts=[];scale_rows=[]
    for sp,rows in wild.groupby("taxon_name"):
        f=rows.dropna(subset=list(VARIABLES))
        sn=f.loc[f.range_state=="native"];si=f.loc[f.range_state=="introduced"]
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        en=env[env._range=="native"];ei=env[env._range=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(en,ei,sn,si,return_grid=True)
        assert abs(np.minimum(nat,inv).sum()-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        n=make_cells(sn,True);full=make_cells(si,True);ex=EXPECTED[sp]
        assert (len(n),len(full))==ex[:2]
        a=n[T].to_numpy(float);mu=a.mean(0);sd=a.std(0)
        for k,name in enumerate(AXES_TRAIT):
            scale_rows.append(dict(species=sp,space="trait",axis=name,center=mu[k],sd=sd[k]))
        def trait(design,group,frame):
            add(out,cohorts,sp,"trait",design,group,frame,T,AXES_TRAIT,mu,sd)
        trait("native_all","native",n);trait("native_all","introduced",full)
        coords=transform(si)
        uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
        for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
        inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
        lab=si.loc[inside].copy()
        lab["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
        c=make_cells(lab);c=c.loc[c.novel!=.5].reset_index(drop=True)
        for label,frame in [("native",n),("novel",c[c.novel==1])]:
            trait("native_novel",label,frame)
        def matched(c,definition,expected):
            pairs=match_cells(c)
            assert sum(len(l) for l,r in pairs)==expected
            left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
            assert len(np.unique(np.r_[left,right]))==2*expected
            assert c.iloc[left].block.value_counts().sort_index().equals(c.iloc[right].block.value_counts().sort_index())
            trait(definition,"analogous",c.iloc[left]);trait(definition,"novel",c.iloc[right])
        matched(c,"matched_native_core95",ex[2])
        hc=head_cells(si,transform,edges)[["i","j","u","v"]]
        c=full.merge(hc,on=["i","j"],validate="one_to_one")
        u,v=c.u.to_numpy(),c.v.to_numpy()
        c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
        c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
        matched(c,"matched_density_surplus",ex[3])
        # Five-dimensional continuous R climate comparison is raw observed
        # climate of the identical complete-head cells, not background corrected.
        x,y=make_joint_cells(sn),make_joint_cells(si)
        assert (len(x),len(y))==ex[:2]
        z=x[list(VARIABLES)].to_numpy(float);cm=z.mean(0);cs=z.std(0)
        for k,name in enumerate(AXES_CLIMATE):
            scale_rows.append(dict(species=sp,space="climate_raw",axis=name,center=cm[k],sd=cs[k]))
        for label,frame in [("native",x),("introduced",y)]:
            add(out,cohorts,sp,"climate_raw","native_all",label,frame,VARIABLES,AXES_CLIMATE,cm,cs)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    pd.concat(out,ignore_index=True).to_csv(OUT,index=False)
    public=ROOT/"part2_r_hypervolume";public.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(cohorts).to_csv(public/"cohorts.csv",index=False)
    pd.DataFrame(scale_rows).to_csv(public/"scales.csv",index=False)
    print("Verified exact complete-head and geographic-pair cohort sizes; wrote private analysis inputs.",flush=True)
if __name__=="__main__":
    main()
