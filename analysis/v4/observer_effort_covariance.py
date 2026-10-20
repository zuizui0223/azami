"""Retrospective native-range thistle observer-effort sensitivity.

Input: 1 record per iNaturalist observation, original observer ID, historical
image measurement, original-resolution CHELSA. Outputs only species-level and
aggregate statistics, never observer IDs, photos, coordinates or observation IDs.

This analysis quantifies observed *reporting effort*, not unobserved abundance
or a causal environmental effect.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

ENV = ("chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15",
       "chelsa_rsds_mean","chelsa_vpd_mean","chelsa_sfcwind_mean","chelsa_gsp","chelsa_npp")
RESPONSES = (("colour","corolla_lab_chroma_median","chelsa_rsds_mean"),
             ("presentation","orientation_angle_degrees_median","chelsa_bio12"))
SCHEMES = ("observation","equal_taxon","equal_taxon_region","equal_taxon_region_observer")


def digest(path):
    import hashlib
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def equal_area_cells(lat,lon):
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    if not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("Missing or invalid coordinates in native roster")
    return (np.floor((np.sin(np.deg2rad(lat))+1)*45).clip(0,89).astype(int)*180+
            np.floor(((lon+180)%360)/2).astype(int))


def observation_weights(f,scheme):
    if scheme not in SCHEMES:raise ValueError("Unknown weighting scheme")
    if scheme=="observation":
        return np.ones(len(f))/len(f)
    n=f.groupby("taxon").obs_id.transform("size").to_numpy(float)
    if scheme=="equal_taxon":raw=1/n
    else:
        keys=["taxon","cell"]
        nc=f.groupby(keys).obs_id.transform("size").to_numpy(float)
        kc=f.groupby("taxon").cell.transform("nunique").to_numpy(float)
        if scheme=="equal_taxon_region":raw=1/(nc*kc)
        else:
            if f.user_id.isna().any():
                raise ValueError("Missing observer ID cannot be treated as independent")
            no=f.groupby(keys).user_id.transform("nunique").to_numpy(float)
            nu=f.groupby(keys+["user_id"]).obs_id.transform("size").to_numpy(float)
            raw=1/(kc*no*nu)
    return raw/raw.sum()


def covariance_parts(f,scheme,response,predictor):
    """Exact between- and within-taxon contribution to ONE pooled coefficient."""
    if f[["taxon","cell","user_id",response,predictor]].isna().any().any():
        raise ValueError("Missing observer, taxon or phenotype")
    w=observation_weights(f,scheme)
    ids,taxa=pd.factorize(f.taxon,sort=True)
    n=len(taxa)
    q=np.bincount(ids,weights=w,minlength=n)
    if (q<=0).any():raise ValueError("No weight for eligible species")
    x=f[predictor].to_numpy(float);y=f[response].to_numpy(float)
    ax=np.bincount(ids,weights=w*x,minlength=n)/q
    ay=np.bincount(ids,weights=w*y,minlength=n)/q
    mx=np.dot(w,x);my=np.dot(w,y)
    vx=np.dot(w,(x-mx)**2);vy=np.dot(w,(y-my)**2)
    if vx<=0 or vy<=0:raise ValueError("Constant predictor/response")
    denom=np.sqrt(vx*vy)
    total=np.dot(w,(x-mx)*(y-my))
    between=np.dot(q,(ax-mx)*(ay-my))
    within=np.dot(w,(x-ax[ids])*(y-ay[ids]))
    if not np.isclose(total,between+within,atol=1e-10,rtol=1e-9):
        raise ValueError("Covariance identity violated")
    beta=np.corrcoef(ax,ay)[0,1]
    return {
        "scheme":scheme, "n_observations":len(f), "n_taxa":n,
        "pooled_standardized_coefficient":float(total/denom),
        "between_contribution":float(between/denom),
        "within_contribution":float(within/denom),
        "equal_taxon_species_mean_correlation":float(beta),
        "identity_residual":float(abs(total-between-within)/denom),
        "max_taxon_weight":float(q.max()),
        "top_two_taxon_weight_mass":float(sum(sorted(q,reverse=True)[:2])),
    }


def load_matched(native,observer,traits,climate):
    n=pd.read_csv(native,dtype={"obs_id":str})
    o=pd.read_csv(observer,dtype={"obs_id":str,"user_id":str},usecols=["obs_id","user_id"])
    t=pd.read_csv(traits,dtype={"obs_id":str},usecols=["obs_id","taxon_name",*[r[1] for r in RESPONSES]])
    c=pd.read_csv(climate,dtype={"obs_id":str},usecols=["obs_id",*ENV])
    for name,d in [("native",n),("observer",o),("traits",t),("climate",c)]:
        if not d.obs_id.is_unique:raise ValueError(name+" has duplicate observation IDs")
    if len(n)!=313350 or n.accepted_key.nunique()!=306:
        raise ValueError("Full native master denominator changed")
    if len(c)!=len(n) or set(c.obs_id)!=set(n.obs_id):
        raise ValueError("Native climate table must include every source observation")
    f=n.merge(o,on="obs_id",validate="one_to_one")
    if f.user_id.isna().any():raise ValueError("Some native records lack original observer")
    f["taxon"]=f.accepted_key.astype("int64").astype(str)
    f["cell"]=equal_area_cells(f.latitude,f.longitude)
    if f.taxon_name.isna().any():raise ValueError("Unresolved source taxon")
    if t.taxon_name.isna().any():raise ValueError("Historical phenotype taxon missing")
    frame=f.merge(t,on="obs_id",how="left",validate="one_to_one",suffixes=("","_measured"))
    mismatch=(frame.taxon_name_measured.notna()&
              frame.taxon_name.ne(frame.taxon_name_measured))
    if mismatch.any():raise ValueError("Historical image trait taxon differs from source")
    frame=frame.drop(columns="taxon_name_measured")
    frame=frame.merge(c,on="obs_id",validate="one_to_one")
    return f,frame


def analyze(native,observer,traits,climate,outdir):
    root=Path(outdir);root.mkdir(exist_ok=True,parents=True)
    source,frame=load_matched(native,observer,traits,climate)
    nt=len(source);taxon_n=source.taxon.value_counts();observer_n=source.user_id.value_counts()
    n_top=max(1,round(len(observer_n)*.01))
    high_users=observer_n.head(n_top).index
    high_ten_percent_users=observer_n.head(max(1,round(len(observer_n)*.10))).index
    groups=source.groupby(["taxon","cell","user_id"]).obs_id.transform("size")
    stats={
       "native_observations":nt,"accepted_taxa":source.taxon.nunique(),
       "original_observers":source.user_id.nunique(),
       "top_two_taxon_fraction":float(taxon_n.head(2).sum()/nt),
       "top_five_taxon_fraction":float(taxon_n.head(5).sum()/nt),
       "simpson_effective_taxa":float(1/np.sum((taxon_n/nt)**2)),
       "top_one_percent_observer_fraction":float(observer_n.head(n_top).sum()/nt),
       "top_ten_percent_observer_fraction":float(observer_n.head(max(1,round(len(observer_n)*.10))).sum()/nt),
       "repeated_observer_taxon_cell_fraction":float((groups>1).mean()),
       "observers_in_top_one_percent":int(n_top)
    }
    rows=[]
    for label,response,predictor in RESPONSES:
        original=frame.dropna(subset=[response,*ENV]).copy()
        cases={
           "all":original,
           "omit_top_two_taxa":original.loc[~original.taxon.isin(taxon_n.head(2).index)],
           "omit_top_five_taxa":original.loc[~original.taxon.isin(taxon_n.head(5).index)],
           "omit_top_one_percent_observers":original.loc[~original.user_id.isin(high_users)],
           "omit_top_ten_percent_observers":original.loc[~original.user_id.isin(high_ten_percent_users)]
        }
        for case,sub in cases.items():
            sub=sub.reset_index(drop=True)
            sub=sub.loc[sub.groupby("taxon").obs_id.transform("size")>=5].reset_index(drop=True)
            for scheme in SCHEMES:
                row=covariance_parts(sub,scheme,response,predictor)
                row.update(response_group=label,omission=case)
                rows.append(row)
    out={"scope":"posthoc native photo reporting-effort sensitivity",
         "source":stats,"models":rows,
         "input_digests":{"native":digest(native),"original_observers":digest(observer),
                          "archived_image_traits":digest(traits),"original_resolution_climate":digest(climate)},
         "limits":["The within/between split is an exact descriptive covariance identity, not a causal decomposition.",
                   "Observer/cell weights represent sensitivity to sampling effort, not true species abundance.",
                   "No new original flower photographs are measured; historical phenotype selection remains.",
                   "Native means WCVP/TDWG regional compatibility, not individually certified local provenance.",
                   "Taxon mean correlation is standardized and not a spatially adjusted causal slope.",
                   "No private IDs, source coordinates or observer identifiers are exported."]}
    pd.DataFrame(rows).to_csv(root/"observer_effort_covariance.csv",index=False)
    (root/"observer_effort_receipt.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    return out


def main():
    p=argparse.ArgumentParser()
    for k in ["native-roster","observer-metadata","archived-traits","native-chelsa","out-dir"]:
        p.add_argument("--"+k,required=True)
    a=p.parse_args()
    out=analyze(a.native_roster,a.observer_metadata,a.archived_traits,a.native_chelsa,a.out_dir)
    print(json.dumps({"source":out["source"],"rows":len(out["models"])},indent=2))

if __name__=="__main__":main()
