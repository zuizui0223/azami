#!/usr/bin/env python3
"""Observer and dominant-taxon sensitivity for the all-source native azami v4.
Private observer IDs and exact coordinates NEVER leave local working inputs.
Only aggregate coefficients and denominator audits can be published.
This reuses historically measured photo traits: it does NOT remeasure photographs.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import hashlib
import zipfile
import numpy as np
import pandas as pd
import statsmodels.api as sm

ORIGINAL_SHA = "2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae"
ENV = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
       "chelsa_rsds_mean", "chelsa_vpd_mean", "chelsa_sfcwind_mean",
       "chelsa_gsp", "chelsa_npp"]
TRAITS = {"colour": ("corolla_lab_chroma_median", "chelsa_rsds_mean"),
          "presentation": ("orientation_angle_degrees_median", "chelsa_bio12")}
COLOUR_COMPONENTS = ["corolla_lab_lightness_median", "corolla_lab_chroma_median",
                     "corolla_hue_sin_median", "corolla_hue_cos_median"]
ALL_TRAITS = ["orientation_angle_degrees_median", *COLOUR_COMPONENTS]
SCHEMES = ("observation", "equal_taxon", "equal_taxon_cell", "equal_taxon_cell_observer")

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def observer_source(metadata_zip, native_roster, output):
    """Extract one original observer per observation, only into a PRIVATE table."""
    if digest(metadata_zip) != ORIGINAL_SHA:
        raise ValueError("Wrong original pre-detector metadata snapshot")
    with zipfile.ZipFile(metadata_zip) as z, z.open("photo_metadata_merged.csv") as member:
        rows = pd.read_csv(member, usecols=["obs_id", "photo_id", "user_id"],
                           dtype={"obs_id": str, "photo_id": str, "user_id": str})
    if len(rows) != 1122854 or rows.obs_id.nunique() != 665115:
        raise ValueError("Unexpected original 665115-observation source")
    if rows.photo_id.nunique() != len(rows):
        raise ValueError("Duplicate photo identifier")
    if rows.groupby("obs_id").user_id.nunique(dropna=False).max() != 1:
        raise ValueError("Observer identity conflicts among photographs")
    u = rows.drop_duplicates("obs_id")[["obs_id", "user_id"]]
    n = pd.read_csv(native_roster, dtype={"obs_id": str})
    if len(n) != 313350 or not n.obs_id.is_unique:
        raise ValueError("Expected 313350 native observations")
    out = n[["obs_id"]].merge(u, on="obs_id", how="left", validate="one_to_one")
    if out.user_id.isna().any():
        raise ValueError("Missing observer ID: cannot invent independent observers")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return {"original_observations": 665115, "native_observations": len(out),
            "native_observers": int(out.user_id.nunique()),
            "private_table_sha256": digest(output)}

def cells(lat, lon):
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    if not (np.isfinite(lat).all() and np.isfinite(lon).all()):
        raise ValueError("Missing native coordinates")
    return (np.floor((np.sin(np.deg2rad(lat)) + 1) * 45).clip(0, 89).astype(int) * 180
            + np.floor(((lon + 180) % 360) / 2).astype(int))

def spatial(lat, lon):
    lat = np.deg2rad(np.asarray(lat, float))
    lon = np.deg2rad(np.asarray(lon, float))
    x, y, z = np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)
    return np.column_stack((x, y, z, x*y, x*z, y*z, x*x-y*y, 3*z*z-1))

def weights(f, method):
    if method not in SCHEMES:
        raise ValueError("Unknown observer weighting rule")
    if method == "observation":
        return np.full(len(f), 1/len(f))
    taxon_n = f.groupby("taxon").obs_id.transform("size").to_numpy(float)
    if method == "equal_taxon":
        raw = 1/taxon_n
    else:
        n_cell = f.groupby(["taxon", "cell"]).obs_id.transform("size").to_numpy(float)
        cells_per_taxon = f.groupby("taxon").cell.transform("nunique").to_numpy(float)
        if method == "equal_taxon_cell":
            raw = 1/(n_cell*cells_per_taxon)
        else:
            if f.user_id.isna().any():
                raise ValueError("Missing original user ID")
            repeated = f.groupby(["taxon", "cell", "user_id"]).obs_id.transform("size").to_numpy(float)
            observers_per_cell = f.groupby(["taxon", "cell"]).user_id.transform("nunique").to_numpy(float)
            raw = 1/(repeated*observers_per_cell*cells_per_taxon)
    return raw/raw.sum()

def means(a, codes, w, n):
    a = np.asarray(a, float)
    if a.ndim == 1:
        a = a[:, None]
    den = np.bincount(codes, weights=w, minlength=n)
    if (den<=0).any():
        raise ValueError("Zero taxon weight")
    return np.column_stack([np.bincount(codes, weights=w*a[:,k], minlength=n)/den
                            for k in range(a.shape[1])])

def standardized(a, w=None):
    a = np.asarray(a, float)
    if a.ndim == 1:
        a = a[:, None]
    w = np.full(len(a), 1/len(a)) if w is None else np.asarray(w)/sum(w)
    mean = np.average(a, weights=w, axis=0)
    sd = np.sqrt(np.average((a-mean)**2, weights=w, axis=0))
    if (sd<=1e-12).any():
        raise ValueError("Constant response or predictor")
    return (a-mean)/sd

def fit(f, method, response, predictor, scale, adjustment="marginal"):
    names = [predictor] if adjustment == "marginal" else [predictor] + [x for x in ENV if x!=predictor]
    z = f.dropna(subset=[response, *names]).copy()
    z = z.loc[z.groupby("taxon").obs_id.transform("size")>=5].reset_index(drop=True)
    if z.taxon.nunique()<20:
        raise ValueError("Fewer than 20 measured taxa")
    w = weights(z, method)
    X = z[names].to_numpy(float)
    if adjustment == "joint_spatial":
        X = np.column_stack((X, spatial(z.latitude, z.longitude)))
    y = z[response].to_numpy(float)
    if scale == "among":
        codes, taxa = pd.factorize(z.taxon, sort=True)
        values = means(np.column_stack((y, X)), codes, w, len(taxa))
        yy = standardized(values[:, 0])[:, 0]
        xx = standardized(values[:, 1:])
        result = sm.OLS(yy, sm.add_constant(xx)).fit(cov_type="HC3")
        index = 1
    elif scale == "within":
        codes, taxa = pd.factorize(z.taxon, sort=True)
        mx = means(X, codes, w, len(taxa))
        my = means(y, codes, w, len(taxa))[:, 0]
        xx = X-mx[codes]
        yy = y-my[codes]
        sd = np.sqrt(np.average(xx**2, weights=w, axis=0))
        sy = np.sqrt(np.average(yy**2, weights=w))
        if np.any(sd<=1e-12) or sy<=1e-12:
            raise ValueError("No within-taxon predictor support")
        xx = xx/sd
        yy = yy/sy
        result = sm.WLS(yy, xx, weights=w).fit(cov_type="cluster",
                                               cov_kwds={"groups":z.taxon})
        index = 0
    elif scale == "pooled" and adjustment=="marginal":
        yy = standardized(y, w)[:, 0]
        xx = standardized(X, w)
        result = sm.WLS(yy, xx, weights=w).fit(cov_type="cluster",
                                               cov_kwds={"groups":z.taxon})
        index = 0
    else:
        raise ValueError("Unknown analysis scale")
    conf = np.asarray(result.conf_int())
    return {"scale":scale, "weighting":method, "model":adjustment,
            "n_observations":len(z), "n_taxa":int(z.taxon.nunique()),
            "beta":float(result.params[index]),
            "ci95":[float(t) for t in conf[index]],
            "r2":float(result.rsquared),
            "condition_number":float(result.condition_number)}

def build_inputs(roster, observer, measured, climate):
    if len(roster)!=313350 or roster.obs_id.nunique()!=len(roster):
        raise ValueError("Invalid 313350-record source")
    if observer.obs_id.nunique()!=len(observer):
        raise ValueError("Duplicate observer observation")
    for d in [measured, climate]:
        if not d.obs_id.is_unique:
            raise ValueError("Duplicate measurement/environment")
    r = roster.merge(observer[["obs_id", "user_id"]], on="obs_id",
                     how="left", validate="one_to_one")
    if r.user_id.isna().any():
        raise ValueError("Unknown observer in native roster")
    r["taxon"] = r.accepted_key.astype("int64").astype(str)
    r["cell"] = cells(r.latitude, r.longitude)
    m = r.merge(measured[["obs_id", "taxon_name", *ALL_TRAITS]], on="obs_id",
                how="left", validate="one_to_one", suffixes=("","_measured"))
    wrong = m.taxon_name_measured.notna() & (m.taxon_name != m.taxon_name_measured)
    if wrong.any():
        raise ValueError("Taxon name disagrees with original measurement")
    m = m.drop(columns=["taxon_name_measured"])
    m = m.merge(climate[["obs_id",*ENV]], on="obs_id", validate="one_to_one")
    return r, m

def analysis(roster, observer, measured, climate):
    source, frame = build_inputs(roster,observer,measured,climate)
    total = len(source)
    taxa = source.taxon.value_counts()
    users = source.user_id.value_counts()
    if source.taxon.nunique()!=306 or source.user_id.nunique()<90000:
        raise ValueError("Unexpected native source observer/taxon coverage")
    group_repeat = source.groupby(["taxon","cell","user_id"]).obs_id.transform("size")
    top1 = users.head(round(len(users)*.01)).index
    top10 = users.head(round(len(users)*.1)).index
    audit = {"all_native_observations":total, "accepted_taxa":int(len(taxa)),
             "unique_observers":int(len(users)),
             "top2_species_observation_fraction":float(taxa.iloc[:2].sum()/total),
             "top5_species_observation_fraction":float(taxa.iloc[:5].sum()/total),
             "simpson_effective_species_observation_weighted":float(1/sum((taxa/total)**2)),
             "fraction_repeated_taxon_cell_observer":float((group_repeat>1).mean()),
             "top1pct_observers_count":int(len(top1)),
             "top1pct_observers_observation_fraction":float(users.head(len(top1)).sum()/total),
             "top10pct_observers_observation_fraction":float(users.head(len(top10)).sum()/total)}
    rows=[]
    for trait,(response,predictor) in TRAITS.items():
        for weighting in SCHEMES:
            for scale in ["pooled","within","among"]:
                for adjust in (["marginal"] if scale=="pooled"
                               else ["marginal","joint_environment","joint_spatial"]):
                    a=fit(frame,weighting,response,predictor,scale,adjust)
                    rows.append(dict(a,trait=trait,omission="none"))
            for omit, targets in [("top2_species",taxa.head(2).index),
                                  ("top5_species",taxa.head(5).index)]:
                d=frame.loc[~frame.taxon.isin(targets)]
                for scale in ["pooled","within","among"]:
                    a=fit(d,weighting,response,predictor,scale,"marginal")
                    rows.append(dict(a,trait=trait,omission=omit))
        # Fix the *species support* before comparing observer omission.
        for percent,selected in [(1,top1),(10,top10)]:
            clean=frame.loc[~frame.user_id.isin(selected)].copy()
            valid=clean.dropna(subset=[response,*ENV])
            eligible=valid.loc[valid.groupby("taxon").obs_id.transform("size")>=5,"taxon"].unique()
            control=frame.loc[frame.taxon.isin(eligible)]
            for sub,name in [(control,"same_taxa_before_observer_omission"),
                             (clean,"removed_top_observers")]:
                for adjust in ["marginal","joint_spatial"]:
                    a=fit(sub,"equal_taxon_cell_observer",response,predictor,"among",adjust)
                    rows.append(dict(a,trait=trait,omission=f"top{percent}pct_{name}"))
    return audit, rows

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="command", required=True)
    get=sub.add_parser("prepare")
    for arg in ["metadata-zip","native-roster","out-private-csv"]:
        get.add_argument("--"+arg, required=True)
    go=sub.add_parser("analyse")
    for arg in ["native-roster","observer-private-csv","archived-traits",
                "native-climate","out-dir"]:
        go.add_argument("--"+arg, required=True)
    a=p.parse_args()
    if a.command=="prepare":
        print(json.dumps(observer_source(a.metadata_zip,a.native_roster,a.out_private_csv)))
        return
    native=pd.read_csv(a.native_roster, dtype={"obs_id":str})
    observer=pd.read_csv(a.observer_private_csv,dtype={"obs_id":str,"user_id":str})
    measured=pd.read_csv(a.archived_traits,dtype={"obs_id":str},
                         usecols=["obs_id","taxon_name",*ALL_TRAITS])
    climate=pd.read_csv(a.native_climate,dtype={"obs_id":str},
                        usecols=["obs_id",*ENV])
    receipt, estimates=analysis(native,observer,measured,climate)
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    (out/"observer_bias_summary.json").write_text(
        json.dumps({"source":receipt,"fits":estimates,
        "limitations":["Exploratory, original photo observations not ecological census",
                       "Previously measured traits only; missing photos remain missing",
                       "Native means WCVP/TDWG3 regional agreement, not point ancestry",
                       "Species-level versus observation-level effects are different estimands",
                       "Photographer and region balance is a sensitivity, not formal detection correction",
                       "Joint space/environment coefficients remain confounded",
                       "No individual user IDs or coordinates are in this aggregate receipt"]},
                   ensure_ascii=False,indent=2)+"\n")
    pd.DataFrame(estimates).to_csv(out/"observer_bias_fits.csv",index=False)
    print(json.dumps(receipt,indent=2))
if __name__=="__main__":
    main()
