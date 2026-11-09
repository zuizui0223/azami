"""Monthly normals at dated V5 observations; export only ecological aggregates."""
import hashlib
import io
import os
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import rasterio
import requests
import statsmodels.api as sm

from analysis.v4.native_introduced_opportunity import (
    decrypt_checked, ORIGINAL_SHA, ORIGINAL_SEALED_SHA, ORIGINAL_AAD)
from analysis.v5.run import blocks, cluster_fit, stats_row, TRAITS

ROOT = Path("work/submission")
HASHES = {
    "tavg":"5e567dcfe94379b94229492849ce91078b1c6e5210aaf435fba449fae6b95405",
    "tmax":"599d56a7745a2cd0a50dbe77e15fcb19d06a2c57e5c6922a18155388a58f4f20",
    "prec":"1090f578f6402672a2b6438e252156f6cf867a25002a8ad6ffcbfa19e516ddc6",
    "srad":"c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026",
}


def dates():
    path = ROOT/"private/original-metadata.aesgcm"
    raw = decrypt_checked(path.read_bytes(), os.environ["AZAMI_METADATA_KEY"],
                          ORIGINAL_SEALED_SHA, ORIGINAL_SHA, ORIGINAL_AAD)
    with ZipFile(io.BytesIO(raw)) as archive, archive.open("photo_metadata_merged.csv") as f:
        dates = pd.read_csv(f, usecols=["obs_id","observed_on"],dtype={"obs_id":str}).drop_duplicates()
    if dates.obs_id.duplicated().any():
        raise ValueError("Conflicting observation dates")
    return dates


def normals(frame, folder):
    lat, lon = frame.latitude.to_numpy(), frame.longitude.to_numpy()
    result = {}
    for axis, digest in HASHES.items():
        target = folder/(axis+".zip")
        if not target.exists():
            url = f"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_{axis}.zip"
            response = requests.get(url,timeout=(30,180))
            response.raise_for_status()
            target.write_bytes(response.content)
        if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
            raise ValueError(f"WorldClim {axis} checksum mismatch")
        values = []
        for month in range(1,13):
            name = f"wc2.1_10m_{axis}_{month:02d}.tif"
            with rasterio.open(f"zip://{target.resolve()}!{name}") as src:
                ri, ci = rasterio.transform.rowcol(src.transform,lon,lat)
                r = np.clip(ri,0,src.height-1)
                c = np.clip(ci,0,src.width-1)
                v = src.read(1,masked=True).astype("float32").filled(np.nan)[r,c]
                v[(r!=ri)|(c!=ci)] = np.nan
                values.append(v)
        result[axis] = np.column_stack(values)
    return result


def collect(frame, monthly):
    m = frame.month.fillna(1).astype(int).to_numpy()-1
    quarter = np.stack([np.roll(monthly["tavg"],-i,axis=1) for i in range(3)]).mean(axis=0)
    valid = np.isfinite(quarter).all(axis=1)
    start = np.nan_to_num(quarter,nan=-1e6).argmax(axis=1)
    ids = np.arange(len(frame))
    for window, indices in {
        "warm_quarter":(start[:,None]+np.arange(3))%12,
        "observation_month":m[:,None],
    }.items():
        for axis, values in monthly.items():
            selected = np.take_along_axis(values,indices,axis=1)
            for kind in ["mean","peak"]:
                v = selected.mean(1) if kind=="mean" else selected.max(1)
                frame[f"{window}_{kind}_{axis}"] = np.where(valid & frame.month.notna(),v,np.nan)
    frame["warm_quarter_bio05"] = np.where(valid & frame.month.notna(),monthly["tmax"].max(axis=1),np.nan)
    # All tied warmest quarters count as qualifying in the alignment audit.
    ties = np.isclose(quarter,quarter.max(1,keepdims=True),rtol=0,atol=1e-6)
    membership = np.logical_or.reduce([ties[ids,(m-i)%12] for i in range(3)])
    frame["in_warm"] = np.where(valid & frame.month.notna(),membership,np.nan)
    return frame


def evaluate(frame):
    summaries, coeffs, diagnostics = [], [], []
    frame["cell"] = blocks(frame,5)
    for trait, col in TRAITS.items():
        eligible = frame[col+"__eligible"].eq(True) & np.isfinite(frame[col])
        if trait=="angle":
            eligible &= frame[col].between(0,180)
        d = frame[eligible].copy()
        aligned = d.dropna(subset=["in_warm"])
        fractions = aligned.groupby(["accepted_key","cell"]).in_warm.mean().groupby(level=0).mean()
        summaries.append(dict(trait=trait,eligible=len(d),dated=d.month.notna().sum(),matched=len(aligned),
                              taxa=len(fractions),raw_share=aligned.in_warm.mean(),balanced_share=fractions.mean(),
                              taxon_p25=fractions.quantile(.25),taxon_median=fractions.median(),taxon_p75=fractions.quantile(.75)))
        for window in ["warm_quarter","observation_month"]:
            specs = {"mean":["mean_tavg","mean_prec","mean_srad"],
                     "peak":["peak_tmax","peak_prec","peak_srad"]}
            if window=="observation_month":
                specs["peak"]=["mean_tmax","mean_prec","mean_srad"]
            hottest = "bio05" if window=="warm_quarter" else "mean_tmax"
            specs["hot_with_mean_rain_sun"] = [hottest,"mean_prec","mean_srad"]
            specs["four"] = ["mean_tavg",hottest,"mean_prec","mean_srad"]
            columns = list(dict.fromkeys(sum(specs.values(),[])))
            fields = [window+"_"+v for v in columns]
            complete = d.dropna(subset=fields)
            cells = complete.groupby(["accepted_key","cell"])[[col,"latitude","longitude",*fields]].mean().reset_index()
            cells = cells[cells.accepted_key.map(complete.accepted_key.value_counts())>=5].copy()
            cells.rename(columns=dict(zip(fields,columns)),inplace=True)
            weights = 1/cells.accepted_key.map(cells.accepted_key.value_counts()).to_numpy()
            for key in [col,*columns]:
                v = cells[key]
                cells[key] = (v-np.average(v,weights=weights))/np.sqrt(np.average((v-np.average(v,weights=weights))**2,weights=weights))
            avg = cells.groupby("accepted_key")[columns].transform("mean")
            levels = {"between":avg,"within":cells[columns]-avg}
            for spec,axes in {**specs,"combined":columns}.items():
                for scale, values in levels.items():
                    cov = np.cov(values[axes].T,aweights=weights)
                    corr = cov/np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
                    vif = []
                    for j in range(len(axes)):
                        keep = np.arange(len(axes)) != j
                        r = corr[j,keep]
                        unique = 1-r @ np.linalg.pinv(corr[np.ix_(keep,keep)]) @ r
                        vif.append(1/unique if unique > 1e-10 else np.inf)
                    for i,a in enumerate(axes):
                        diagnostics.append(dict(trait=trait,window=window,spec=spec,scale=scale,a=a,b=a,stat="VIF",value=vif[i]))
                        for j in range(i):
                            diagnostics.append(dict(trait=trait,window=window,spec=spec,scale=scale,a=a,b=axes[j],stat="correlation",value=corr[i,j]))
                if spec=="combined":
                    continue
                design = sm.add_constant(np.column_stack([v[axes] for v in levels.values()]),has_constant="add")
                fit = cluster_fit(cells[col],design,weights,blocks(cells,20))
                for k,scale in enumerate(levels):
                    for j,axis in enumerate(axes):
                        coeffs.append(dict(trait=trait,window=window,spec=spec,axis=axis,scale=scale,
                                           taxa=cells.accepted_key.nunique(),cells=len(cells),**stats_row(fit,1+k*len(axes)+j)))
    return {"seasonal_alignment":summaries,"seasonal_slopes":coeffs,"seasonal_collinearity":diagnostics}


def main():
    data = pd.read_csv(ROOT/"private/decoded/observations.csv.gz",dtype={"obs_id":str},low_memory=False)
    if "observed_on" not in data:
        data = data.merge(dates(),on="obs_id",how="left",validate="one_to_one")
    data["month"] = pd.to_datetime(data.observed_on,errors="coerce",format="mixed").dt.month
    data = data[np.isfinite(data.latitude)&np.isfinite(data.longitude)].copy()
    folder = ROOT/"temp"
    outputs = evaluate(collect(data,normals(data,folder)))
    for name, rows in outputs.items():
        pd.DataFrame(rows).to_csv(ROOT/"v5"/(name+".csv"),index=False)

if __name__=="__main__":
    main()
