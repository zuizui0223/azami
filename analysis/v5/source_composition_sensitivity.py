"""How sensitive is an introduced pooled mean to the choice of native source?"""
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

ROOT=Path(__file__).parent/"results"
OUT=ROOT/"adjusted_exposures"
TRAITS=["angle","chroma","outline"]


def main():
    source=pd.read_csv(ROOT/"regional_climate_shifts/native_baseline_regions.csv")
    records=[];geometry=[]
    for sp,rows in source[source.axis.isin(TRAITS)].groupby("species"):
        means=rows.pivot(index=["range","region_index"],columns="axis",values="signed_shift")[TRAITS]
        native=means.loc["native"].to_numpy()
        n=len(native)
        # Published native shifts use a leave-one-native-out reference.
        native *= (n-1)/n
        np.testing.assert_allclose(native.mean(axis=0),0,atol=1e-7)
        target=means.loc["introduced"].mean(axis=0).to_numpy()
        for axes in [[0,1],[0,1,2]]:
            hull=ConvexHull(native[:,axes])
            inside=np.all(hull.equations[:,:-1]@target[axes]+hull.equations[:,-1]<=1e-9)
            geometry.append(dict(species=sp,traits="+".join(TRAITS[k] for k in axes),inside_native_mean_hull=inside))
        for k in [1,2,3]:
            subset=np.array(list(combinations(range(n),k)))
            scenarios=native[subset].mean(axis=1)
            diff=target-scenarios
            for j,trait in enumerate(TRAITS):
                changes=diff[:,j]
                records.append(dict(species=sp,trait=trait,n_native_regions=n,
                    n_intro_regions=len(means.loc["introduced"]),sources=k,scenarios=len(subset),
                    introduced_minus_full_native=target[j],difference_min=changes.min(),
                    difference_p05=np.quantile(changes,.05),difference_p95=np.quantile(changes,.95),
                    difference_max=changes.max(),positive_count=(changes>0).sum(),
                    negative_count=(changes<0).sum()))
    pd.DataFrame(records).to_csv(OUT/"source_composition.csv",index=False)
    pd.DataFrame(geometry).to_csv(OUT/"source_composition_joint.csv",index=False)
if __name__=="__main__":main()
