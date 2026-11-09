"""Compare trait departures from native geographic variation with climate shifts."""
from itertools import permutations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).parent/"results"
OUT=ROOT/"regional_climate_shifts"
AXES=["chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15","chelsa_rsds_mean"]


def correlation(x,y,perms):
    rx=rankdata(x);rx-=rx.mean()
    ry=rankdata(y);ry-=ry.mean()
    rho=(rx@ry)/(np.linalg.norm(rx)*np.linalg.norm(ry))
    p=np.mean(abs(ry[perms]@rx)>=abs(rx@ry)-1e-9)
    return rho,p


def main():
    table=pd.read_csv(OUT/"native_baseline_regions.csv")
    old=pd.read_csv(ROOT/"invasion_design/native_baseline_regions.csv")
    keys=["species","axis","range","region_index"]
    repeat=old.merge(table,on=keys,validate="one_to_one")
    assert len(repeat)==len(old)
    np.testing.assert_allclose(repeat.dissimilarity_x,repeat.dissimilarity_y,atol=1e-10)
    correlations=[];descriptive=[];profile=[]
    for sp,rows in table.groupby("species"):
        native=rows.query("range=='native'").pivot(index="region_index",columns="axis",values="dissimilarity")
        target=rows.query("range=='introduced'")
        delta=target.pivot(index="region_index",columns="axis",values="signed_shift")
        dif=target.pivot(index="region_index",columns="axis",values="dissimilarity")
        perms=np.asarray(list(permutations(range(len(delta)))))
        for trait in ["angle","chroma"]:
            departure=dif[trait]-native[trait].quantile(.90)
            flag=departure>0
            for idx in delta.index:
                record=dict(species=sp,trait=trait,region=int(idx),exceeds=bool(flag[idx]),
                            excess_over_native_p90=departure[idx],trait_shift=delta.loc[idx,trait])
                record.update({axis:delta.loc[idx,axis] for axis in AXES})
                profile.append(record)
            for axis in AXES:
                for test,x,y in [("magnitude",dif[axis],departure),
                                 ("climate_direction",delta[axis],departure),
                                 ("signed",delta[axis],delta[trait])]:
                    rho,p=correlation(x,y,perms)
                    loo=[x.drop(k).corr(y.drop(k),method="spearman") for k in x.index]
                    correlations.append(dict(species=sp,trait=trait,climate=axis,test=test,rho=rho,
                                             exact_p=p,loo_min=min(loo),loo_max=max(loo)))
                if 0<flag.sum()<len(flag):
                    a=delta[axis].to_numpy()
                    idx=flag.to_numpy()
                    contrast=a[idx].mean()-a[~idx].mean()
                    perm_delta=a[perms][:,idx].mean(1)-a[perms][:,~idx].mean(1)
                    descriptive.append(dict(species=sp,trait=trait,climate=axis,n_above=idx.sum(),
                        n_other=(~idx).sum(),mean_climate_shift_above=a[idx].mean(),mean_climate_shift_other=a[~idx].mean(),
                        contrast=contrast,p=np.mean(abs(perm_delta)>=abs(contrast)-1e-9)))
    summary=pd.DataFrame(correlations)
    summary["holm"]=summary.groupby("test").exact_p.transform(lambda p:multipletests(p,method="holm")[1])
    summary.to_csv(OUT/"regional_tests.csv",index=False)
    desc=pd.DataFrame(descriptive)
    desc["holm"]=multipletests(desc.p,method="holm")[1]
    desc.to_csv(OUT/"exceedance_groups.csv",index=False)
    pd.DataFrame(profile).to_csv(OUT/"regional_profiles.csv",index=False)
if __name__=="__main__":
    main()
