import pandas as pd
from analysis.v5.phenotype_temporal_support import temporal_groups,supported_strata
def test_temporal_support():
    rows=[]
    for yr in [2016,2017,2021,2022]:
        for user in ["a","b","c"]:
            rows.append(dict(obs_id=f"{yr}{user}",latitude=20.1,longitude=30.1,month=6,year=yr,user_id=user,novel=1))
    d=temporal_groups(pd.DataFrame(rows),.25)
    q=supported_strata(d,"middle","late",3)
    assert len(q)==1 and q.climate_group.iloc[0]=="novel"
    d=d.loc[d.year!=2022]
    assert supported_strata(d,"middle","late",3).empty
