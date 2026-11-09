"""Native-referenced joint precipitation contrasts, cell-weighted by range."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from analysis.v5.novel_head_overlap import T


def run(native, introduced, core, out):
    cols = ["chelsa_bio12", "chelsa_bio15"]
    data = pd.concat([native.assign(inv=0), introduced.assign(inv=1)], ignore_index=True)
    data["outside"] = ~core[data.u, data.v]
    results, summaries, target = [], [], []
    for q in [.5, 1/3]:
        rain = native[cols[0]].quantile(q)
        season = native[cols[1]].quantile(1-q)
        for scope, d in [("whole", data.copy()), ("shared_background", data[data.common].copy())]:
            dry = d[cols[0]] <= rain
            seasonal = d[cols[1]] >= season
            d["group"] = 2*dry.astype(int) + seasonal.astype(int)
            groups = d.groupby(["inv","group"])
            for (inv,k),g in groups:
                summaries.append(dict(q=q,rain=rain,season=season,scope=scope,introduced=inv,group=k,
                    cells=len(g),blocks=g.block.nunique(),rare_mean=g.climate_rarity.mean(),
                    **{col:g[col].mean() for col in T}))
            for term, subset, variable in [
                ("introduced_minus_native_within_target", d[d.group==3], "inv"),
                ("outside_minus_inside_in_introduced_target", d[(d.group==3)&(d.inv==1)], "outside"),
            ]:
                counts = subset.groupby(variable).size()
                for col in T:
                    row = dict(scope=scope,q=q,term=term,trait=col,
                               first=counts.get(0,0) if variable=="inv" else int((~subset.outside).sum()),
                               second=counts.get(1,0) if variable=="inv" else int(subset.outside.sum()))
                    if min(row["first"],row["second"]) >= 5 and subset.block.nunique()>2:
                        fit = sm.OLS(subset[col],sm.add_constant(subset[variable].astype(float))).fit(
                            cov_type="cluster",cov_kwds={"groups":subset.block},use_t=True)
                        row.update(effect=fit.params.iloc[1],lo=fit.conf_int().iloc[1,0],
                                   hi=fit.conf_int().iloc[1,1],p=fit.pvalues.iloc[1])
                    target.append(row)
            # Four climates, independent group means in each range. A spatial block
            # shared by both ranges gets one cluster ID.
            category = d.inv*4 + d.group
            design = np.eye(8)[category]
            if np.linalg.matrix_rank(design) != 8:
                raise ValueError("Missing precipitation combination")
            for col in T:
                fit = sm.OLS(d[col], design).fit(cov_type="cluster",cov_kwds={"groups":d.block},use_t=True)
                for ref in [0,1,2]:
                    cn = np.eye(8)[3] - np.eye(8)[ref]
                    ci = np.eye(8)[7] - np.eye(8)[4+ref]
                    for name, c in [("native",cn),("introduced",ci),("difference",ci-cn)]:
                        test = fit.t_test(c)
                        lo,hi = test.conf_int().ravel()
                        results.append(dict(q=q,scope=scope,ref=ref,trait=col,term=name,
                            effect=test.effect.item(),lo=lo,hi=hi,p=test.pvalue.item()))
    t = pd.DataFrame(target)
    t["holm"] = t.groupby(["scope","q","term"]).p.transform(lambda v:multipletests(v,method="holm")[1] if v.notna().all() else np.nan)
    t.to_csv(out/"vulgare_target_comparison.csv",index=False)
    r = pd.DataFrame(results)
    r["holm"] = r.groupby(["scope","q","ref","term"]).p.transform(lambda v:multipletests(v,method="holm")[1])
    r.to_csv(out/"vulgare_rain_contrasts.csv",index=False)
    pd.DataFrame(summaries).to_csv(out/"vulgare_rain_groups.csv",index=False)
