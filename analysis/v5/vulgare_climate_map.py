"""Aggregate introduced C. vulgare cells into climate space without geographic identifiers."""
import numpy as np
import pandas as pd
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v5.novel_head_overlap import T


def export(cells, core, out):
    d = cells.copy()
    temp, _, rain, _, _ = VARIABLES
    d["rare"] = ~core[d.u, d.v]
    descriptors = [*VARIABLES, T[1], "climate_rarity"]
    profiles = []
    tiles = []
    for scope, g in [("all", d), ("common", d[d.common])]:
        for label, block in g.groupby("rare"):
            row = dict(scope=scope, rare=label, n=len(block))
            for col in descriptors:
                for name, q in [("q25", .25), ("median", .5), ("q75", .75)]:
                    row[col + "_" + name] = block[col].quantile(q)
            profiles.append(row)
        for view, x, y in [
            ("annual", np.floor(g[temp] / 3) * 3, np.floor(g[rain] / 300) * 300),
            ("pca", g.u // 4, g.v // 4),
        ]:
            for (a, b), block in g.groupby([x.rename("x"), y.rename("y")]):
                row = dict(scope=scope, view=view, x=a, y=b, n=len(block), rare_n=int(block.rare.sum()))
                if len(block) >= 5:
                    row.update(block[descriptors].median().to_dict())
                tiles.append(row)
    pd.DataFrame(tiles).to_csv(out / "vulgare_tiles.csv", index=False)
    pd.DataFrame(profiles).to_csv(out / "vulgare_profiles.csv", index=False)
