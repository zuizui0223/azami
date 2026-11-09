"""Separate resemblance in marginal trait differences from their joint arrangement."""
from pathlib import Path

import numpy as np
import pandas as pd
from analysis.v4 import s1_h3_climate_analogue_alignment as analogue
from analysis.v4.s1_climate_trait_hypervolume import load_source, geo_aggregate, standardise_cells

ROOT = Path("work/bio18_geo")
# Laplace integral of sqrt(sum of squares) under independent empirical marginals.
nodes, weights = np.polynomial.legendre.leggauss(160)
times = np.exp(36 * nodes)
weights = 36 * weights * np.exp(-18 * nodes) / (2 * np.sqrt(np.pi))


def expected_distances(values, check=False):
    transforms = np.exp(-values[:, :, None]**2 * times).mean(axis=0)
    product = np.cumprod(transforms, axis=0)
    expected = (1 - product[1:]) @ weights
    observed = np.sqrt(np.cumsum(values**2, axis=1)[:, 1:]).mean(axis=0)
    if check:
        for dimension in (2, 3):
            grid = np.ix_(*(values[:, j]**2 for j in range(dimension)))
            exact = np.sqrt(sum(grid)).mean()
            np.testing.assert_allclose(expected[dimension-2], exact, atol=5e-6)
    return np.column_stack([observed, expected])


def main():
    source = ROOT / "private/decoded"
    df, _ = load_source(source / "observations.csv.gz", source / "measurement_receipt.json")
    rows = []
    for degrees in (1, 5):
        cells, _ = geo_aggregate(df, degrees)
        z, _ = standardise_cells(cells)
        taxa = analogue.eligible_taxa(z)
        rng = np.random.default_rng(20261109)
        draws, leave_one = [], []
        for rep in range(24):
            x, y, sector = analogue._sample_taxa(taxa, rng)
            pairs, support = analogue.matched_pairs(x, y, sector)
            if support["status"] != "qualified":
                continue
            margins = np.abs(y[pairs["i"], :, None, :] - y[pairs["j"], None, :, :])
            scores = []
            for k, (d, mask) in enumerate(zip(margins, pairs["mask"])):
                check = rep == 0 and k == 0
                scores.append(expected_distances(d.reshape(-1, 3), check) - expected_distances(d[mask], check))
            scores = np.stack(scores)
            draws.append(scores.mean(axis=0))
            excess = scores[:, :, 0] - scores[:, :, 1]
            leave_one.append([excess[(pairs["i"] != t) & (pairs["j"] != t)].mean(axis=0) for t in range(len(taxa))])
        obs, recombined = np.mean(draws, axis=0).T
        baseline = pd.read_csv(ROOT / "public/v5_trait_convergence.csv")
        original = baseline[(baseline.degrees == degrees) & (baseline.metric == "joint")].observed.item()
        np.testing.assert_allclose(obs[1], original)
        loo = np.mean(leave_one, axis=0)
        for j in range(2):
            rows.append(dict(degrees=degrees, dimensions=j+2, taxa=len(taxa), draws=len(draws),
                             observed=obs[j], independent_marginals=recombined[j], excess=obs[j]-recombined[j],
                             leave_one_taxon_min=loo[:, j].min(), leave_one_taxon_max=loo[:, j].max()))
    pd.DataFrame(rows).to_csv(ROOT / "public/v5_configuration_convergence.csv", index=False)


if __name__ == "__main__":
    main()
