import numpy as np

from analysis.v3.hierarchical_ecology import (
    common_spatial_residualize,
    estimate_taxon_slopes,
    random_effects_reml,
    spherical_basis,
)


def test_common_spatial_residualization_recovers_within_taxon_signal():
    rng = np.random.default_rng(20260908)
    taxa = np.repeat(["a", "b", "c", "d"], 40)
    lat = np.concatenate([rng.normal(10+i*10, 1.0, 40) for i in range(4)])
    lon = np.concatenate([rng.normal(30+i*15, 1.0, 40) for i in range(4)])
    x = rng.normal(size=len(taxa))
    basis = spherical_basis(lat, lon)
    taxon_intercept = np.repeat([4.0, -3.0, 8.0, 1.0], 40)
    y = taxon_intercept + 1.8*x + 2.5*basis[:, 0] - 1.2*basis[:, 3] + rng.normal(0, .15, len(taxa))
    yr, xr, info = common_spatial_residualize(y, x[:, None], taxa, lat, lon)
    beta = float(np.linalg.lstsq(xr, yr, rcond=None)[0][0])
    assert abs(beta - 1.8) < .08
    assert info["taxa"] == 4
    assert info["spatial_basis_rank"] >= 1


def test_supported_taxon_slopes_and_reasons_are_explicit():
    rng = np.random.default_rng(8)
    taxa = np.repeat(["a", "b", "c"], [20, 20, 6])
    cells = np.array([f"a{i%5}" for i in range(20)] + [f"b{i%5}" for i in range(20)] + [f"c{i%2}" for i in range(6)])
    x1 = rng.normal(size=len(taxa)); x2 = rng.normal(size=len(taxa))
    y = np.where(taxa=="a", 1.0*x1-.4*x2, np.where(taxa=="b", .2*x1+.8*x2, .5*x1)) + rng.normal(0,.2,len(taxa))
    slopes, ledger = estimate_taxon_slopes(y, np.column_stack([x1,x2]), taxa, cells, ["p1","p2"])
    assert {s.taxon for s in slopes} == {"a", "b"}
    row = ledger.set_index("taxon").loc["c"]
    assert not bool(row["slope_estimable"])
    assert "below_minimum_observations" in row["reasons"]
    assert "below_minimum_cells" in row["reasons"]


def test_random_effects_reml_shrinks_noisy_taxon_slopes():
    result = random_effects_reml(
        estimates=[0.9, 1.1, 1.0, 3.5],
        variances=[0.01, 0.01, 0.02, 4.0],
        taxa=["a", "b", "c", "noisy"],
    )
    assert 0.8 < result.hypermean < 1.3
    noisy = result.shrinkage.set_index("taxon").loc["noisy"]
    assert abs(noisy["shrunken_slope"] - result.hypermean) < abs(noisy["raw_slope"] - result.hypermean)
    assert result.tau2 >= 0
