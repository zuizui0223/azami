import numpy as np
import pytest

from analysis.v3.hierarchical_ecology import (
    common_spatial_residualize,
    estimate_taxon_slopes,
    random_effects_reml,
    spherical_basis,
    joint_spatial_taxon_slopes,
    joint_nuisance_taxon_slopes,
    morans_i,
)


@pytest.mark.parametrize('nuisance_count',[0,4,13])
def test_shared_calendar_imaging_nuisance_matches_full_joint_interaction_reference(nuisance_count):
    rng=np.random.default_rng(81037)
    groups=np.repeat(['a','b','c'],45)
    n=len(groups)
    x=rng.normal(size=(n,3))
    nuisance=rng.normal(size=(n,nuisance_count))
    if nuisance_count:
        nuisance[:,0]+=.6*x[:,0]
        # Retain redundant and structurally constant nuisance terms without
        # making extra environmental slopes or silently changing membership.
        nuisance=np.column_stack([nuisance,nuisance[:,0],np.ones(n)])
    y=rng.normal(size=n)+x[:,0]*(1+np.repeat([-.8,.2,1.3],45))
    if nuisance_count: y+=nuisance[:,:nuisance_count]@rng.normal(size=nuisance_count)
    fit=joint_nuisance_taxon_slopes(y,x,groups,nuisance)
    indicators=np.column_stack([groups==t for t in fit.taxa])
    interactions=np.column_stack([x*(groups==t)[:,None] for t in fit.taxa])
    design=np.column_stack([indicators,nuisance,interactions])
    inverse=np.linalg.pinv(design)
    beta=inverse@y
    residual=y-design@beta
    leverage=np.einsum('ij,ji->i',design,inverse)
    covariance=(inverse*(residual/(1-leverage))**2)@inverse.T
    np.testing.assert_allclose(fit.slopes.ravel(),beta[-9:],atol=1e-10)
    np.testing.assert_allclose(fit.covariance,covariance[-9:,-9:],atol=1e-10)
    np.testing.assert_allclose(fit.residuals,residual,atol=1e-10)
    np.testing.assert_allclose(fit.leverage,leverage,atol=1e-10)
    assert fit.residual_df==n-np.linalg.matrix_rank(design)


def test_shared_nuisance_alias_or_invalid_rows_are_not_silently_repaired():
    rng=np.random.default_rng(431)
    n=90
    x=rng.normal(size=(n,2))
    y=rng.normal(size=n)
    groups=np.repeat(['a','b','c'],30)
    with pytest.raises(ValueError,match='aliased'):
        joint_nuisance_taxon_slopes(y,x,groups,x[:,:1])
    bad=np.ones((n,1)); bad[0,0]=np.nan
    with pytest.raises(ValueError,match='finite and aligned'):
        joint_nuisance_taxon_slopes(y,x,groups,bad)


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


def heterogeneous_fixture(noise=0.0):
    rng = np.random.default_rng(809)
    taxa = np.repeat(["a", "b"], 200)
    lat = rng.uniform(-70, 70, len(taxa))
    lon = rng.uniform(-180, 180, len(taxa))
    b = spherical_basis(lat, lon)
    x = 3 * b[:, 0] + rng.normal(0, .3, len(taxa))
    y = np.where(taxa == "a", 2.0, -1.0) * x + 1.3 * b[:, 2] + np.where(taxa == "a", 4., -3.)
    y += rng.normal(0, noise, len(y))
    return y, x[:, None], taxa, lat, lon


def test_joint_spatial_fit_recovers_heterogeneous_slopes_without_noise():
    y, x, g, lat, lon = heterogeneous_fixture()
    fit = joint_spatial_taxon_slopes(y, x, g, lat, lon)
    np.testing.assert_allclose(fit.slopes[:, 0], [2, -1], atol=1e-11)
    np.testing.assert_allclose(fit.residuals, 0, atol=1e-11)
    # Regression evidence for the rejected split-after-common-FWL shortcut:
    # the original commit produced 1.979819 and -1.265935, not 2 and -1.
    yr, xr, _ = common_spatial_residualize(y, x, g, lat, lon)
    rejected = [np.linalg.lstsq(xr[g == t], yr[g == t], rcond=None)[0][0] for t in fit.taxa]
    assert max(abs(np.asarray(rejected) - [2, -1])) > .25


@pytest.mark.parametrize("predictor_count", [1, 3])
def test_block_solver_and_full_hc3_covariance_match_dense_joint_reference(predictor_count):
    y, x, g, lat, lon = heterogeneous_fixture(noise=.6)
    if predictor_count > 1:
        extra = np.random.default_rng(8).normal(size=(len(y), predictor_count - 1))
        x = np.column_stack([x, extra])
        y += extra @ np.array([.7, -.3])
    fit = joint_spatial_taxon_slopes(y, x, g, lat, lon)
    indicators = np.column_stack([g == t for t in fit.taxa])
    interactions = np.column_stack([x * (g == t)[:, None] for t in fit.taxa])
    design = np.column_stack([indicators, spherical_basis(lat, lon), interactions])
    inverse = np.linalg.pinv(design)
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - design @ beta
    leverage = np.einsum("ij,ji->i", design, inverse)
    covariance = (inverse * (residual / (1 - leverage))**2) @ inverse.T
    size = interactions.shape[1]
    np.testing.assert_allclose(fit.slopes.ravel(), beta[-size:], atol=1e-10)
    np.testing.assert_allclose(fit.covariance, covariance[-size:, -size:], atol=1e-10)
    np.testing.assert_allclose(fit.leverage, leverage, atol=1e-11)
    np.testing.assert_allclose(fit.residuals, residual, atol=1e-10)
    assert fit.residual_df == len(y) - np.linalg.matrix_rank(design)
    assert abs(fit.covariance[0, -1]) > 1e-8  # shared nuisance induces cross-taxon dependence


def test_joint_solver_rejects_slope_spatial_alias_and_missing_taxa():
    y, x, g, lat, lon = heterogeneous_fixture()
    with pytest.raises(ValueError, match="aliased"):
        joint_spatial_taxon_slopes(y, spherical_basis(lat, lon)[:, :1], g, lat, lon)
    missing = g.astype(object)
    missing[0] = None
    with pytest.raises(ValueError, match="nonmissing"):
        joint_spatial_taxon_slopes(y, x, missing, lat, lon)


@pytest.mark.parametrize("constant", [0.1, 1e12])
def test_joint_solver_rejects_exact_constant_predictors(constant):
    y, x, g, lat, lon = heterogeneous_fixture()
    with pytest.raises(ValueError, match="not estimable"):
        joint_spatial_taxon_slopes(y, np.full_like(x, constant), g, lat, lon)


def test_centering_costs_an_intercept_degree_of_freedom():
    rng = np.random.default_rng(11)
    y = rng.normal(size=10)
    x = rng.normal(size=(10, 1))
    slopes, _ = estimate_taxon_slopes(y, x, np.repeat("a", 10), np.arange(10), ["x"])
    assert slopes[0].residual_df == 8


def test_moran_excludes_self_when_coordinates_coincide():
    # k=n-1 must give the complete graph without diagonal self-edges, even
    # when every point has identical coordinates (query order is then tied).
    n = 10
    assert morans_i(np.arange(n), np.zeros(n), np.zeros(n), k=n-1) == pytest.approx(-1/(n-1))


def test_moran_rejects_misaligned_coordinates():
    with pytest.raises(ValueError, match="length"):
        morans_i(np.arange(10), np.zeros(11), np.zeros(11))


@pytest.mark.parametrize("vary_longitude", [False, True])
def test_constant_or_redundant_spatial_terms_do_not_create_extra_rank(vary_longitude):
    rng = np.random.default_rng(915)
    n = 60
    taxa = np.repeat(["a", "b"], 30)
    lat = np.full(n, 37.2)
    lon = rng.uniform(-160, 160, n) if vary_longitude else np.full(n, 135.7)
    x = rng.normal(size=(n, 1))
    y = x[:, 0] + rng.normal(size=n)
    fit = joint_spatial_taxon_slopes(y, x, taxa, lat, lon)
    indicators = np.column_stack([taxa == t for t in fit.taxa])
    design = np.column_stack([indicators, spherical_basis(lat, lon),
                              *[x * (taxa == t)[:, None] for t in fit.taxa]])
    inverse = np.linalg.pinv(design)
    residual = y - design @ (inverse @ y)
    leverage = np.einsum("ij,ji->i", design, inverse)
    covariance = (inverse * (residual / (1-leverage))**2) @ inverse.T
    assert fit.residual_df == n - np.linalg.matrix_rank(design)
    np.testing.assert_allclose(fit.covariance, covariance[-2:, -2:], atol=1e-10)
    np.testing.assert_allclose(fit.leverage, leverage, atol=1e-11)
    if not vary_longitude:
        assert fit.spatial_rank == 0
