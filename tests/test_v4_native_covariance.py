import numpy as np
import pytest
from analysis.v4.native_covariance import kernels, solve_components, psd


def test_expected_moments_recover_full_covariances_unequal_n():
    a = np.array([[1, .8, .1, .1], [.8, 1, .1, .1],
                  [.1, .1, 1, .4], [.1, .1, .4, 1]])
    _, gram, noise, _ = kernels(a, [5, 20, 60, 8])
    truth = np.array([[[.4, -.1], [-.1, .3]], [[.2, .05], [.05, .5]]])
    r = np.array([[1, .3], [.3, 2]])
    expectation = np.einsum("lk,kij->lij", gram, truth)+noise[:, None, None]*r
    fitted, unconstrained, _ = solve_components(expectation-noise[:, None, None]*r, gram)
    np.testing.assert_allclose(fitted, truth, atol=1e-7)
    np.testing.assert_allclose(unconstrained, truth, atol=1e-7)


def test_proportional_kernels_stop():
    with pytest.raises(ValueError, match="not separable"):
        kernels(np.eye(5), np.ones(5)*10)


def test_psd_solver_boundary():
    gram = np.array([[3., 1.], [1., 2.]])
    s = np.array([[[-1., .2], [.2, .3]], [[1., .4], [.4, 2.]]])
    fit, _, _ = solve_components(s, gram)
    assert min(np.linalg.eigvalsh(fit).ravel()) > -1e-10
    # Fixed-point optimality for each PSD block.
    np.testing.assert_allclose(fit[0], psd((s[0]-fit[1])/3), atol=1e-7)
    np.testing.assert_allclose(fit[1], psd((s[1]-fit[0])/2), atol=1e-7)
