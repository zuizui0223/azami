import numpy as np
from Bio import Phylo
from io import StringIO
from analysis.v4.native_phylogeny import statistics, evaluate, tree_covariance


def test_exact_crossproducts_match_dense_observation_reml():
    rng = np.random.default_rng(42)
    groups = np.repeat(np.arange(6), [3, 7, 4, 9, 2, 5])
    x = np.column_stack([np.ones(len(groups)), rng.normal(size=(len(groups), 2))])
    y = rng.normal(size=len(groups))
    a = .7*np.eye(6) + .3*np.ones((6, 6))
    s = statistics(x, y, groups)
    for ratios in ([0, 0], [0, .6], [.8, .2], [.8, 0]):
        obj, beta, cov, variance = evaluate(s, a, ratios)
        z = np.eye(6)[groups]
        v = np.eye(len(y)) + z @ (ratios[0]*a+ratios[1]*np.eye(6)) @ z.T
        info = x.T @ np.linalg.solve(v, x)
        expected = np.linalg.solve(info, x.T @ np.linalg.solve(v, y))
        residual = y-x@expected
        var = residual @ np.linalg.solve(v, residual)/(len(y)-x.shape[1])
        expected_obj = ((len(y)-x.shape[1])*np.log(var)
                        +np.linalg.slogdet(v)[1]+np.linalg.slogdet(info)[1])
        np.testing.assert_allclose(beta, expected, atol=1e-10)
        np.testing.assert_allclose(cov, var*np.linalg.inv(info), atol=1e-10)
        np.testing.assert_allclose(obj, expected_obj, atol=1e-10)


def test_tree_covariance_tip_order_and_shared_branch():
    tree = Phylo.read(StringIO("((A:1,B:1):2,C:3);"), "newick")
    a = tree_covariance(tree, ["C", "B", "A"])
    np.testing.assert_allclose(a, [[1, 0, 0], [0, 1, 2/3], [0, 2/3, 1]])
