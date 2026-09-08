import numpy as np

import analysis.v3.joint_partial_pooling as pooling


def _fixture(seed=2901):
    rng = np.random.default_rng(seed)
    groups = np.repeat(['a','b','c','d','e','f'], 24)
    x = rng.normal(size=(len(groups), 3))
    nuisance = np.column_stack([
        0.3 * x[:, 0] + rng.normal(size=len(groups)),
        rng.normal(size=len(groups)),
    ])
    slopes = rng.normal([0.6, -0.2, 0.4], [0.25, 0.2, 0.3], size=(6, 3))
    lookup = {label: i for i, label in enumerate(sorted(set(groups)))}
    y = np.array([x[i] @ slopes[lookup[groups[i]]] for i in range(len(groups))])
    y += nuisance @ np.array([0.35, -0.15]) + rng.normal(0, 0.6, len(groups))
    return y, x, groups, nuisance


def test_abnormal_lbfgs_flag_requires_independent_score_root_verification(monkeypatch):
    original = pooling.minimize

    def abnormal_after_candidate(*args, **kwargs):
        result = original(*args, **kwargs)
        result.success = False
        result.status = 2
        result.message = 'ABNORMAL: synthetic regression test'
        return result

    monkeypatch.setattr(pooling, 'minimize', abnormal_after_candidate)
    result = pooling.fit_joint_partial_pooling(*_fixture())
    verified = [attempt for attempt in result.optimizer_attempts if attempt['accepted']]
    assert verified
    assert all(not attempt['success'] for attempt in verified)
    assert all(attempt['termination_verified_by'] == 'score_root_polish' for attempt in verified)
    assert all(attempt['stationarity_polishing']['accepted'] for attempt in verified)
    assert all(attempt['stationarity_polishing']['minimum_free_curvature'] > 0 for attempt in verified)
    assert np.isfinite(result.hypermean).all()
