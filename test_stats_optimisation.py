import numpy as np
from adc12pem import optimisation as O, stats as S


def test_holm_matches_known_values():
    assert np.allclose(S.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])


def test_factorial_additive_data_has_no_interaction():
    # cells exactly additive: 10, 12, 13, 15
    r = S.factorial_2x2([10, 12, 13, 15], [0.5] * 4, 3)
    assert abs(r["interaction"]["estimate"]) < 1e-12 and r["interaction"]["p"] > 0.99


def test_hypervolume_and_nondominated():
    assert abs(O.hypervolume(np.array([[0.5, 0.0], [0.0, 0.5]]), [1, 1]) - 0.75) < 1e-12
    assert abs(O.hypervolume(np.zeros((1, 4)), [1.05] * 4) - 1.05 ** 4) < 1e-12
    F = np.array([[0, 1], [1, 0], [1, 1], [0.5, 0.5]])
    assert list(O.non_dominated(F)) == [True, True, False, True]


def test_hypervolume_matches_monte_carlo():
    rng = np.random.default_rng(3); F = rng.random((25, 4)); P = F[O.non_dominated(F)]
    mc = rng.random((200000, 4)) * 1.05; dom = np.zeros(len(mc), bool)
    for p in P:
        dom |= np.all(mc >= p, axis=1)
    assert abs(O.hypervolume(P, [1.05] * 4) - dom.mean() * 1.05 ** 4) < 0.01


def test_synthetic_surface_reproduces_anchors_and_training_set_shape():
    a = O.anchors()
    t = O.synthetic_truth(a[O.DESIGN_VARS].values, noise=False)
    assert np.allclose(t.values, a[O.RESPONSES].values, atol=1e-9)
    df, origin = O.load_training_set()
    assert len(df) == 44 and df.origin.eq("measured").sum() == 4
    X = df[O.DESIGN_VARS].values
    assert (X >= O.LOWER - 1e-9).all() and (X <= O.UPPER + 1e-9).all()
