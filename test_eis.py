import numpy as np
from adc12pem import eis


def synthetic(Rs=4.8, Rct=1800.0, Y0=2.1e-5, n=0.88, noise=0.002, seed=1):
    f = np.logspace(5, -2, 71)
    z = eis.z_randles_cpe(f, Rs, Rct, Y0, n)
    rng = np.random.default_rng(seed)
    return f, z * (1 + noise * (rng.standard_normal(len(f)) + 1j * rng.standard_normal(len(f))))


def test_fit_recovers_parameters():
    f, z = synthetic()
    r = eis.fit_randles_cpe(f, z)
    assert abs(r["Rct"] / 1800 - 1) < 0.02 and abs(r["Rs"] / 4.8 - 1) < 0.05 and abs(r["n"] - 0.88) < 0.01
    assert abs(r["Ceff"] * 1e6 - eis.brug_ceff(4.8, 1800, 2.1e-5, 0.88) * 1e6) < 0.3


def test_kk_passes_valid_and_flags_invalid():
    f, z = synthetic(noise=0.0005)
    assert eis.lin_kk(f, z)["passes_1pct"]
    z_bad = z.copy(); z_bad[30:40] *= 1.15          # non-causal distortion
    assert not eis.lin_kk(f, z_bad)["passes_1pct"]


def test_arrhenius_roundtrip():
    T = np.array([25, 40, 60, 80.0]); Ea = 30e3
    R = 2.0 * np.exp(Ea / 8.314462618 * (1 / (T + 273.15) - 1 / 298.15))
    assert abs(eis.arrhenius_fit(T, R)["Ea_kJ_mol"] - 30.0) < 1e-6
    assert abs(eis.ea_two_point(25, R[0], 80, R[-1]) - 30.0) < 1e-6
