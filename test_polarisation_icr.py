import numpy as np
from adc12pem import icr, polarisation as pol


def test_tafel_and_lpr_recover_butler_volmer():
    icorr, ba, bc, ec = 2e-6, 0.060, 0.120, -0.65
    E = np.linspace(ec - 0.25, ec + 0.3, 2000)
    i = pol.butler_volmer(E, ec, icorr, ba, bc)
    ec_hat = pol.find_ecorr(E, i)
    assert abs(ec_hat - ec) < 1e-3
    tf = pol.tafel_cathodic(E, i, ec_hat)
    assert abs(tf["icorr"] / icorr - 1) < 0.05 and abs(tf["beta_c"] / bc - 1) < 0.03
    rp = pol.lpr_rp(E, i, ec_hat)
    B = pol.stern_geary_B(ba, bc)
    assert abs(rp - B / icorr) / (B / icorr) < 0.03


def test_implied_beta_a_roundtrip():
    B = pol.stern_geary_B(0.022, 0.120)
    assert abs(pol.implied_beta_a(B, 0.120) - 0.022) < 1e-9


def test_icr_power_law_recovered():
    P = np.array([0.5, 1.0, 1.5, 2.0]); y = 4.0 * P ** -0.6
    assert abs(icr.secant_exponent(P, y) - 0.6) < 1e-9
    assert abs(icr.quad_loglog_predict(icr.quad_loglog_fit(P, y), 1.4) - 4.0 * 1.4 ** -0.6) < 1e-9
