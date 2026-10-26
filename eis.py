"""EIS tools: Rs(Rct||CPE) model, Brug effective capacitance, fitting, Lin-KK test, Arrhenius."""
from __future__ import annotations
import numpy as np
from scipy.optimize import least_squares
from .constants import R_GAS


def z_randles_cpe(freq, Rs, Rct, Y0, n):
    """Eq. (1): Z = Rs + Rct / (1 + Rct*Y0*(j w)^n)."""
    jw = 1j * 2 * np.pi * np.asarray(freq, float)
    return Rs + Rct / (1.0 + Rct * Y0 * jw ** n)


def brug_ceff(Rs, Rct, Y0, n):
    """Eq. (2): Ceff = Y0^(1/n) * (1/Rs + 1/Rct)^((n-1)/n)."""
    return Y0 ** (1.0 / n) * (1.0 / Rs + 1.0 / Rct) ** ((n - 1.0) / n)


def fit_randles_cpe(freq, z, p0=None):
    """Complex NLLS fit with 1/|Z| weighting (manuscript Table 7 footnote).

    Parameters are fitted in log space (positivity). Returns dict with
    Rs, Rct, Y0, n, chi2 (weighted, per d.o.f.), Ceff, standard errors.
    """
    freq = np.asarray(freq, float)
    z = np.asarray(z, complex)
    w = 1.0 / np.abs(z)
    if p0 is None:
        Rs0 = float(np.real(z[np.argmax(freq)]))
        Rct0 = float(np.real(z[np.argmin(freq)]) - Rs0)
        p0 = [Rs0, max(Rct0, 1.0), 1e-5, 0.9]

    def unpack(p):
        return np.exp(p[0]), np.exp(p[1]), np.exp(p[2]), p[3]

    def resid(p):
        zm = z_randles_cpe(freq, *unpack(p))
        d = (zm - z) * w
        return np.concatenate([d.real, d.imag])

    x0 = [np.log(p0[0]), np.log(p0[1]), np.log(p0[2]), p0[3]]
    sol = least_squares(resid, x0, bounds=([-10, -5, -20, 0.3], [10, 20, 0, 1.0]), x_scale="jac")
    Rs, Rct, Y0, n = unpack(sol.x)
    dof = max(1, 2 * len(freq) - 4)
    chi2 = float(np.sum(sol.fun ** 2) / dof)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * chi2
        se = np.sqrt(np.diag(cov))
        se_phys = dict(Rs=Rs * se[0], Rct=Rct * se[1], Y0=Y0 * se[2], n=se[3])
    except np.linalg.LinAlgError:
        se_phys = {}
    return dict(Rs=Rs, Rct=Rct, Y0=Y0, n=n, chi2=chi2, Ceff=brug_ceff(Rs, Rct, Y0, n), se=se_phys)


def lin_kk(freq, z, c_threshold=0.85, m_max=None):
    """Linear Kramers-Kronig validity test (Boukamp 1995; Schoenleber 2014).

    Fits Rs + sum_k R_k/(1+j w tau_k) + j w L with M voigt elements (tau log-spaced over the
    measured range) for M = 1..m_max.  Among the fits that are not over-fitted -
    c = 1 - |sum(neg R_k)| / sum(pos R_k) >= c_threshold - the one with the smallest
    RMS relative residual is returned.  Residuals are relative to |Z|.
    """
    freq = np.asarray(freq, float)
    z = np.asarray(z, complex)
    w = 2 * np.pi * freq
    m_max = m_max or min(len(freq), 60)
    tau_min, tau_max = 1.0 / w.max(), 1.0 / w.min()
    wt = 1.0 / np.abs(z)
    cands = []
    for M in range(1, m_max + 1):
        tau = np.logspace(np.log10(tau_min), np.log10(tau_max), M) if M > 1 else np.array([np.sqrt(tau_min * tau_max)])
        A_re = np.column_stack([np.ones_like(w)] + [1 / (1 + (w * t) ** 2) for t in tau] + [np.zeros_like(w)])
        A_im = np.column_stack([np.zeros_like(w)] + [-(w * t) / (1 + (w * t) ** 2) for t in tau] + [w])
        A = np.vstack([A_re * wt[:, None], A_im * wt[:, None]])
        b = np.concatenate([z.real * wt, z.imag * wt])
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        Rk = x[1:-1]
        pos, neg = Rk[Rk > 0].sum(), Rk[Rk < 0].sum()
        c = 1 - abs(neg) / pos if pos > 0 else 0.0
        zfit = (A_re @ x) + 1j * (A_im @ x)
        res = (z - zfit) / np.abs(z)
        rms = float(np.sqrt(np.mean(np.concatenate([res.real, res.imag]) ** 2)))
        cands.append(dict(M=M, c=c, rms=rms, res_re=res.real, res_im=res.imag, zfit=zfit))
    ok = [d for d in cands if d["c"] >= c_threshold] or cands[:1]
    best = min(ok, key=lambda d: d["rms"])
    best["max_abs_residual_pct"] = 100 * float(np.max(np.abs(np.concatenate([best["res_re"], best["res_im"]]))))
    best["passes_1pct"] = best["max_abs_residual_pct"] <= 1.0
    return best


def arrhenius_fit(T_C, Rct):
    """Eq. (3): ln(1/Rct) = ln A - Ea/(R T). Returns Ea (kJ/mol), lnA, R2."""
    T = np.asarray(T_C, float) + 273.15
    y = np.log(1.0 / np.asarray(Rct, float))
    x = 1.0 / T
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    r2 = 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    return dict(Ea_kJ_mol=-slope * R_GAS / 1000.0, lnA=intercept, R2=r2)


def ea_two_point(T1_C, R1, T2_C, R2):
    """Apparent Ea (kJ/mol) from two Rct values."""
    T1, T2 = T1_C + 273.15, T2_C + 273.15
    return R_GAS * np.log(R1 / R2) / (1 / T1 - 1 / T2) / 1000.0
