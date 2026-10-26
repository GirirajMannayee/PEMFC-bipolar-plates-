"""Potentiodynamic polarisation / LPR analysis (ASTM G59, G102; Stern-Geary)."""
from __future__ import annotations
import numpy as np


def find_ecorr(E, i):
    """Ecorr = potential of the minimum |i| (simple, robust for clean scans)."""
    E, i = np.asarray(E, float), np.asarray(i, float)
    return float(E[np.argmin(np.abs(i))])


def tafel_cathodic(E, i, ecorr, window=(-0.25, -0.12)):
    """Cathodic Tafel extrapolation (manuscript Sec. 2.7).

    Fits log10|i| vs E over [Ecorr+window[0], Ecorr+window[1]] and extrapolates to
    Ecorr. Returns Icorr (same unit as i) and beta_c (V/decade, positive).
    """
    E, i = np.asarray(E, float), np.asarray(i, float)
    mask = (E >= ecorr + window[0]) & (E <= ecorr + window[1]) & (i < 0)
    if mask.sum() < 5:
        raise ValueError("Too few points in the cathodic Tafel window")
    slope, intercept = np.polyfit(E[mask], np.log10(np.abs(i[mask])), 1)
    icorr = 10 ** (slope * ecorr + intercept)
    return dict(icorr=icorr, beta_c=-1.0 / slope, n_points=int(mask.sum()))


def lpr_rp(E, i, ecorr, half_window=0.005):
    """Polarisation resistance from the linear region Ecorr +/- 5 mV (ohm cm2 if i in A/cm2)."""
    E, i = np.asarray(E, float), np.asarray(i, float)
    mask = np.abs(E - ecorr) <= half_window
    if mask.sum() < 3:
        raise ValueError("Too few points in the LPR window")
    slope, _ = np.polyfit(i[mask], E[mask], 1)
    return float(slope)


def stern_geary_B(beta_a, beta_c):
    """B = ba*bc / (2.303 (ba+bc)) (same unit as the betas)."""
    return beta_a * beta_c / (2.303 * (beta_a + beta_c))


def icorr_from_rp(Rp, B):
    """Eq. (4): Icorr = B / Rp."""
    return B / Rp


def implied_beta_a(B, beta_c):
    """Anodic slope implied by B and beta_c: solves B = ba bc /(2.303 (ba+bc))."""
    return 2.303 * B * beta_c / (beta_c - 2.303 * B)


def butler_volmer(E, ecorr, icorr, beta_a, beta_c):
    """Pure activation-controlled Butler-Volmer current (test / demo utility)."""
    return icorr * (10 ** ((E - ecorr) / beta_a) - 10 ** (-(E - ecorr) / beta_c))
