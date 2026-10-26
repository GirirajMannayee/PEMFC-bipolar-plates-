"""Interfacial contact resistance analysis (manuscript Sec. 2.8)."""
from __future__ import annotations
import numpy as np


def quad_loglog_fit(P, icr):
    """ln ICR = a0 + a1 ln P + a2 (ln P)^2 -> coefficient vector [a2, a1, a0]."""
    return np.polyfit(np.log(P), np.log(icr), 2)


def quad_loglog_predict(coef, P):
    return np.exp(np.polyval(coef, np.log(P)))


def secant_exponent(P, icr):
    """b in ICR = k P^-b from the end points of the pressure range (secant)."""
    P, icr = np.asarray(P, float), np.asarray(icr, float)
    return float(np.log(icr[0] / icr[-1]) / np.log(P[-1] / P[0]))


def local_exponents(P, icr):
    P, icr = np.asarray(P, float), np.asarray(icr, float)
    return np.log(icr[:-1] / icr[1:]) / np.log(P[1:] / P[:-1])


def two_interface_basis(icr_per_interface):
    """Two-surface literature convention: double the per-interface value."""
    return 2.0 * np.asarray(icr_per_interface)
