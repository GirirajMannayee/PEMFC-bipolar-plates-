"""Surrogate-assisted multi-objective optimisation (manuscript Sec. 2.10 / 3.11).

Pipeline
--------
1. Training set = 4 *measured* anchor compositions (S1-S4) + 40 Latin-hypercube
   design points (virtual).  The manuscript does not publish how the 40 virtual
   responses were generated; ``load_training_set`` therefore uses the authors'
   file ``data/optimisation/virtual_lhs_40.csv`` if present, otherwise a clearly
   labelled SYNTHETIC DEMO surface (see ``synthetic_truth``) so that the whole
   pipeline can be executed end-to-end.  Numbers obtained from the demo surface
   are NOT the manuscript's results.
2. RF / GBT / GP surrogates, repeated 5-fold CV (nested 3-fold tuning for GBT).
3. Pareto screening of 2^14 scrambled Sobol designs on GP means.
4. ParEGO-EI Bayesian loop vs random search (hypervolume, ref = 1.05).
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import qmc
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GridSearchCV, KFold, RepeatedKFold

from .constants import DESIGN_BOUNDS, PROCESSED, DATA

DESIGN_VARS = ["cu_wt", "cnt_wt", "rpm"]
RESPONSES = ["sigma", "ln_Rct", "ln_Icorr", "ln_ICR1"]
# True = larger is better.  Objectives are converted to minimisation internally.
MAXIMISE = {"sigma": True, "ln_Rct": True, "ln_Icorr": False, "ln_ICR1": False}
LOWER = np.array([DESIGN_BOUNDS[v][0] for v in DESIGN_VARS])
UPPER = np.array([DESIGN_BOUNDS[v][1] for v in DESIGN_VARS])
HV_REF = 1.05


# ----------------------------------------------------------------- data
def to_unit(X):
    return (np.asarray(X, float) - LOWER) / (UPPER - LOWER)


def from_unit(U):
    return LOWER + np.asarray(U, float) * (UPPER - LOWER)


def anchors() -> pd.DataFrame:
    """The four experimentally characterised compositions (S1-S4) from Tables 1, 6-9."""
    comp = pd.read_csv(PROCESSED / "table01_nominal_composition.csv").set_index("sample")
    t6 = pd.read_csv(PROCESSED / "table06_conductivity_hardness.csv").set_index("sample")
    t7 = pd.read_csv(PROCESSED / "table07_eis_25C.csv").set_index("sample")
    t8 = pd.read_csv(PROCESSED / "table08_polarisation_lpr.csv").set_index("sample")
    t9 = pd.read_csv(PROCESSED / "table09_icr_vs_pressure.csv").set_index("sample")
    df = pd.DataFrame({
        "cu_wt": comp.cu_added_wt_pct, "cnt_wt": comp.cnt_wt_pct, "rpm": 500.0,
        "sigma": t6.sigma_MS_per_m_mean,
        "ln_Rct": np.log(t7.Rct_kohm_cm2_mean),
        "ln_Icorr": np.log(t8.Icorr_uA_cm2_mean),
        "ln_ICR1": np.log(t9.icr_1p0MPa_mean),
    })
    df["origin"] = "measured"
    return df


def lhs_design(n=40, seed=42) -> np.ndarray:
    return from_unit(qmc.LatinHypercube(d=3, seed=seed).random(n))


def synthetic_truth(X, rng=None, noise=True):
    """SYNTHETIC DEMO response surface (NOT experimental data).

    Smooth saturating Cu and (dispersion-corrected) CNT effects with an
    interaction term, calibrated so that it reproduces the four measured anchor
    compositions exactly at 500 rpm, plus penalties for Cu > 6.5 wt.% and
    effective CNT > 1.1 wt.% (galvanic / agglomeration), mimicking the trends
    qualitatively described in Sec. 3.11.  Replace with the authors' dataset.
    """
    X = np.atleast_2d(np.asarray(X, float))
    cu, cnt, rpm = X[:, 0], X[:, 1], X[:, 2]
    a = anchors()
    g = lambda r: np.exp(-((r - 520.0) / 180.0) ** 2)   # dispersion efficiency vs rpm
    cnt_eff = cnt * g(rpm) / g(500.0)
    u = (1 - np.exp(-cu / 3.0)) / (1 - np.exp(-5.62 / 3.0))
    v = (1 - np.exp(-cnt_eff / 0.5)) / (1 - np.exp(-0.7 / 0.5))
    pen = 1.2 * np.maximum(0, cu - 5.8) ** 2 + 6.0 * np.maximum(0, cnt_eff - 1.0) ** 2
    scale = {"sigma": 0.3, "ln_Rct": 0.06, "ln_Icorr": 0.06, "ln_ICR1": 0.04}
    scale_pen = {"sigma": 1.0, "ln_Rct": 0.25, "ln_Icorr": 0.25, "ln_ICR1": 0.2}
    out = {}
    for r in RESPONSES:
        y1, y2, y3, y4 = (a.loc[s, r] for s in ["S1", "S2", "S3", "S4"])
        A, B, C = y2 - y1, y3 - y1, y4 - y2 - y3 + y1
        y = y1 + A * u + B * v + C * u * v
        y = y - pen * scale_pen[r] if MAXIMISE[r] else y + pen * scale_pen[r]
        if noise:
            y = y + (rng or np.random.default_rng(0)).normal(0, scale[r], size=len(y))
        out[r] = y
    return pd.DataFrame(out)


def load_training_set(seed=42):
    """Return (DataFrame of 44 rows, origin string).

    Looks for authors' ``data/optimisation/virtual_lhs_40.csv`` with columns
    cu_wt, cnt_wt, rpm, sigma, ln_Rct, ln_Icorr, ln_ICR1 (or Rct, Icorr, ICR1 in
    natural units - they are log-transformed automatically).
    """
    anc = anchors()
    f = DATA / "optimisation" / "virtual_lhs_40.csv"
    if f.exists():
        v = pd.read_csv(f)
        for raw, ln in [("Rct", "ln_Rct"), ("Icorr", "ln_Icorr"), ("ICR1", "ln_ICR1")]:
            if ln not in v.columns and raw in v.columns:
                v[ln] = np.log(v[raw])
        v = v[DESIGN_VARS + RESPONSES].copy()
        v["origin"] = "virtual (author-supplied)"
        origin = "AUTHOR_SUPPLIED"
    else:
        X = lhs_design(40, seed)
        v = pd.concat([pd.DataFrame(X, columns=DESIGN_VARS),
                       synthetic_truth(X, np.random.default_rng(seed))], axis=1)
        v["origin"] = "virtual (SYNTHETIC DEMO)"
        origin = "SYNTHETIC_DEMO"
    df = pd.concat([anc.reset_index(drop=True), v], ignore_index=True)
    return df, origin


# ----------------------------------------------------------------- models
def make_gp(seed=0):
    k = ConstantKernel(1.0, (1e-2, 1e2)) * Matern(length_scale=[1.0] * 3, length_scale_bounds=(0.05, 20.0), nu=2.5) \
        + WhiteKernel(1e-2, (1e-6, 1.0))
    return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=5, random_state=seed)


def make_rf(seed=0):
    return RandomForestRegressor(n_estimators=150, random_state=seed)


def make_gbt(seed=0):
    grid = {"max_depth": [2, 3], "learning_rate": [0.05, 0.1]}
    return GridSearchCV(GradientBoostingRegressor(n_estimators=200, random_state=seed), grid,
                        cv=KFold(3, shuffle=True, random_state=seed))


def cross_validate(U, Y: pd.DataFrame, seed=42, n_splits=5, n_repeats=2):
    """Pooled out-of-fold metrics + predictions for RF/GBT/GP. U = unit-scaled design."""
    rkf = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    makers = {"RF": make_rf, "GBT": make_gbt, "GP": make_gp}
    rows, preds = [], {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for r in Y.columns:
            y = Y[r].to_numpy()
            for name, mk in makers.items():
                P, T = [], []
                for tr, te in rkf.split(U):
                    m = mk(seed).fit(U[tr], y[tr])
                    P.append(m.predict(U[te])); T.append(y[te])
                P, T = np.concatenate(P), np.concatenate(T)
                sse, sst = np.sum((T - P) ** 2), np.sum((T - T.mean()) ** 2)
                rows.append(dict(response=r, model=name, R2=1 - sse / sst,
                                 MAE=mean_absolute_error(T, P), RMSE=float(np.sqrt(mean_squared_error(T, P)))))
                preds[(r, name)] = (T, P)
    return pd.DataFrame(rows), preds


def permutation_importance_oof(U, Y, seed=42, n_splits=5, n_repeats=20):
    """Out-of-fold permutation importance of the RF (mean, SD over folds)."""
    kf = KFold(n_splits, shuffle=True, random_state=seed)
    recs = []
    for r in Y.columns:
        y = Y[r].to_numpy()
        per_fold = []
        for tr, te in kf.split(U):
            m = make_rf(seed).fit(U[tr], y[tr])
            pi = permutation_importance(m, U[te], y[te], n_repeats=n_repeats, random_state=seed)
            per_fold.append(pi.importances_mean)
        per_fold = np.array(per_fold)
        for j, v in enumerate(DESIGN_VARS):
            recs.append(dict(response=r, variable=v, mean=per_fold[:, j].mean(), sd=per_fold[:, j].std(ddof=1)))
    return pd.DataFrame(recs)


def fit_gps(U, Y, seed=0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {r: make_gp(seed).fit(U, Y[r].to_numpy()) for r in Y.columns}


# ------------------------------------------------------- multi-objective
def to_min(F: np.ndarray) -> np.ndarray:
    """Convert [sigma, ln_Rct, ln_Icorr, ln_ICR1] to all-minimise form."""
    F = np.asarray(F, float).copy()
    for j, r in enumerate(RESPONSES):
        if MAXIMISE[r]:
            F[:, j] = -F[:, j]
    return F


def non_dominated(F: np.ndarray, chunk=512) -> np.ndarray:
    """Boolean mask of non-dominated rows (minimisation)."""
    n = len(F)
    mask = np.ones(n, bool)
    for s in range(0, n, chunk):
        blk = F[s:s + chunk]
        le = np.all(F[None, :, :] <= blk[:, None, :], axis=2)
        lt = np.any(F[None, :, :] < blk[:, None, :], axis=2)
        mask[s:s + chunk] = ~np.any(le & lt, axis=1)
    return mask


def hypervolume(P: np.ndarray, ref) -> float:
    """Exact dominated hypervolume (minimisation) by dimension-sweep (HSO)."""
    P = np.asarray(P, float)
    ref = np.asarray(ref, float)
    P = P[np.all(P < ref, axis=1)]
    if len(P) == 0:
        return 0.0
    return _hv(P, ref)


def _hv(P, ref):
    if P.shape[1] == 1:
        return float(ref[0] - P.min())
    P = P[np.argsort(P[:, -1])]
    vol = 0.0
    for i in range(len(P)):
        z_hi = P[i + 1, -1] if i + 1 < len(P) else ref[-1]
        depth = z_hi - P[i, -1]
        if depth > 0:
            sub = P[: i + 1, :-1]
            sub = sub[non_dominated(sub)]
            vol += depth * _hv(sub, ref[:-1])
    return float(vol)


def chebyshev(Fn: np.ndarray, lam, rho=0.05):
    """Eq. (6): g = max_i lam_i f_i + rho * sum_i lam_i f_i  (Fn normalised, minimise)."""
    w = Fn * np.asarray(lam)
    return w.max(axis=1) + rho * w.sum(axis=1)


PREFERENCES = {   # assumed weight vectors (sigma, ln Rct, ln Icorr, ln ICR) - not listed in the manuscript
    "Balanced": (0.25, 0.25, 0.25, 0.25),
    "Conductivity-priority": (0.55, 0.15, 0.15, 0.15),
    "Corrosion-priority": (0.10, 0.30, 0.50, 0.10),
}


def gp_predict(gps, U):
    mu = np.column_stack([gps[r].predict(U) for r in RESPONSES])
    sd = np.column_stack([gps[r].predict(U, return_std=True)[1] for r in RESPONSES])
    return mu, sd


def expected_improvement(mu, sd, best):
    """Eq. (7) for minimisation."""
    sd = np.maximum(sd, 1e-12)
    z = (best - mu) / sd
    return (best - mu) * stats.norm.cdf(z) + sd * stats.norm.pdf(z)


def parego(oracle, norm_lo, norm_hi, n_init=12, n_iter=20, seed=42, pool=2048, rho=0.05):
    """ParEGO-EI loop on the unit cube. ``oracle(U)`` -> raw responses (n x 4).

    Returns hypervolume traces for ParEGO and random search sharing the same
    initial design (normalised objectives, ref = 1.05).
    """
    rng = np.random.default_rng(seed)
    U0 = qmc.LatinHypercube(d=3, seed=seed).random(n_init)

    def nf(U):
        return (to_min(oracle(U)) - norm_lo) / (norm_hi - norm_lo)

    def trace(Fn_all):
        return [hypervolume(Fn_all[: k + 1][non_dominated(Fn_all[: k + 1])], [HV_REF] * 4) for k in range(len(Fn_all))]

    # Bayesian
    U, F = U0.copy(), nf(U0)
    for it in range(n_iter):
        lam = rng.dirichlet(np.ones(4))
        g = chebyshev(F, lam, rho)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gp = make_gp(seed + it).fit(U, g)
        cand = qmc.Sobol(d=3, scramble=True, seed=seed * 1000 + it).random(pool)
        mu, sd = gp.predict(cand, return_std=True)
        ei = expected_improvement(mu, sd, g.min())
        xn = cand[np.argmax(ei)][None, :]
        U, F = np.vstack([U, xn]), np.vstack([F, nf(xn)])
    hv_bo = trace(F)
    # Random
    Ur = np.vstack([U0, rng.random((n_iter, 3))])
    hv_rs = trace(nf(Ur))
    return np.array(hv_bo), np.array(hv_rs)
