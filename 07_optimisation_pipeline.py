#!/usr/bin/env python3
"""AI-guided multi-objective optimisation (manuscript Sec. 2.10 / 3.11): Figs 21-24, Tables 14-15.

IMPORTANT - dataset provenance
------------------------------
The manuscript trains on 4 measured anchors + 40 *virtual* Latin-hypercube points but does not publish the
virtual responses.  Place the authors' file at ``data/optimisation/virtual_lhs_40.csv`` (columns: cu_wt, cnt_wt, rpm,
sigma, ln_Rct, ln_Icorr, ln_ICR1  or  Rct, Icorr, ICR1) to reproduce Tables 14-15.  Without it the script uses
the SYNTHETIC DEMO surface from ``adc12pem.optimisation.synthetic_truth``; all outputs are then tagged
``SYNTHETIC_DEMO`` and figures carry a red watermark.  Demo numbers are NOT manuscript results.

Usage: python scripts/07_optimisation_pipeline.py [--quick]   (--quick: 1 CV repeat, 2 BO seeds)
"""
import _bootstrap  # noqa: F401
import argparse
import json
import platform
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.stats import qmc

from adc12pem import optimisation as O, style
from adc12pem.constants import COLOURS, MARKERS, RANDOM_SEEDS, RES_OPT, SAMPLES
from adc12pem.io import load, save_fig, save_table
from adc12pem.plotting import watermark

ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); args = ap.parse_args()
style.apply()
RES_OPT.mkdir(parents=True, exist_ok=True)
SEED = 42
df, origin = O.load_training_set(SEED)
demo = origin == "SYNTHETIC_DEMO"
print(f"Training set: {len(df)} rows, origin = {origin}")
U = O.to_unit(df[O.DESIGN_VARS].values); Y = df[O.RESPONSES]
df.assign(dataset_origin=origin).to_csv(RES_OPT / "training_set_used.csv", index=False)

# ------------------------------------------------------------- 1. surrogate cross-validation (Table 14, Fig 21)
cv, preds = O.cross_validate(U, Y, seed=SEED, n_repeats=1 if args.quick else 2)
t14 = cv.pivot(index="response", columns="model", values=["R2", "RMSE", "MAE"]).reindex(O.RESPONSES)
t14.columns = [f"{m}_{k}" for k, m in t14.columns]
t14["dataset_origin"] = origin
save_table(t14, "table14_recomputed")
print(t14.round(3).to_string())

lab = {"sigma": "\u03c3 (MS m$^{-1}$)", "ln_Rct": "ln $R_{ct}$", "ln_Icorr": "ln $I_{corr}$", "ln_ICR1": "ln ICR (1.0 MPa)"}
fig, axs = plt.subplots(2, 2, figsize=(style.COL2 * 0.8, 5.0))
for ax, r, pl in zip(axs.ravel(), O.RESPONSES, "abcd"):
    T, P = preds[(r, "GP")]
    ax.scatter(T, P, s=14, facecolor="none", edgecolor="#0072B2", lw=0.8)
    lo, hi = min(T.min(), P.min()), max(T.max(), P.max()); pad = 0.05 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8)
    ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad); ax.set_aspect("equal")
    row = cv[(cv.response == r) & (cv.model == "GP")].iloc[0]
    ax.text(0.05, 0.93, f"R$^2$ = {row.R2:.2f}\nRMSE = {row.RMSE:.3f}", transform=ax.transAxes, va="top", fontsize=7)
    ax.set_xlabel(f"Observed {lab[r]}"); ax.set_ylabel(f"Predicted {lab[r]}"); style.panel_label(ax, f"({pl})", x=-0.25)
if demo:
    fig.suptitle("SYNTHETIC DEMO DATA \u2013 not manuscript results", color="#B00020", fontsize=7, fontweight="bold", y=0.995)
fig.tight_layout(); save_fig(fig, "fig21_gp_cross_validation")

# ------------------------------------------------------------- 2. final GPs, Sobol screening, Pareto set
gps = O.fit_gps(U, Y, seed=SEED)
pool_U = qmc.Sobol(d=3, scramble=True, seed=SEED).random(2 ** 14)
pool_X = O.from_unit(pool_U)
mu, sd = O.gp_predict(gps, pool_U)
Fmin = O.to_min(mu)
nd = O.non_dominated(Fmin)
f_lo, f_hi = Fmin.min(axis=0), Fmin.max(axis=0)
Fn = (Fmin - f_lo) / (f_hi - f_lo)
s4_mu, _ = O.gp_predict(gps, O.to_unit([[5.62, 0.70, 500.0]]))
s4_min = O.to_min(s4_mu)
dominates_s4 = np.all(Fmin <= s4_min, axis=1) & np.any(Fmin < s4_min, axis=1)
stats = dict(dataset_origin=origin, n_candidates=len(pool_X), n_pareto=int(nd.sum()), n_dominating_S4_all_four=int(dominates_s4.sum()),
             n_pareto_dominating_S4=int((dominates_s4 & nd).sum()))
print(stats)
par = pd.DataFrame(pool_X[nd], columns=O.DESIGN_VARS)
for j, r in enumerate(O.RESPONSES):
    par[r + "_pred"] = mu[nd, j]
par["dominates_S4"] = dominates_s4[nd]
save_table(par, "pareto_set", index=False)

# candidates (Table 15)
rows = []
for name, lam in O.PREFERENCES.items():
    g = O.chebyshev(Fn[nd], lam); k = np.where(nd)[0][np.argmin(g)]
    rows.append((name, k))
cand_idx = dict(rows)

def describe(name, X, mu_, sd_):
    d = dict(candidate=name, cu_wt=X[0], cnt_wt=X[1], rpm=X[2])
    for j, r in enumerate(O.RESPONSES):
        lo, hi = mu_[j] - 1.96 * sd_[j], mu_[j] + 1.96 * sd_[j]
        f = (lambda v: v) if r == "sigma" else np.exp
        nm = {"sigma": "sigma_MS_per_m", "ln_Rct": "Rct_kohm_cm2", "ln_Icorr": "Icorr_uA_cm2", "ln_ICR1": "ICR_1p0MPa_mOhm_cm2"}[r]
        d[nm], d[nm + "_pi95_lo"], d[nm + "_pi95_hi"] = float(f(mu_[j])), float(f(lo)), float(f(hi))
    return d

t15 = [describe(n, pool_X[k], mu[k], sd[k]) for n, k in cand_idx.items()]
m4, s4 = O.gp_predict(gps, O.to_unit([[5.62, 0.70, 500.0]]))
t15.append(describe("S4 (GP prediction)", [5.62, 0.70, 500.0], m4[0], s4[0]))
t15 = pd.DataFrame(t15); t15["dataset_origin"] = origin
save_table(t15, "table15_recomputed", index=False)
print(t15.round(3).to_string(index=False))

# ------------------------------------------------------------- 3. permutation importance (Fig 23)
imp = O.permutation_importance_oof(U, Y, seed=SEED, n_repeats=5 if args.quick else 20)
imp["dataset_origin"] = origin
save_table(imp, "permutation_importance_rf", index=False)
fig, ax = plt.subplots(figsize=(style.COL2 * 0.75, 2.8)); w = 0.26
for k, v in enumerate(O.DESIGN_VARS):
    sub = imp[imp.variable == v].set_index("response").loc[O.RESPONSES]
    ax.bar(np.arange(4) + (k - 1) * w, sub["mean"], w * 0.92, yerr=sub["sd"], capsize=2.5, label={"cu_wt": "Cu", "cnt_wt": "CNT", "rpm": "Stirring"}[v],
           color=["#0072B2", "#009E73", "#E69F00"][k], edgecolor="black", lw=0.5)
ax.set_xticks(range(4)); ax.set_xticklabels([lab[r] for r in O.RESPONSES]); ax.set_ylabel("Permutation importance (\u0394R$^2$)")
ax.legend(fontsize=6.5, ncol=3, loc="upper right")
if demo: watermark(ax, y=1.12)
fig.tight_layout(); save_fig(fig, "fig23_permutation_importance")

# ------------------------------------------------------------- 4. response surfaces at 500 rpm (Fig 22)
gc, gn = np.meshgrid(np.linspace(0, 8, 80), np.linspace(0, 1.5, 60))
grid = np.column_stack([gc.ravel(), gn.ravel(), np.full(gc.size, 500.0)])
gm, _ = O.gp_predict(gps, O.to_unit(grid))
panels = [("sigma", "\u03c3 (MS m$^{-1}$)", lambda v: v, "viridis"), ("ln_Rct", "$R_{ct}$ (k\u03a9 cm$^2$)", np.exp, "viridis"),
          ("ln_Icorr", "$I_{corr}$ (\u00b5A cm$^{-2}$)", np.exp, "viridis_r"), ("ln_ICR1", "ICR @1.0 MPa (m\u03a9 cm$^2$)", np.exp, "viridis_r")]
fig, axs = plt.subplots(2, 2, figsize=(style.COL2, 5.4))
virt = df[df.origin != "measured"]; anc = df[df.origin == "measured"].reset_index(drop=True)
for ax, (r, ttl, f, cm), pl in zip(axs.ravel(), panels, "abcd"):
    j = O.RESPONSES.index(r)
    cs = ax.contourf(gc, gn, f(gm[:, j]).reshape(gc.shape), 18, cmap=cm)
    fig.colorbar(cs, ax=ax, label=ttl, pad=0.02)
    ax.scatter(virt.cu_wt, virt.cnt_wt, s=10, facecolor="none", edgecolor="white", lw=0.6, label="Virtual LHS (all rpm)")
    for i, s_ in enumerate(SAMPLES):
        ax.scatter(anc.cu_wt[i], anc.cnt_wt[i], marker=MARKERS[s_], s=34, color=COLOURS[s_] if s_ != "S1" else "black", edgecolor="white", lw=0.8, zorder=5)
        ax.annotate(s_, (anc.cu_wt[i], anc.cnt_wt[i]), xytext=(4, 4), textcoords="offset points", fontsize=6.5, color="white", fontweight="bold")
    ax.scatter(t15.cu_wt[:3], t15.cnt_wt[:3], marker="*", s=70, color="#F0E442", edgecolor="black", lw=0.6, zorder=6, label="Candidates")
    ax.set_xlabel("Added Cu (wt.%)"); ax.set_ylabel("CNT (wt.%)"); style.panel_label(ax, f"({pl})", x=-0.2)
axs[0, 0].legend(fontsize=5.5, loc="lower right", framealpha=0.6)
if demo: fig.suptitle("SYNTHETIC DEMO DATA \u2013 not manuscript results", color="#B00020", fontsize=7, fontweight="bold", y=1.0)
fig.tight_layout(); save_fig(fig, "fig22_response_surfaces_500rpm")

# ------------------------------------------------------------- 5. Pareto projections + ParEGO (Fig 24)
oracle = lambda Uq: np.column_stack([gps[r].predict(np.atleast_2d(Uq)) for r in O.RESPONSES])
seeds = RANDOM_SEEDS[:2] if args.quick else RANDOM_SEEDS
bo, rs = [], []
for sd_ in seeds:
    a, b = O.parego(oracle, f_lo, f_hi, n_init=12, n_iter=20, seed=sd_)
    bo.append(a); rs.append(b); print(f"seed {sd_}: HV start {a[11]:.3f} -> BO {a[-1]:.3f} | random {b[-1]:.3f}")
bo, rs = np.array(bo), np.array(rs)
n_eval = np.arange(1, bo.shape[1] + 1)
hv = pd.DataFrame({"evaluation": n_eval, "parego_mean": bo.mean(0), "parego_sd": bo.std(0, ddof=1), "random_mean": rs.mean(0), "random_sd": rs.std(0, ddof=1)})
hv["dataset_origin"] = origin
save_table(hv, "hypervolume_traces", index=False)
hv_max = O.HV_REF ** 4
summary = dict(stats, hv_max=hv_max, hv_start_after_12=float(bo[:, 11].mean()), hv_parego_final_mean=float(bo[:, -1].mean()), hv_parego_final_sd=float(bo[:, -1].std(ddof=1)),
               hv_random_final_mean=float(rs[:, -1].mean()), hv_random_final_sd=float(rs[:, -1].std(ddof=1)),
               parego_pct_of_max=float(100 * bo[:, -1].mean() / hv_max), n_seeds=len(seeds), seeds=seeds)
fig, axs = plt.subplots(1, 3, figsize=(style.COL2, 2.6))
for ax, (jx, jy, xl, yl), pl in zip(axs[:2], [(0, 2, "\u03c3 (MS m$^{-1}$)", "ln $I_{corr}$"), (1, 3, "ln $R_{ct}$", "ln ICR (1.0 MPa)")], "ab"):
    ax.scatter(mu[~nd, jx], mu[~nd, jy], s=1.5, color="lightgrey", rasterized=True)
    sc = ax.scatter(mu[nd, jx], mu[nd, jy], s=6, c=pool_X[nd, 1], cmap="viridis", rasterized=True)
    for i, s_ in enumerate(SAMPLES):
        ax.scatter(anc[O.RESPONSES[jx]][i], anc[O.RESPONSES[jy]][i], marker=MARKERS[s_], s=26, color=COLOURS[s_], edgecolor="white", lw=0.6, zorder=5)
    ci = list(cand_idx.values())
    ax.scatter(mu[ci, jx], mu[ci, jy], marker="*", s=60, color="#F0E442", edgecolor="black", lw=0.6, zorder=6)
    ax.set_xlabel(xl); ax.set_ylabel(yl); style.panel_label(ax, f"({pl})", x=-0.28)
fig.colorbar(sc, ax=axs[1], label="CNT (wt.%)", pad=0.02)
ax = axs[2]
ax.fill_between(n_eval, bo.mean(0) - bo.std(0, ddof=1), bo.mean(0) + bo.std(0, ddof=1), color="#D55E00", alpha=0.25)
ax.fill_between(n_eval, rs.mean(0) - rs.std(0, ddof=1), rs.mean(0) + rs.std(0, ddof=1), color="grey", alpha=0.25)
ax.plot(n_eval, bo.mean(0), color="#D55E00", label="ParEGO-EI"); ax.plot(n_eval, rs.mean(0), color="k", ls="--", label="Random")
ax.axhline(hv_max, color="grey", ls=":", lw=0.8); ax.axvline(12.5, color="grey", ls=":", lw=0.8)
ax.set_xlabel("Evaluations"); ax.set_ylabel("Hypervolume (ref = 1.05)"); ax.legend(fontsize=6.5, loc="lower right"); style.panel_label(ax, "(c)", x=-0.28)
if demo: fig.suptitle("SYNTHETIC DEMO DATA \u2013 not manuscript results", color="#B00020", fontsize=7, fontweight="bold", y=1.0)
fig.tight_layout(); save_fig(fig, "fig24_pareto_and_bayesian_optimisation")

summary["environment"] = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__)
(RES_OPT / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
if sklearn.__version__ != "1.5.2":
    print(f"WARNING: manuscript used scikit-learn 1.5.2; running {sklearn.__version__} - numbers may differ slightly.", file=sys.stderr)
print("07 done")
