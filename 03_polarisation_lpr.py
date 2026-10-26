#!/usr/bin/env python3
"""Polarisation / LPR analysis: EIS-LPR cross-validation (Figs 15, 16) and raw-curve processing (Fig 14).

Always: Rct/Rp ratios, Rp*Icorr (Stern-Geary B), implied anodic slopes from tabulated data.
If raw curves exist in data/raw/polarisation/*.csv (columns: sample, replicate, potential_V,
current_density_A_cm2): Ecorr, cathodic-Tafel Icorr/beta_c (-0.25 to -0.12 V vs Ecorr), LPR Rp (+/-5 mV),
and Fig 14.
"""
import _bootstrap  # noqa: F401
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from adc12pem import polarisation as pol, stats as S, style
from adc12pem.constants import COLOURS, DOE_ICORR_UA_CM2, LABELS, MARKERS, RAW, SAMPLES
from adc12pem.io import load, save_fig, save_table
from adc12pem.plotting import sample_bars

style.apply()
t7, t8 = load("table07_eis_25C"), load("table08_polarisation_lpr")

# ---------------------------------------------------------------- cross-validation table
x = pd.DataFrame({
    "Rct_ohm_cm2": t7.Rct_kohm_cm2_mean * 1e3,
    "Rct_sd_ohm_cm2": t7.Rct_kohm_cm2_sd * 1e3,
    "Rp_ohm_cm2": t8.Rp_ohm_cm2,
    "Icorr_uA_cm2": t8.Icorr_uA_cm2_mean,
})
x["Rct_over_Rp"] = x.Rct_ohm_cm2 / x.Rp_ohm_cm2
x["B_mV_(Rp*Icorr)"] = x.Rp_ohm_cm2 * x.Icorr_uA_cm2 * 1e-3
x["beta_a_implied_mV_dec"] = [pol.implied_beta_a(b, 120.0) for b in x["B_mV_(Rp*Icorr)"]]
x["within_10pct"] = (x.Rct_over_Rp - 1).abs() <= 0.10
x["Icorr_if_B_plus30pct"] = x.Icorr_uA_cm2 * 1.3
x["Icorr_if_B_minus30pct"] = x.Icorr_uA_cm2 * 0.7
save_table(x, "eis_lpr_crossvalidation")
print(x.round(3).to_string())

# ---------------------------------------------------------------- Fig 15
_, _, let_rct = S.tukey_from_summary(SAMPLES, t7.Rct_kohm_cm2_mean.values, t7.Rct_kohm_cm2_sd.values, 3)
m, s = np.log(t8.Icorr_uA_cm2_mean), t8.Icorr_uA_cm2_sd / t8.Icorr_uA_cm2_mean
_, _, let_ic = S.tukey_from_summary(SAMPLES, m.values, s.values, 3)
fig, axs = plt.subplots(1, 2, figsize=(style.COL2, 2.7))
sample_bars(axs[0], t7.Rct_kohm_cm2_mean, t7.Rct_kohm_cm2_sd, let_rct, "$R_{ct}$ (k\u03a9 cm$^2$)", ylim=(0, 5))
sample_bars(axs[1], t8.Icorr_uA_cm2_mean, t8.Icorr_uA_cm2_sd, let_ic, "$I_{corr}$ (\u00b5A cm$^{-2}$)", ylim=(0, 6))
axs[1].axhline(DOE_ICORR_UA_CM2, color="#B00020", ls="--", lw=1, label="DOE target (1 \u00b5A cm$^{-2}$)")
axs[1].legend(loc="upper right", fontsize=6.5)
for ax, lab in zip(axs, ["(a)", "(b)"]):
    style.panel_label(ax, lab, x=-0.2)
fig.tight_layout(); save_fig(fig, "fig15_rct_icorr")

# ---------------------------------------------------------------- Fig 16
fig, axs = plt.subplots(1, 2, figsize=(style.COL2, 2.9), gridspec_kw=dict(width_ratios=[1, 1]))
ax = axs[0]
lim = np.array([1.2, 4.9])
ax.fill_between(lim, lim * 0.9, lim * 1.1, color="grey", alpha=0.2, label="\u00b110% band")
ax.plot(lim, lim, "k--", lw=0.9, label="1:1")
for s_ in SAMPLES:
    ax.errorbar(x.loc[s_, "Rp_ohm_cm2"] / 1e3, x.loc[s_, "Rct_ohm_cm2"] / 1e3, yerr=x.loc[s_, "Rct_sd_ohm_cm2"] / 1e3,
                fmt=MARKERS[s_], color=COLOURS[s_], capsize=2.5, label=f"{s_}")
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
ax.set_xlabel("$R_p$ (LPR) (k\u03a9 cm$^2$)"); ax.set_ylabel("$R_{ct}$ (EIS) (k\u03a9 cm$^2$)"); ax.legend(fontsize=6, loc="upper left")
ax = axs[1]
ax.axhspan(0.9, 1.1, color="grey", alpha=0.2); ax.axhline(1, color="k", ls="--", lw=0.9)
ax.bar(range(4), x.Rct_over_Rp, color=[COLOURS[s_] for s_ in SAMPLES], edgecolor="black", lw=0.6, alpha=0.85)
for i, v in enumerate(x.Rct_over_Rp):
    ax.text(i, v + 0.01, f"{v:.2f}", ha="center", fontsize=7)
ax.set_ylim(0.8, 1.2); ax.set_xticks(range(4)); ax.set_xticklabels(SAMPLES); ax.set_ylabel("$R_{ct}$ / $R_p$")
for a_, lab in zip(axs, ["(a)", "(b)"]):
    style.panel_label(a_, lab, x=-0.2)
fig.tight_layout(); save_fig(fig, "fig16_eis_lpr_crossvalidation")

# ---------------------------------------------------------------- raw curves
files = [f for f in sorted(glob.glob(str(RAW / "polarisation" / "*.csv"))) if "TEMPLATE" not in f.upper()]
if not files:
    print("No raw polarisation files in data/raw/polarisation - skipping curve analysis.")
else:
    recs, fig, ax = [], *plt.subplots(figsize=(style.COL1 * 1.3, 3.0))
    for fn in files:
        d = pd.read_csv(fn)
        for (s_, rep), g in d.groupby(["sample", "replicate"]):
            E, i = g.potential_V.values, g.current_density_A_cm2.values
            ec = pol.find_ecorr(E, i)
            tf = pol.tafel_cathodic(E, i, ec)
            rp = pol.lpr_rp(E, i, ec)
            recs.append(dict(sample=s_, replicate=rep, Ecorr_V=ec, Icorr_uA_cm2=tf["icorr"] * 1e6, beta_c_mV_dec=tf["beta_c"] * 1e3,
                             Rp_ohm_cm2=rp, B_mV=rp * tf["icorr"] * 1e3))
            if rep == g.replicate.min():
                ax.semilogx(np.abs(i) * 1e6 + 1e-3, E - ec, color=COLOURS[s_], label=s_ if rep == g.replicate.min() else None)
    ax.axvline(DOE_ICORR_UA_CM2, color="#B00020", ls="--", lw=0.9)
    ax.set_xlabel("|i| (\u00b5A cm$^{-2}$)"); ax.set_ylabel("E \u2212 $E_{corr}$ (V vs Ag/AgCl)"); ax.legend(fontsize=6)
    save_fig(fig, "fig14_polarisation_curves")
    r = pd.DataFrame(recs); save_table(r, "polarisation_fits_raw", index=False)
    save_table(r.groupby("sample").agg(["mean", "std"]), "polarisation_summary_raw")
    print(r.round(4).to_string(index=False))
print("03 done")
