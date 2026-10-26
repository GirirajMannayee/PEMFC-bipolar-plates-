#!/usr/bin/env python3
"""EIS analysis: Brug capacitance check, model Nyquist/Bode curves (Figs 12, 13) and raw-data fitting.

Two modes
---------
* Always: recompute Ceff (Eq. 2) and film thickness from the tabulated Table 7 parameters, and draw
  *model* Nyquist/Bode curves from the tabulated Rs, Rct, Y0, n (NOT measured spectra).
* If raw spectra exist in data/raw/eis/*.csv (see data/raw/eis/README.md): Lin-KK validation (Sec. 2.6),
  complex NLLS fit of Rs(Rct||CPE) with 1/|Z| weighting, per-sample summary.
"""
import _bootstrap  # noqa: F401
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

from adc12pem import eis, style
from adc12pem.constants import COLOURS, LABELS, MARKERS, RAW, SAMPLES
from adc12pem.io import load, save_fig, save_table

style.apply()
t7 = load("table07_eis_25C")
EPS0 = 8.854e-14  # F/cm

# ---------------------------------------------------------------- tabulated-parameter checks
rows = []
for s, r in t7.iterrows():
    Rs, Rct, Y0, n = r.Rs_ohm_cm2_mean, r.Rct_kohm_cm2_mean * 1e3, r["Y0_1e-5_S_sn_cm2"] * 1e-5, r.n_cpe
    ceff = eis.brug_ceff(Rs, Rct, Y0, n)
    f = np.logspace(5, -2, 400)
    z = eis.z_randles_cpe(f, Rs, Rct, Y0, n)
    phase = np.degrees(np.arctan2(z.imag, z.real))
    rows.append(dict(sample=s, Ceff_uF_cm2_from_Brug=ceff * 1e6, Ceff_uF_cm2_reported=r.Ceff_uF_cm2_mean,
                     film_thickness_nm_eps10=10 * EPS0 / ceff * 1e7, peak_phase_deg_model=float(np.max(-phase)),
                     f_at_peak_phase_Hz=float(f[np.argmax(-phase)])))
chk = pd.DataFrame(rows).set_index("sample")
save_table(chk, "eis_table7_consistency")
print(chk.round(3))

# ---------------------------------------------------------------- model Nyquist / Bode (Figs 12, 13)
f = np.logspace(5, -2, 500)
Z = {s: eis.z_randles_cpe(f, t7.loc[s, "Rs_ohm_cm2_mean"], t7.loc[s, "Rct_kohm_cm2_mean"] * 1e3,
                          t7.loc[s, "Y0_1e-5_S_sn_cm2"] * 1e-5, t7.loc[s, "n_cpe"]) for s in SAMPLES}
NOTE = "Model curves from tabulated Table 7 parameters (not measured spectra)"
fig, axs = plt.subplots(1, 2, figsize=(style.COL2, 3.0))
for ax, sel, lab in [(axs[0], SAMPLES, "(a)"), (axs[1], ["S1", "S2", "S3"], "(b)")]:
    for s in sel:
        ax.plot(Z[s].real / 1e3, -Z[s].imag / 1e3, color=COLOURS[s], label=f"{s} {LABELS[s]}")
    ax.set_xlabel("$Z'$ (k\u03a9 cm$^2$)"); ax.set_ylabel("\u2212$Z''$ (k\u03a9 cm$^2$)")
    xmax = max(Z[s].real.max() for s in sel) / 1e3 * 1.04
    ax.set_xlim(0, xmax); ax.set_ylim(0, xmax * 0.6)
    ax.set_aspect("equal", adjustable="box"); style.panel_label(ax, lab, x=-0.2)
axs[0].legend(loc="upper right", fontsize=6)
fig.text(0.5, -0.03, NOTE, ha="center", fontsize=6.5, color="#B00020")
fig.tight_layout(); save_fig(fig, "fig12_nyquist_model")

fig, axs = plt.subplots(1, 2, figsize=(style.COL2, 2.7))
for s in SAMPLES:
    axs[0].loglog(f, np.abs(Z[s]), color=COLOURS[s], label=f"{s}")
    axs[1].semilogx(f, -np.degrees(np.angle(Z[s])), color=COLOURS[s])
axs[0].set_xlabel("Frequency (Hz)"); axs[0].set_ylabel("|Z| (\u03a9 cm$^2$)")
axs[1].set_xlabel("Frequency (Hz)"); axs[1].set_ylabel("\u2212Phase angle (\u00b0)")
axs[0].legend(fontsize=6)
for ax, lab in zip(axs, ["(a)", "(b)"]):
    style.panel_label(ax, lab, x=-0.2)
fig.text(0.5, -0.04, NOTE, ha="center", fontsize=6.5, color="#B00020")
fig.tight_layout(); save_fig(fig, "fig13_bode_model")

# ---------------------------------------------------------------- raw-data pipeline
files = sorted(glob.glob(str(RAW / "eis" / "*.csv")))
files = [x for x in files if "TEMPLATE" not in x.upper()]
if not files:
    print("No raw EIS files in data/raw/eis - skipping fitting (see data/raw/eis/README.md).")
else:
    recs = []
    for fn in files:
        d = pd.read_csv(fn)
        d["Z"] = d.Zre_ohm_cm2 + 1j * d.Zim_ohm_cm2
        for (s, rep, T), g in d.groupby(["sample", "replicate", "temperature_C"]):
            g = g.sort_values("freq_Hz", ascending=False)
            kk = eis.lin_kk(g.freq_Hz.values, g.Z.values)
            fit = eis.fit_randles_cpe(g.freq_Hz.values, g.Z.values)
            recs.append(dict(sample=s, replicate=rep, temperature_C=T, KK_M=kk["M"], KK_max_resid_pct=kk["max_abs_residual_pct"],
                             KK_pass=kk["passes_1pct"], Rs=fit["Rs"], Rct_ohm_cm2=fit["Rct"], Y0=fit["Y0"], n=fit["n"],
                             chi2=fit["chi2"], Ceff_uF_cm2=fit["Ceff"] * 1e6, source=fn.split("/")[-1]))
    fits = pd.DataFrame(recs)
    save_table(fits, "eis_fits_raw", index=False)
    ok = fits[fits.KK_pass]
    summ = ok[ok.temperature_C == 25].groupby("sample").agg(["mean", "std"])
    save_table(summ, "eis_summary_25C_raw")
    print(fits.round(4).to_string(index=False))
print("02 done")
