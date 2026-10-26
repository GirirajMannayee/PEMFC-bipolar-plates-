#!/usr/bin/env python3
"""Interfacial contact resistance (Table 9, Fig 17): quadratic log-log interpolation to 1.4 MPa,
secant and local exponents, per-interface vs two-interface basis.
Raw ICR data (data/raw/icr/*.csv: sample, replicate, pressure_MPa, icr_mOhm_cm2) are processed if present.
"""
import _bootstrap  # noqa: F401
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from adc12pem import icr, style
from adc12pem.constants import COLOURS, DOE_ICR_MOHM_CM2, DOE_PRESSURE_MPA, LABELS, MARKERS, RAW, SAMPLES
from adc12pem.io import load, save_fig, save_table

style.apply()
t9 = load("table09_icr_vs_pressure")
P = np.array([0.5, 1.0, 1.5, 2.0])
cols = ["icr_0p5MPa", "icr_1p0MPa", "icr_1p5MPa", "icr_2p0MPa"]
rows = []
fig, ax = plt.subplots(figsize=(style.COL1 * 1.35, 3.0))
Pg = np.linspace(0.45, 2.05, 200)
for s in SAMPLES:
    y = np.array([t9.loc[s, c + "_mean"] for c in cols]); sd = np.array([t9.loc[s, c + "_sd"] for c in cols])
    coef = icr.quad_loglog_fit(P, y)
    loc = icr.local_exponents(P, y)
    v14 = float(icr.quad_loglog_predict(coef, DOE_PRESSURE_MPA))
    rows.append(dict(sample=s, a2=coef[0], a1=coef[1], a0=coef[2], ICR_1p4_fit=v14, ICR_1p4_reported=t9.loc[s, "icr_1p4MPa_interpolated"],
                     ICR_1p4_two_interface=2 * v14, secant_b=icr.secant_exponent(P, y), local_b_0p5_1p0=loc[0], local_b_1p0_1p5=loc[1],
                     local_b_1p5_2p0=loc[2], ratio_to_S1_at_each_P=";".join(f"{a / b:.3f}" for a, b in zip(y, [t9.loc['S1', c + '_mean'] for c in cols]))))
    ax.errorbar(P, y, yerr=sd, fmt=MARKERS[s], color=COLOURS[s], capsize=2.5, label=f"{s} {LABELS[s]}")
    ax.plot(Pg, icr.quad_loglog_predict(coef, Pg), color=COLOURS[s], lw=1)
ax.axvline(DOE_PRESSURE_MPA, color="grey", ls=":", lw=1)
ax.text(DOE_PRESSURE_MPA + 0.02, ax.get_ylim()[1] * 0.97, "1.4 MPa\n(138 N cm$^{-2}$)", fontsize=6.5, va="top")
ax.text(0.97, 0.52, f"DOE target {DOE_ICR_MOHM_CM2:.0f} m\u03a9 cm$^2$\n(off scale)", transform=ax.transAxes, ha="right", fontsize=6.5, color="#B00020")
ax.set_xlabel("Compaction pressure (MPa)"); ax.set_ylabel("ICR per specimen\u2013GDL interface (m\u03a9 cm$^2$)")
ax.set_xlim(0.4, 2.1); ax.set_ylim(0, 7); ax.legend(fontsize=6, loc="upper right")
fig.tight_layout(); save_fig(fig, "fig17_icr_vs_pressure")
out = pd.DataFrame(rows).set_index("sample")
save_table(out, "icr_analysis")
print(out.round(3).to_string())
print("NOTE: ICR(S2..S4)/ICR(S1) is the same at every pressure -> identical exponents (flagged in manuscript Sec. 4).")

files = [f for f in sorted(glob.glob(str(RAW / "icr" / "*.csv"))) if "TEMPLATE" not in f.upper()]
if files:
    d = pd.concat([pd.read_csv(f) for f in files])
    recs = []
    for (s, rep), g in d.groupby(["sample", "replicate"]):
        g = g.sort_values("pressure_MPa"); c = icr.quad_loglog_fit(g.pressure_MPa.values, g.icr_mOhm_cm2.values)
        recs.append(dict(sample=s, replicate=rep, ICR_1p4=float(icr.quad_loglog_predict(c, 1.4)),
                         secant_b=icr.secant_exponent(g.pressure_MPa.values, g.icr_mOhm_cm2.values)))
    r = pd.DataFrame(recs); save_table(r, "icr_fits_raw", index=False); save_table(r.groupby("sample").agg(["mean", "std"]), "icr_summary_raw")
print("04 done")
