#!/usr/bin/env python3
"""Temperature dependence (Table 10, Fig 18).

Only Rct at 25 and 80 degC is tabulated in the manuscript, so Fig 18 is drawn from those two points
(two-point apparent Ea).  The 40 and 60 degC points - and therefore the 4-point regression Ea and R2 -
need the raw EIS fits (results/tables/eis_fits_raw.csv, produced by 02_eis_analysis.py).
"""
import _bootstrap  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from adc12pem import eis, style
from adc12pem.constants import COLOURS, LABELS, MARKERS, RES_TABLES, SAMPLES
from adc12pem.io import load, save_fig, save_table

style.apply()
t10 = load("table10_temperature_arrhenius")
rows = []
for s, r in t10.iterrows():
    rows.append(dict(sample=s, Ea_reported_kJ_mol=r.Ea_kJ_mol_mean, Ea_two_point_kJ_mol=eis.ea_two_point(25, r.Rct_25C_kohm_cm2, 80, r.Rct_80C_kohm_cm2),
                     ratio_80_over_25=r.Rct_80C_kohm_cm2 / r.Rct_25C_kohm_cm2))
tab = pd.DataFrame(rows).set_index("sample")
save_table(tab, "temperature_dependence")
adv = pd.DataFrame({"temperature_C": [25, 80], "S4_vs_S1_advantage_pct": [
    100 * (t10.loc["S4", "Rct_25C_kohm_cm2"] / t10.loc["S1", "Rct_25C_kohm_cm2"] - 1),
    100 * (t10.loc["S4", "Rct_80C_kohm_cm2"] / t10.loc["S1", "Rct_80C_kohm_cm2"] - 1)]})
save_table(adv, "temperature_advantage_S4_vs_S1", index=False)
print(tab.round(2).to_string()); print(adv.round(1).to_string(index=False))

fig, axs = plt.subplots(1, 2, figsize=(style.COL2, 2.9))
Tg = np.linspace(25, 80, 50)
for s, r in t10.iterrows():
    T = np.array([25.0, 80.0]); R = np.array([r.Rct_25C_kohm_cm2, r.Rct_80C_kohm_cm2])
    axs[0].plot(T, R, MARKERS[s], color=COLOURS[s], label=f"{s} {LABELS[s]}")
    ea = tab.loc[s, "Ea_two_point_kJ_mol"] * 1000
    Rg = R[0] * np.exp(ea / 8.314462618 * (1 / (Tg + 273.15) - 1 / 298.15))
    axs[0].plot(Tg, Rg, color=COLOURS[s], lw=1)
    x = 1000 / (T + 273.15); y = np.log(1 / (R * 1e3))
    axs[1].plot(x, y, MARKERS[s], color=COLOURS[s])
    axs[1].plot(1000 / (Tg + 273.15), np.log(1 / (Rg * 1e3)), color=COLOURS[s], lw=1,
                label=f"{s}: $E_a$ = {r.Ea_kJ_mol_mean:.1f} kJ mol$^{{-1}}$")
axs[0].set_xlabel("Temperature (\u00b0C)"); axs[0].set_ylabel("$R_{ct}$ (k\u03a9 cm$^2$)"); axs[0].legend(fontsize=6)
axs[1].set_xlabel("1000 / T (K$^{-1}$)"); axs[1].set_ylabel("ln(1/$R_{ct}$) [\u03a9$^{-1}$ cm$^{-2}$]"); axs[1].legend(fontsize=6, loc="lower left")
for ax, lab in zip(axs, ["(a)", "(b)"]):
    style.panel_label(ax, lab, x=-0.2)
fig.text(0.5, -0.04, "Markers: tabulated 25 and 80 \u00b0C values; lines: two-point Arrhenius interpolation (40/60 \u00b0C data not tabulated).",
         ha="center", fontsize=6.5, color="#B00020")
fig.tight_layout(); save_fig(fig, "fig18_temperature_arrhenius")

f = RES_TABLES / "eis_fits_raw.csv"
if f.exists():
    d = pd.read_csv(f); d = d[d.KK_pass]
    recs = []
    for (s, rep), g in d.groupby(["sample", "replicate"]):
        if g.temperature_C.nunique() >= 3:
            a = eis.arrhenius_fit(g.temperature_C.values, g.Rct_ohm_cm2.values); recs.append(dict(sample=s, replicate=rep, **a))
    if recs:
        r = pd.DataFrame(recs); save_table(r, "arrhenius_fits_raw", index=False); save_table(r.groupby("sample").agg(["mean", "std"]), "arrhenius_summary_raw")
print("05 done")
