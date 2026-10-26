#!/usr/bin/env python3
"""Integrated benchmark (Tables 12-13), conductivity physical-bound check (Sec. 3.3) and Fig 20 scorecard.

Scorecard method (inferred - the manuscript does not state the formula): each indicator is divided by the
best value among S1-S4 (1/Icorr and 1/ICR used so that larger is better); composite = geometric mean.
This reproduces the reported composite scores (0.49 / 0.63 / 0.61 / 1.00) within rounding.
"""
import _bootstrap  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from adc12pem import style, verify
from adc12pem.constants import LABELS, SAMPLES
from adc12pem.io import load, save_fig, save_table

style.apply()
t6, t7, t8, t9, t10 = (load(n) for n in ["table06_conductivity_hardness", "table07_eis_25C", "table08_polarisation_lpr",
                                         "table09_icr_vs_pressure", "table10_temperature_arrhenius"])
# ---------------------------------------------------------------- Table 12 (assembled)
t12 = pd.DataFrame({
    "sigma_MS_per_m": t6.sigma_MS_per_m_mean, "hardness_HV": t6.hardness_HV_mean,
    "Eocp_V": t8.Ecorr_V_vs_AgAgCl_mean, "Rct_25C_kohm_cm2": t7.Rct_kohm_cm2_mean, "Rct_80C_kohm_cm2": t10.Rct_80C_kohm_cm2,
    "Ea_kJ_mol": t10.Ea_kJ_mol_mean, "Icorr_uA_cm2": t8.Icorr_uA_cm2_mean, "ICR_1p4MPa_mOhm_cm2": t9.icr_1p4MPa_interpolated}).T
abs_rows = {"Eocp_V": ("mV", 1000.0), "Ea_kJ_mol": ("kJ/mol", 1.0)}   # differences, not percentages
t12["S4_vs_S1"] = [f"{(t12.loc[r, 'S4'] - t12.loc[r, 'S1']) * abs_rows[r][1]:+.1f} {abs_rows[r][0]}" if r in abs_rows
                   else f"{100 * (t12.loc[r, 'S4'] / t12.loc[r, 'S1'] - 1):+.0f}%" for r in t12.index]
save_table(t12, "table12_assembled")
print(t12.round(2))

# ---------------------------------------------------------------- scorecard (Fig 20)
raw, norm, comp = verify.scorecard(t6, t7, t8, t9)
sc = norm.copy(); sc["composite"] = comp
save_table(sc, "scorecard_normalised")
print(sc.round(3))
names = {"sigma": "\u03c3", "hardness": "Hardness", "Rct": "$R_{ct}$", "inv_Icorr": "1/$I_{corr}$", "inv_ICR14": "1/ICR\n(1.4 MPa)"}
delta = 100 * (raw / raw.loc["S1"] - 1)
fig, ax = plt.subplots(figsize=(style.COL2 * 0.85, 2.6))
data = np.column_stack([norm.values, comp.values])
im = ax.imshow(data, cmap="YlGnBu", vmin=0.2, vmax=1.0, aspect="auto")
for i, s in enumerate(SAMPLES):
    for j in range(6):
        txt = f"{data[i, j]:.2f}" + (f"\n({delta.iloc[i, j]:+.0f}%)" if j < 5 and s != "S1" else "")
        ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="white" if data[i, j] > 0.7 else "black")
ax.set_xticks(range(6)); ax.set_xticklabels([names[c] for c in norm.columns] + ["Composite\n(geo. mean)"], fontsize=7)
ax.set_yticks(range(4)); ax.set_yticklabels([f"{s} {LABELS[s]}" for s in SAMPLES], fontsize=7)
ax.axvline(4.5, color="white", lw=3); fig.colorbar(im, ax=ax, label="Score (value / best)", fraction=0.04, pad=0.02)
fig.text(0.5, -0.02, "Cells: normalised score; in brackets, change relative to S1.", ha="center", fontsize=6.5)
fig.tight_layout(); save_fig(fig, "fig20_scorecard")

# ---------------------------------------------------------------- conductivity physical bound (Sec. 3.3)
wb, f = verify.wiener_bound()
bound = pd.DataFrame([dict(
    quantity=q, value=v) for q, v in [
    ("Cu volume fraction (%) [assumed rho_Cu=8.96, rho_ADC12=2.74 g/cm3]", 100 * f[1]),
    ("CNT volume fraction (%) [assumed rho_CNT=1.75 g/cm3]", 100 * f[2]),
    ("Wiener upper bound sigma (MS/m) for S4", wb),
    ("Bound relative to S1 (%)", 100 * (wb / 16.8 - 1)),
    ("Measured S4 sigma (MS/m)", 24.6),
    ("Measured S4 relative to S1 (%)", 100 * (24.6 / 16.8 - 1)),
    ("Measured S4 minus bound (MS/m)", 24.6 - wb)]])
save_table(bound, "conductivity_physical_bound", index=False)
print(bound.round(2).to_string(index=False))
print("06 done")
