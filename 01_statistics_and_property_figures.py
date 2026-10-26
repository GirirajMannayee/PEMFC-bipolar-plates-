#!/usr/bin/env python3
"""Statistics (ANOVA, Tukey HSD, 2x2 factorial) and Figs 10, 19.

Inputs : data/processed/table06/07/08/09 (mean +/- SD, n = 3)
Outputs: results/tables/anova_tukey.csv, factorial_effects.csv, tukey_pairwise_*.csv
         results/figures/fig10_conductivity_hardness.(png|pdf), fig19_factorial.(png|pdf)
"""
import _bootstrap  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter, NullFormatter

from adc12pem import stats as S, style
from adc12pem.constants import COLOURS, LABELS, SAMPLES
from adc12pem.io import load, save_fig, save_table
from adc12pem.plotting import sample_bars

style.apply()
N = 3
t6, t7, t8, t9 = (load(n) for n in ["table06_conductivity_hardness", "table07_eis_25C",
                                     "table08_polarisation_lpr", "table09_icr_vs_pressure"])

# ------------------------------------------------------------------ one-way ANOVA + Tukey
responses = {
    "sigma (MS/m)": (t6.sigma_MS_per_m_mean, t6.sigma_MS_per_m_sd, False),
    "hardness (HV)": (t6.hardness_HV_mean, t6.hardness_HV_sd, False),
    "Rct (kOhm cm2) [linear scale, Table 7 groups]": (t7.Rct_kohm_cm2_mean, t7.Rct_kohm_cm2_sd, False),
    "ln Rct": (t7.Rct_kohm_cm2_mean, t7.Rct_kohm_cm2_sd, True),
    "ln Icorr": (t8.Icorr_uA_cm2_mean, t8.Icorr_uA_cm2_sd, True),
}
rows, letters = [], {}
for name, (m, s, lg) in responses.items():
    mm, ss = (np.log(m), s / m) if lg else (m, s)
    a = S.oneway_from_summary(mm.values, ss.values, N)
    pw, hsd, let = S.tukey_from_summary(SAMPLES, mm.values, ss.values, N)
    letters[name] = let
    rows.append(dict(response=name, F=a["F"], df_between=a["dfb"], df_within=a["dfw"], p=a["p"], eta2=a["eta2"],
                     tukey_HSD=hsd, **{f"group_{k}": v for k, v in let.items()}))
    save_table(pw, f"tukey_pairwise_{name.split(' ')[0].replace('/', '')}_{'ln' if lg else 'lin'}", index=False)
anova = pd.DataFrame(rows)
save_table(anova, "anova_tukey", index=False)
print(anova.round(4).to_string(index=False))

# ------------------------------------------------------------------ 2x2 factorial (Table 11)
fac = {
    "sigma_MS_per_m": (t6.sigma_MS_per_m_mean, t6.sigma_MS_per_m_sd, False),
    "hardness_HV": (t6.hardness_HV_mean, t6.hardness_HV_sd, False),
    "Rct": (t7.Rct_kohm_cm2_mean, t7.Rct_kohm_cm2_sd, True),
    "Icorr": (t8.Icorr_uA_cm2_mean, t8.Icorr_uA_cm2_sd, True),
    "ICR_1p0MPa": (t9.icr_1p0MPa_mean, t9.icr_1p0MPa_sd, True),
}
res, rows = {}, []
for k, (m, s, lg) in fac.items():
    res[k] = S.factorial_2x2(m.values, s.values, N, log=lg)
p_holm = dict(zip(fac, S.holm([res[k]["interaction"]["p"] for k in fac])))
for k, (m, s, lg) in fac.items():
    for eff in ["cu", "cnt", "interaction"]:
        r = res[k][eff]
        rows.append(dict(response=k, scale="ratio" if lg else "additive", effect=eff, estimate=r["estimate"],
                         ci95_lo=r["lo"], ci95_hi=r["hi"], t=r["t"], p=r["p"],
                         p_holm=p_holm[k] if eff == "interaction" else np.nan, eta2=res[k]["eta2"]))
fe = pd.DataFrame(rows)
save_table(fe, "factorial_effects", index=False)
print(fe.round(4).to_string(index=False))

# ------------------------------------------------------------------ Fig 10
fig, axs = plt.subplots(1, 2, figsize=(style.COL2, 2.7))
sample_bars(axs[0], t6.sigma_MS_per_m_mean, t6.sigma_MS_per_m_sd, letters["sigma (MS/m)"],
            "Electrical conductivity, \u03c3 (MS m$^{-1}$)", ylim=(0, 30))
sample_bars(axs[1], t6.hardness_HV_mean, t6.hardness_HV_sd, letters["hardness (HV)"], "Vickers hardness (HV)", ylim=(0, 140))
for ax, lab in zip(axs, "ab"):
    style.panel_label(ax, f"({lab})", x=-0.2)
sec = axs[0].secondary_yaxis("right", functions=(lambda v: v / 58 * 100, lambda v: v * 58 / 100))
sec.set_ylabel("%IACS")
fig.text(0.5, -0.04, "Bars: mean \u00b1 SD (n = 3). Letters: Tukey HSD groups (\u03b1 = 0.05); shared letter = not significantly different.",
         ha="center", fontsize=6.5)
fig.tight_layout()
save_fig(fig, "fig10_conductivity_hardness")

# ------------------------------------------------------------------ Fig 19
fig = plt.figure(figsize=(style.COL2, 5.6))
gs = GridSpec(2, 5, height_ratios=[1, 1.15], hspace=0.6, wspace=0.95, left=0.075, right=0.985, top=0.93, bottom=0.09)
cell = {
    "sigma_MS_per_m": ("\u03c3 (MS m$^{-1}$)", t6.sigma_MS_per_m_mean, t6.sigma_MS_per_m_sd, False),
    "hardness_HV": ("Hardness (HV)", t6.hardness_HV_mean, t6.hardness_HV_sd, False),
    "Rct": ("$R_{ct}$ (k\u03a9 cm$^2$)", t7.Rct_kohm_cm2_mean, t7.Rct_kohm_cm2_sd, True),
    "Icorr": ("$I_{corr}$ (\u00b5A cm$^{-2}$)", t8.Icorr_uA_cm2_mean, t8.Icorr_uA_cm2_sd, True),
    "ICR_1p0MPa": ("ICR @1.0 MPa (m\u03a9 cm$^2$)", t9.icr_1p0MPa_mean, t9.icr_1p0MPa_sd, True),
}
for j, (k, (lab, m, s, lg)) in enumerate(cell.items()):
    ax = fig.add_subplot(gs[0, j])
    for cnt, (a_, b_), ls, col in [("no CNT", ("S1", "S2"), "--", "#0072B2"), ("+ CNT", ("S3", "S4"), "-", "#D55E00")]:
        ax.errorbar([0, 1], [m[a_], m[b_]], yerr=[s[a_], s[b_]], ls=ls, marker="o", color=col, capsize=2.5, label=cnt)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["\u2013Cu", "+Cu"]); ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel(lab, fontsize=6.5, labelpad=2)
    if lg:
        ax.set_yscale("log")
        lo_, hi_ = min(m[["S1", "S2", "S3", "S4"]] - s[["S1", "S2", "S3", "S4"]]), max(m[["S1", "S2", "S3", "S4"]] + s[["S1", "S2", "S3", "S4"]])
        tk = [t for t in [1, 1.5, 2, 3, 4, 5] if lo_ * 0.97 <= t <= hi_ * 1.03]
        ax.set_yticks(tk)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
    if j == 0:
        ax.legend(loc="lower right", fontsize=6)
        style.panel_label(ax, "(a)", x=-0.55)
# (b) forest plot of ratio responses
axb = fig.add_subplot(gs[1, :3])
ypos, yt, yl = 0, [], []
for k in ["Rct", "Icorr", "ICR_1p0MPa"]:
    for eff, ec, mk in [("cu", "#0072B2", "o"), ("cnt", "#009E73", "s"), ("interaction", "#D55E00", "D")]:
        r = res[k][eff]
        axb.errorbar(r["estimate"], ypos, xerr=[[r["estimate"] - r["lo"]], [r["hi"] - r["estimate"]]], fmt=mk, color=ec, capsize=2.5,
                     label={"cu": "Cu main", "cnt": "CNT main", "interaction": "Cu\u00d7CNT interaction"}[eff] if k == "Rct" else None)
        ypos -= 1
    yt.append(ypos + 2); yl.append({"Rct": "$R_{ct}$", "Icorr": "$I_{corr}$", "ICR_1p0MPa": "ICR (1.0 MPa)"}[k]); ypos -= 0.7
axb.axvline(1, color="grey", lw=0.8, ls=":")
axb.set_xscale("log"); axb.set_yticks(yt); axb.set_yticklabels(yl)
axb.set_xlabel("Effect ratio (95% CI); 1 = no effect")
axb.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=6, columnspacing=1.0, handletextpad=0.3)
axb.set_xlim(0.45, 2.0); axb.set_xticks([0.5, 0.7, 1, 1.5, 2])
axb.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}")); axb.xaxis.set_minor_formatter(NullFormatter())
style.panel_label(axb, "(b)", x=-0.14)
# (c) additive responses relative to ADC12
axc = fig.add_subplot(gs[1, 3:])
w = 0.25
for ii, k in enumerate(["sigma_MS_per_m", "hardness_HV"]):
    base = cell[k][1]["S1"]
    for jj, (eff, ec) in enumerate([("cu", "#0072B2"), ("cnt", "#009E73"), ("interaction", "#D55E00")]):
        r = res[k][eff]
        axc.bar(ii + (jj - 1) * w, 100 * r["estimate"] / base, w * 0.9, color=ec, alpha=0.85, edgecolor="black", lw=0.5,
                yerr=[[100 * (r["estimate"] - r["lo"]) / base], [100 * (r["hi"] - r["estimate"]) / base]], capsize=2.5)
from matplotlib.patches import Patch
axc.legend(handles=[Patch(color=c_, label=l_) for c_, l_ in [("#0072B2", "Cu main"), ("#009E73", "CNT main"), ("#D55E00", "Cu\u00d7CNT")]],
           loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=6, columnspacing=1.0, handletextpad=0.3)
axc.axhline(0, color="black", lw=0.8)
axc.set_xticks([0, 1]); axc.set_xticklabels(["\u03c3", "Hardness"]); axc.set_ylabel("Effect (% of ADC12 value)")
style.panel_label(axc, "(c)", x=-0.2)
save_fig(fig, "fig19_factorial")
print("01 done")
