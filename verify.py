"""Recompute the manuscript's derived numbers from the tabulated data and compare.

Every check is (id, description, computed, reported, tolerance).  A check passes
if |computed - reported| <= tolerance.  Run ``python scripts/08_verify_manuscript_numbers.py``
for a markdown report, or ``pytest`` to assert that all checks pass.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from . import eis, icr, polarisation as pol, stats as S
from .constants import PROCESSED, SAMPLES, IACS_MS_PER_M
from .io import load


def scorecard(t6, t7, t8, t9):
    """Normalised indicators (value / best) and composite (geometric mean)."""
    raw = pd.DataFrame({
        "sigma": t6.sigma_MS_per_m_mean, "hardness": t6.hardness_HV_mean,
        "Rct": t7.Rct_kohm_cm2_mean, "inv_Icorr": 1.0 / t8.Icorr_uA_cm2_mean,
        "inv_ICR14": 1.0 / t9.icr_1p4MPa_interpolated})
    norm = raw / raw.max()
    comp = np.exp(np.log(norm).mean(axis=1))
    return raw, norm, comp


def wiener_bound(sigma_matrix=16.8, cu_wt=5.62, cnt_wt=0.70, rho_matrix=2.74, rho_cu=8.96, rho_cnt=1.75,
                 sigma_cu=58.0, sigma_cnt=100.0):
    """Upper (parallel) bound on conductivity for matrix + Cu + CNT, MS/m.
    Densities are ASSUMED (not stated in the manuscript); chosen to reproduce the
    manuscript's ~1.8 vol.% Cu and ~1.1 vol.% CNT."""
    m = 100 - cu_wt - cnt_wt
    v = np.array([m / rho_matrix, cu_wt / rho_cu, cnt_wt / rho_cnt])
    f = v / v.sum()
    return float(f @ np.array([sigma_matrix, sigma_cu, sigma_cnt])), f


def run_checks() -> pd.DataFrame:
    t1 = load("table01_nominal_composition"); t6 = load("table06_conductivity_hardness")
    t7 = load("table07_eis_25C"); t8 = load("table08_polarisation_lpr")
    t9 = load("table09_icr_vs_pressure"); t10 = load("table10_temperature_arrhenius")
    rep11 = pd.read_csv(PROCESSED / "reported_table11_factorial.csv").set_index("response")
    C = []

    def add(cid, desc, comp, rep, tol):
        C.append(dict(id=cid, description=desc, computed=float(comp), reported=float(rep), tolerance=tol,
                      passed=bool(abs(comp - rep) <= tol)))

    # ---- batch composition (Sec. 2.1)
    tot = 1200 + 72 + 9
    add("C01", "S4 batch mass total (g)", tot, 1281, 0)
    add("C02", "S4 ADC12 wt.%", 100 * 1200 / tot, 93.68, 0.01)
    add("C03", "S4 Cu wt.%", 100 * 72 / tot, 5.62, 0.01)
    add("C04", "S4 CNT wt.%", 100 * 9 / tot, 0.70, 0.01)
    add("C05", "Table 1 row sums = 100 (max deviation)", float((t1.adc12_wt_pct + t1.cu_added_wt_pct + t1.cnt_wt_pct - 100).abs().max()), 0, 0.011)

    # ---- one-way ANOVA (Sec. 3.3, 3.5, 3.6)
    a = S.oneway_from_summary(t6.sigma_MS_per_m_mean, t6.sigma_MS_per_m_sd, 3)
    add("A01", "sigma one-way F", a["F"], 139.1, 0.1); add("A02", "sigma eta^2", a["eta2"], 0.981, 0.001)
    a = S.oneway_from_summary(t6.hardness_HV_mean, t6.hardness_HV_sd, 3)
    add("A03", "hardness one-way F", a["F"], 38.0, 0.1); add("A04", "hardness eta^2", a["eta2"], 0.934, 0.001)
    a = S.oneway_from_summary(np.log(t7.Rct_kohm_cm2_mean), t7.Rct_kohm_cm2_sd / t7.Rct_kohm_cm2_mean, 3)
    add("A05", "ln Rct one-way F", a["F"], 70.3, 0.1); add("A06", "ln Rct eta^2", a["eta2"], 0.963, 0.001)
    a = S.oneway_from_summary(np.log(t8.Icorr_uA_cm2_mean), t8.Icorr_uA_cm2_sd / t8.Icorr_uA_cm2_mean, 3)
    add("A07", "ln Icorr one-way F", a["F"], 63.2, 0.1); add("A08", "ln Icorr eta^2", a["eta2"], 0.960, 0.001)
    _, h, _ = S.tukey_from_summary(SAMPLES, t6.sigma_MS_per_m_mean.values, t6.sigma_MS_per_m_sd.values, 3)
    add("A09", "Tukey min. resolvable difference, sigma (MS/m)", h, 1.3, 0.1)
    _, h, _ = S.tukey_from_summary(SAMPLES, t6.hardness_HV_mean.values, t6.hardness_HV_sd.values, 3)
    add("A10", "Tukey min. resolvable difference, hardness (HV)", h, 7, 0.5)
    _, h, _ = S.tukey_from_summary(SAMPLES, t7.Rct_kohm_cm2_mean.values, t7.Rct_kohm_cm2_sd.values, 3)
    add("A11", "Tukey min. resolvable difference, Rct (kOhm cm2)", h, 0.51, 0.01)

    # ---- headline percentages
    pct = lambda new, old: 100 * (new / old - 1)
    add("P01", "Rct S4 vs S1, % (EIS series)", pct(4.17, 1.82), 129, 1)
    add("P02", "Icorr S4 vs S1, %", pct(1.92, 4.79), -60, 1)
    add("P03", "ICR(1.4 MPa) S4 vs S1, %", pct(1.11, 3.55), -69, 1)
    add("P04", "sigma S4 vs S1, %", pct(24.6, 16.8), 46, 1)
    add("P05", "hardness S4 vs S1, %", pct(116, 92), 26, 1)
    add("P06", "Rct S4 vs S1 at 25 C (temperature series), %", pct(4.15, 1.85), 124, 1)
    add("P07", "Rct S4 vs S1 at 80 C (temperature series), %", pct(0.50, 0.39), 28, 1)
    add("P08", "DOE ICR target / S4 ICR (x lower)", 10 / 1.11, 9, 0.5)
    add("P09", "sigma in 1e5 S/cm (S4)", 24.6 * 1e6 / 100 / 1e5, 2.46, 0.005)
    add("P10", "S4 %IACS", 100 * 24.6 / IACS_MS_PER_M, 42, 0.5)
    add("P11", "S4 ICR on two-interface basis (mOhm cm2)", float(icr.two_interface_basis(1.11)), 2.2, 0.05)
    add("P12", "Additive sigma prediction 16.8+3.9+2.3", 16.8 + 3.9 + 2.3, 23.0, 1e-9)

    # ---- factorial (Table 11)
    resp = {
        "sigma_MS_per_m": (t6.sigma_MS_per_m_mean, t6.sigma_MS_per_m_sd, False),
        "hardness_HV": (t6.hardness_HV_mean, t6.hardness_HV_sd, False),
        "Rct": (t7.Rct_kohm_cm2_mean, t7.Rct_kohm_cm2_sd, True),
        "Icorr": (t8.Icorr_uA_cm2_mean, t8.Icorr_uA_cm2_sd, True),
        "ICR_1p0MPa": (t9.icr_1p0MPa_mean, t9.icr_1p0MPa_sd, True)}
    pint = []
    for k, (m, s, lg) in resp.items():
        f = S.factorial_2x2(m.values, s.values, 3, log=lg)
        r = rep11.loc[k]
        tol = 0.011 if lg else 0.06
        for eff, pre in [("cu", "cu"), ("cnt", "cnt"), ("interaction", "int")]:
            t2 = tol if not (not lg and k == "hardness_HV") else 0.06
            add(f"F-{k}-{eff}", f"{k}: {eff} effect", f[eff]["estimate"], r[f"{pre}_effect"], t2)
            add(f"F-{k}-{eff}-lo", f"{k}: {eff} CI lower", f[eff]["lo"], r[f"{pre}_lo"], t2 * 1.5)
            add(f"F-{k}-{eff}-hi", f"{k}: {eff} CI upper", f[eff]["hi"], r[f"{pre}_hi"], t2 * 1.5)
        add(f"F-{k}-eta2", f"{k}: eta^2", f["eta2"], r["eta2"], 0.006)
        pint.append(f["interaction"]["p"])
    padj = S.holm(pint)
    for (k, _), p, pa in zip(resp.items(), pint, padj):
        r = rep11.loc[k]
        if str(r["p_int"]).startswith("<"):
            add(f"F-{k}-p<0.001", f"{k}: interaction p < 0.001", float(p < 0.001), 1, 0)
        else:
            add(f"F-{k}-p", f"{k}: interaction p", p, float(r["p_int"]), 0.0015)
            add(f"F-{k}-pHolm", f"{k}: Holm-adjusted p", pa, float(r["p_holm"]), 0.0025)

    # ---- EIS: Brug Ceff, film thickness, LPR cross-validation
    for s_, r in t7.iterrows():
        c = eis.brug_ceff(r.Rs_ohm_cm2_mean, r.Rct_kohm_cm2_mean * 1e3, r["Y0_1e-5_S_sn_cm2"] * 1e-5, r.n_cpe) * 1e6
        add(f"E-Ceff-{s_}", f"Brug Ceff {s_} (uF/cm2)", c, r.Ceff_uF_cm2_mean, 0.06)
    d_nm = 10 * 8.854e-14 / 6.0e-6 * 1e7
    add("E-thick", "Film thickness at C = 6 uF/cm2, eps = 10 (nm)", d_nm, 1.5, 0.05)
    for s_, r in t8.iterrows():
        rct = t7.loc[s_, "Rct_kohm_cm2_mean"] * 1e3
        add(f"L-ratio-{s_}", f"Rct/Rp {s_}", rct / r.Rp_ohm_cm2, {"S1": 1.09, "S2": 1.02, "S3": 1.09, "S4": 1.01}[s_], 0.006)
        B = r.Rp_ohm_cm2 * r.Icorr_uA_cm2_mean * 1e-6 * 1e3  # mV
        add(f"L-B-{s_}", f"Rp*Icorr {s_} (mV)", B, r.RpIcorr_mV, 0.06)
        add(f"L-ba-{s_}", f"implied beta_a {s_} (mV/dec)", pol.implied_beta_a(r.RpIcorr_mV, 120.0), r.beta_a_mV_dec_derived, 0.6)

    # ---- ICR
    P = np.array([0.5, 1.0, 1.5, 2.0])
    for s_, r in t9.iterrows():
        y = np.array([r.icr_0p5MPa_mean, r.icr_1p0MPa_mean, r.icr_1p5MPa_mean, r.icr_2p0MPa_mean])
        add(f"I-1.4-{s_}", f"ICR at 1.4 MPa {s_}", float(icr.quad_loglog_predict(icr.quad_loglog_fit(P, y), 1.4)), r.icr_1p4MPa_interpolated, 0.006)
        add(f"I-b-{s_}", f"secant exponent b {s_}", icr.secant_exponent(P, y), r.secant_exponent_b, 0.006)

    # ---- temperature
    for s_, r in t10.iterrows():
        add(f"T-Ea-{s_}", f"Ea from 25/80 C two-point {s_} (kJ/mol)", eis.ea_two_point(25, r.Rct_25C_kohm_cm2, 80, r.Rct_80C_kohm_cm2), r.Ea_kJ_mol_mean, 0.5)
        add(f"T-ratio-{s_}", f"Rct(80)/Rct(25) {s_}", r.Rct_80C_kohm_cm2 / r.Rct_25C_kohm_cm2, r.Rct_ratio_80_over_25, 0.006)

    # ---- Wiener bound
    wb, f = wiener_bound()
    add("W01", "Cu volume fraction (%)", 100 * f[1], 1.8, 0.1)
    add("W02", "CNT volume fraction (%)", 100 * f[2], 1.1, 0.1)
    add("W03", "Wiener upper bound (MS/m)", wb, 18.5, 0.1)
    add("W04", "Wiener bound vs S1 (%)", 100 * (wb / 16.8 - 1), 10, 1)

    # ---- scorecard + hypervolume ceiling
    _, _, comp = scorecard(t6, t7, t8, t9)
    for s_, v in zip(SAMPLES, [0.49, 0.63, 0.61, 1.00]):
        add(f"SC-{s_}", f"Composite score {s_} (geometric mean of value/best)", comp[s_], v, 0.015)
    add("H01", "Max hypervolume 1.05^4", 1.05 ** 4, 1.22, 0.005)
    add("H02", "1.19 / 1.05^4 (% of maximum)", 100 * 1.19 / 1.05 ** 4, 98, 0.5)
    return pd.DataFrame(C)


def to_markdown(df: pd.DataFrame) -> str:
    n_pass = int(df.passed.sum())
    lines = ["# Manuscript consistency report", "",
             f"{n_pass} of {len(df)} checks reproduce the manuscript values within tolerance.", "",
             "| ID | Check | Computed | Reported | Tol. | Result |", "|---|---|---:|---:|---:|:--:|"]
    for _, r in df.iterrows():
        lines.append(f"| {r['id']} | {r.description} | {r.computed:.4g} | {r.reported:.4g} | {r.tolerance:g} | {'PASS' if r.passed else '**FAIL**'} |")
    return "\n".join(lines) + "\n"
