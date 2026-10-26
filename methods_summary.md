# Methods summary (as implemented)

| step | manuscript | implementation |
|---|---|---|
| One-way ANOVA, η², Tukey HSD, compact letters | Sec. 2.9 | `stats.oneway_from_summary`, `tukey_from_summary` (studentized-range distribution) |
| 2×2 factorial contrasts | Sec. 2.9, 3.9 | `stats.factorial_2x2`: Cu = ½[(S2−S1)+(S4−S3)], CNT = ½[(S3−S1)+(S4−S2)], interaction = S4−S2−S3+S1; pooled MSE, 8 df; ratio responses analysed on ln scale with SD(ln) ≈ SD/mean |
| Holm correction | Sec. 2.9 | `stats.holm` over the five interaction p-values |
| Rs(Rct‖CPE), 1/\|Z\| weights | Eq. 1 | `eis.fit_randles_cpe` (log-parametrised, bounded NLLS) |
| Brug effective capacitance | Eq. 2 | `eis.brug_ceff` |
| Lin-KK | Sec. 2.6 | `eis.lin_kk` (M Voigt elements, over-fit criterion c ≥ 0.85, 1 % residual limit) |
| Arrhenius | Eq. 3 | `eis.arrhenius_fit`, `eis.ea_two_point` |
| Tafel / LPR / Stern–Geary | Sec. 2.7, Eq. 4 | `polarisation.tafel_cathodic` (−0.25…−0.12 V vs Ecorr), `lpr_rp` (±5 mV) |
| ICR | Sec. 2.8, Eq. 5 | `icr.quad_loglog_fit` → value at 1.4 MPa; secant/local exponents |
| Surrogates | Sec. 2.10 | RF (150 trees), GBT (200 stages, 3-fold grid), GP (Matérn-5/2 ARD + white noise); repeated 5-fold CV ×2 |
| Pareto | Sec. 2.10 | 2¹⁴ scrambled Sobol designs, GP means, non-dominated sorting, Chebyshev scalarisation |
| ParEGO-EI | Eq. 6–7 | `optimisation.parego`; exact 4-D hypervolume (`hypervolume`), ref 1.05 on min–max-normalised objectives |

Design variables: Cu 0–8 wt.%, CNT 0–1.5 wt.%, stirring 300–700 rpm. Objectives: max σ, max ln Rct, min ln Icorr, min ln ICR(1.0 MPa).
