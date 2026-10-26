# Data dictionary - `data/processed/`

All values are transcribed from the manuscript (v8_clean). Sample codes: S1 ADC12, S2 ADC12-Cu, S3 ADC12-CNT, S4 ADC12-Cu-CNT. Replicates n = 3 unless noted.
`*_sd` = sample standard deviation. Electrochemical medium 0.5 M H2SO4, potentials vs Ag/AgCl (3 M KCl).

| file | manuscript | contents |
|---|---|---|
| table01_nominal_composition.csv | Table 1 | nominal wt.% (Cu = *added* Cu only) |
| table02_casting_parameters.csv | Table 2 | melt/casting/dispersion parameters |
| table04_electrochemical_conditions.csv | Table 4 | cell and test settings |
| table06_conductivity_hardness.csv | Table 6 | sigma (MS/m, 1e5 S/cm, %IACS), HV, Tukey groups |
| table07_eis_25C.csv | Table 7 | Rs, Rct, Y0, n, Ceff, chi2 (Y0, n, chi2 = means only) |
| table08_polarisation_lpr.csv | Table 8 | Ecorr, Icorr, beta_c, derived beta_a, Rp, Rp*Icorr, passive current, Eb-Ecorr |
| table09_icr_vs_pressure.csv | Table 9 | ICR (mOhm cm2, **per specimen-GDL interface**) at 0.5/1.0/1.5/2.0 MPa; 1.4 MPa interpolation (no SD); secant exponent b |
| table10_temperature_arrhenius.csv | Table 10 | Ea, R2, Rct at 25 and 80 C (separate temperature series) |
| ocp_stable_values.csv | Sec. 3.4 | stable OCP |
| microstructure_metrics_fig8.csv | Sec. 3.1 | Cu-rich particle size/area fraction, CNT cluster area, porosity (Si-particle values are only given as ranges in the text and omitted) |
| doe_targets_and_status.csv | Table 13 | DOE 2020 benchmark status |
| reported_table11_factorial.csv | Table 11 | **reported** factorial effects - used only to check recomputation |
| reported_table14_surrogate_cv.csv | Table 14 | **reported** CV metrics - comparison target for script 07 |
| reported_table15_candidates.csv | Table 15 | **reported** candidate compositions - comparison target |

Units: Rct kOhm cm2 (Table 7/10) - convert to Ohm cm2 before using with Rs/Y0 (the code does this). Y0 column is in units of 1e-5 S s^n cm-2.
