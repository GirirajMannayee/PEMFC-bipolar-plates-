# Raw EIS spectra

Drop one or more CSV files here (any name except `TEMPLATE*`). Required columns (see `TEMPLATE_eis.csv`):

| column | unit | note |
|---|---|---|
| sample | S1..S4 | |
| replicate | 1..3 | |
| temperature_C | 25, 40, 60, 80 | |
| freq_Hz | Hz | 100 kHz - 10 mHz |
| Zre_ohm_cm2 | ohm cm2 | area-normalised |
| Zim_ohm_cm2 | ohm cm2 | **signed** (capacitive = negative) |

`python scripts/02_eis_analysis.py` then runs the Lin-KK test (1 % residual criterion), fits Rs(Rct||CPE) with 1/|Z| weights and writes
`results/tables/eis_fits_raw.csv`; `05_temperature_arrhenius.py` uses it for the 4-point Arrhenius regression.
