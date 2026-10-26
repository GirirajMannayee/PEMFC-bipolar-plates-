# Known issues, assumptions and open data questions

## Needed from the authors
1. **Raw data** (EIS, polarisation, ICR, XRD, SEM/EDS, hardness/conductivity readings) — see `data/raw/README.md`.
2. **The 40 virtual Latin-hypercube responses** (and the generator used) — without them Tables 14–15 and Figs 21–24 cannot be reproduced; the repo runs a synthetic demo instead.
3. Intermediate-temperature Rct (40, 60 °C) to reproduce the four-point Arrhenius fits and R² of Table 10.
4. Segmentation settings / field counts behind Fig 8 metrics (manuscript says they will be added to Methods).

## Methods inferred by the repository builder (not stated in the manuscript)
* **Scorecard (Fig 20)**: value ÷ best, composite = geometric mean. Reproduces 0.49/0.63/0.61/1.00 (computed 0.495/0.634/0.615/1.000).
* **Conductivity bound (Sec. 3.3)**: densities ρ(ADC12)=2.74, ρ(Cu)=8.96, ρ(CNT)=1.75 g cm⁻³ are *assumed* (chosen to match the stated 1.8 / 1.1 vol.%); σ(Cu)=58, σ(CNT)=100 MS m⁻¹ as in the text → 18.5 MS m⁻¹.
* **Chebyshev preference vectors** for the three candidates (equal; σ-priority 0.55/0.15/0.15/0.15; corrosion-priority 0.10/0.30/0.50/0.10) are assumptions; the manuscript names the strategies only. ρ = 0.05.
* **Film thickness** 1.5 nm uses ε = 10, ε₀ = 8.854×10⁻¹⁴ F cm⁻¹ and Ceff = 6 µF cm⁻².

## Observations from the data audit (already partly acknowledged in the manuscript)
* ICR of S2, S3, S4 are constant multiples of S1 at all four pressures (ratios 0.72, 0.656, 0.312), hence identical exponents — `icr_analysis.csv`. Raw curves should be checked.
* σ(S4) exceeds the Wiener upper bound by ≈ 6 MS m⁻¹ (+46 % measured vs ≈ +10 % possible).
* Two Rct(25 °C) series exist (EIS 1.82/4.17 vs temperature series 1.85/4.15 kΩ cm²); the text uses +129 % (EIS) and +124 % (temperature series) in different places.
* Icorr uncertainty from the Stern–Geary constant is ±30 % (`eis_lpr_crossvalidation.csv`), so the 1.9 µA cm⁻² vs DOE 1 µA cm⁻² gap is not resolved in absolute terms.
* The manuscript's Table 14 R² for ln Rct/ln Icorr (0.51/0.54) are modest; the abstract quotes only the two high values (σ, ICR).
* Fig 10 in the manuscript shows individual specimens; only mean ± SD are tabulated, so regenerated bars show SD only.
