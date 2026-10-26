# ADC12–Cu–CNT bipolar plates for PEMFC — data, code and figures

Companion repository for the manuscript

> **Cu–CNT hybrid ADC12 aluminium for PEMFC bipolar plates with factorial assessment of corrosion, contact resistance and conductivity**
> (manuscript file `29_H_ADC12_Cu_CNT_HPDC_revised_v8_clean`)

[![tests](https://github.com/<USER>/<REPO>/actions/workflows/ci.yml/badge.svg)](https://github.com/<USER>/<REPO>/actions)
[![code: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](DATA_LICENSE.md)

## 1. Purpose

The paper compares four high-pressure-die-cast alloys in a 2×2 factorial design (Cu: absent/present × CNT: absent/present):

| code | alloy | added Cu (wt.%) | CNT (wt.%) |
|---|---|---|---|
| S1 | ADC12 | 0 | 0 |
| S2 | ADC12–Cu | 5.62 | 0 |
| S3 | ADC12–CNT | 0 | 0.70 |
| S4 | ADC12–Cu–CNT | 5.62 | 0.70 |

and evaluates conductivity, hardness, EIS (25–80 °C), polarisation/LPR in 0.5 M H₂SO₄ and interfacial contact resistance (ICR) against DOE bipolar-plate targets,
followed by a Gaussian-process / Pareto / ParEGO design loop for the next casting run.

This repository exists so that a reader, reviewer or the authors themselves can

1. **see every number** behind the paper's tables (`data/processed/`),
2. **re-run every analysis** — ANOVA/Tukey, 2×2 factorial contrasts with Holm correction, Rs(Rct‖CPE) fitting, Lin-KK validation, Tafel/LPR, ICR interpolation, Arrhenius analysis, benchmark scorecard, surrogate optimisation (`src/`, `scripts/`),
3. **regenerate the graphs and result tables** (`results/`), and
4. **audit the manuscript's arithmetic** — 130 automatic checks recompute the quoted percentages, F-values, CIs, p-values, capacitances, exponents and activation energies (`results/tables/manuscript_consistency_report.md`).

## 2. What is — and is not — in this repository (please read)

| item | status |
|---|---|
| Tabulated results (Tables 1–10, 13, OCP, microstructure metrics) | ✅ transcribed from the manuscript in `data/processed/` |
| Statistics, circuit/ICR/Arrhenius analysis, scorecard, DOE benchmark | ✅ fully reproducible from the tables |
| Raw instrument files (EIS spectra, polarisation curves, ICR sweeps, XRD, SEM/EDS, indents) | ⚠️ **not available to the repo builder** – templates and loaders are provided in `data/raw/`; the authors must deposit the files (the manuscript's Data-availability statement promises this) |
| Measured data for the 40 "virtual" Latin-hypercube points used for surrogates | ⚠️ **not published in the manuscript** → `data/optimisation/virtual_lhs_40.csv` must be supplied. Until then script 07 runs on a **clearly labelled synthetic demo surface**; its figures carry a red *SYNTHETIC DEMO* banner and its tables a `dataset_origin` column. Demo numbers are *not* the manuscript's Table 14/15 values. |
| Figs 1–9, 11, 14 (photographs, schematics, SEM, XRD, OCP, polarisation curves) | the original images are in `figures/manuscript/` (Fig01–Fig24 as embedded in the docx); they cannot be regenerated without raw data |
| Figs 12, 13 | drawn as **model curves** from the tabulated Table 7 fit parameters, not measured spectra |
| Fig 18 | only 25 °C and 80 °C Rct are tabulated, so lines are two-point Arrhenius interpolations |

## 3. Quick start

```bash
git clone https://github.com/<USER>/<REPO>.git && cd <REPO>
python -m venv .venv && source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt && pip install -e .
python scripts/run_all.py            # everything (a few minutes; add --quick for a shorter optimiser run)
pytest -q                            # 13 unit/regression tests
```

or `make all`, `make verify`, `make test`, `make figures`, `make optimise`. Scripts also work without `pip install -e .` (they add `src/` to `sys.path`).

## 4. Repository layout

```
├── README.md  LICENSE  DATA_LICENSE.md  CITATION.cff  CHANGELOG.md  CONTRIBUTING.md
├── requirements.txt  requirements-tested.txt  environment.yml  pyproject.toml  Makefile
├── .github/workflows/ci.yml            tests + manuscript-number verification on every push
├── data/
│   ├── processed/                      manuscript tables as tidy CSV (+ "reported_*" comparison targets)
│   ├── raw/{eis,polarisation,icr,xrd,sem_eds,hardness_conductivity,ocp}/   templates + READMEs (raw files go here)
│   ├── optimisation/                   anchors are built from processed/; put virtual_lhs_40.csv here
│   ├── casting_records/                melt-by-melt template
│   └── DATA_DICTIONARY.md
├── src/adc12pem/                       importable analysis library
│   ├── stats.py  eis.py  polarisation.py  icr.py  optimisation.py  verify.py  plotting.py  style.py  io.py  constants.py
├── scripts/                            numbered, runnable pipeline (see §5)
├── tests/                              pytest suite (parameter-recovery tests for fitting code, HV vs Monte-Carlo, all 130 manuscript checks)
├── results/
│   ├── tables/                         recomputed tables + consistency report
│   ├── figures/                        PNG (600 dpi) + PDF (vector, editable text)
│   └── optimisation/                   training set used, summary.json (seeds, versions)
├── figures/manuscript/                 Fig01–Fig24 exactly as embedded in the manuscript
└── docs/                               methods summary, figure/table provenance, reproducibility, known issues
```

## 5. Manuscript item → script → output

| manuscript item | script | output |
|---|---|---|
| Table 6, Fig 10 (σ, HV, ANOVA, Tukey) | `01_statistics_and_property_figures.py` | `anova_tukey.csv`, `fig10_*` |
| Table 11, Fig 19 (2×2 factorial, Holm) | `01_…` | `factorial_effects.csv`, `fig19_factorial.*` |
| Table 7, Figs 12–13, Eq. 1–2, Lin-KK | `02_eis_analysis.py` | `eis_table7_consistency.csv`, `fig12_*`, `fig13_*`, `eis_fits_raw.csv` (raw) |
| Table 8, Figs 14–16, Eq. 4 | `03_polarisation_lpr.py` | `eis_lpr_crossvalidation.csv`, `fig15_*`, `fig16_*` |
| Table 9, Fig 17, Eq. 5 | `04_icr_analysis.py` | `icr_analysis.csv`, `fig17_*` |
| Table 10, Fig 18, Eq. 3 | `05_temperature_arrhenius.py` | `temperature_dependence.csv`, `fig18_*` |
| Tables 12–13, Fig 20, Sec. 3.3 bound | `06_benchmark_and_scorecard.py` | `table12_assembled.csv`, `scorecard_normalised.csv`, `conductivity_physical_bound.csv`, `fig20_*` |
| Tables 14–15, Figs 21–24, Eq. 6–7 | `07_optimisation_pipeline.py` | `table14_recomputed.csv`, `table15_recomputed.csv`, `pareto_set.csv`, `hypervolume_traces.csv`, `fig21–24_*` |
| every derived number in the text | `08_verify_manuscript_numbers.py` | `manuscript_consistency_report.md` |

## 6. Headline results reproduced from the tables

| quantity (S4 vs S1) | recomputed | manuscript |
|---|---|---|
| σ | +46 % | +46 % |
| hardness | +26 % | +26 % |
| Rct (25 °C, EIS series) | +129 % | +129 % |
| Icorr | −60 % | −60 % |
| ICR @ 1.4 MPa (per interface) | 1.11 vs 3.55 mΩ cm² (−69 %) | same |
| Rct advantage 25 → 80 °C | +124 % → +28 % | same |
| Cu×CNT interaction ratios (Rct / Icorr / ICR) | 1.34 / 0.65 / 0.66; Holm p = 0.025 / 0.010 / <0.001 | same |
| Rct/Rp (EIS vs LPR) | 1.09, 1.02, 1.09, 1.01 | same |

## 7. Reproducibility notes

* Seeds: Latin-hypercube/Sobol/CV seed 42; ParEGO and random-search seeds 42–47 (as in the manuscript).
* The manuscript used scikit-learn 1.5.2 (`requirements.txt`); the code was developed and tested with the versions in `requirements-tested.txt`. Script 07 prints a warning if the version differs, since GP optimiser details can shift numbers slightly.
* Statistics are computed from mean ± SD (n = 3), which is exact for balanced one-/two-way ANOVA; the Brown-Forsythe/Welch branch of Sec. 2.9 needs raw replicates and is not run.
* See `docs/known_issues.md` for inferred methods (scorecard formula, assumed densities for the conductivity bound, preference vectors for the Chebyshev candidates) and open data questions.

## 8. How to cite, licence, contact

Cite the paper and this repository (`CITATION.cff`; add the Zenodo DOI after the first tagged release — the manuscript says the DOI is inserted at submission).
Code: MIT. Data and figures: CC BY 4.0 (`DATA_LICENSE.md`).
Corresponding author: Giriraj Mannayee, Department of Design and Automation, Vellore Institute of Technology (VIT), Vellore.
