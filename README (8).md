# Optimisation data

* The four measured anchors (S1-S4) are built on the fly from `data/processed/`.
* **`virtual_lhs_40.csv` (authors' 40 virtual Latin-hypercube points) is not part of the manuscript and must be supplied by the authors.**
  Columns: `cu_wt, cnt_wt, rpm, sigma, ln_Rct, ln_Icorr, ln_ICR1` (or `Rct, Icorr, ICR1` in natural units). See `TEMPLATE_virtual_lhs_40.csv`.
* Without it, `scripts/07_optimisation_pipeline.py` falls back to a clearly labelled SYNTHETIC DEMO surface so the workflow can be exercised;
  those outputs are tagged `SYNTHETIC_DEMO` and are **not** the manuscript's Table 14/15 numbers.
