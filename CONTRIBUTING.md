# Contributing / adding raw data
1. Put raw files in `data/raw/<technique>/` using the column layout of the `TEMPLATE_*.csv` there (never edit `data/processed/` by hand without updating the manuscript).
2. Add the authors' 40 virtual LHS responses as `data/optimisation/virtual_lhs_40.csv`.
3. Run `python scripts/run_all.py` and `pytest -q`; commit regenerated `results/`.
4. Large images (SEM/EDS): Git LFS or a Zenodo deposit; keep a checksum list in the folder README.
