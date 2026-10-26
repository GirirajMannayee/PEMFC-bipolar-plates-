# Reproducibility checklist
- [x] Every derived number in the text recomputed (`results/tables/manuscript_consistency_report.md`, 130 checks).
- [x] Fixed seeds (42; 42–47) and recorded environment (`results/optimisation/summary.json`).
- [x] Fitting code validated by parameter-recovery tests on synthetic data (`tests/`).
- [ ] Raw data deposited (authors).
- [ ] Authors' virtual-LHS file added → rerun `07` → compare with `data/processed/reported_table14/15`.
- [ ] Tag release, archive on Zenodo, put the DOI in `CITATION.cff` and in the manuscript's Data-availability statement.
