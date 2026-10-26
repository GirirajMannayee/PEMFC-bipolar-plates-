PY ?= python
.PHONY: all verify test figures optimise quick clean
all: ; $(PY) scripts/run_all.py
quick: ; $(PY) scripts/run_all.py --quick
verify: ; $(PY) scripts/08_verify_manuscript_numbers.py
test: ; $(PY) -m pytest -q
figures: ; for s in 01 02 03 04 05 06; do $(PY) scripts/$$s*.py || exit 1; done
optimise: ; $(PY) scripts/07_optimisation_pipeline.py
clean: ; rm -rf results/figures/* results/tables/* results/optimisation/*
