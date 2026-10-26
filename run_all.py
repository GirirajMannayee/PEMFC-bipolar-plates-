#!/usr/bin/env python3
"""Run the whole pipeline in order.  ``python scripts/run_all.py [--quick]``"""
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
steps = ["08_verify_manuscript_numbers.py", "01_statistics_and_property_figures.py", "02_eis_analysis.py", "03_polarisation_lpr.py",
         "04_icr_analysis.py", "05_temperature_arrhenius.py", "06_benchmark_and_scorecard.py", "07_optimisation_pipeline.py"]
extra = ["--quick"] if "--quick" in sys.argv else []
for s in steps:
    print(f"\n=== {s} ===", flush=True)
    r = subprocess.run([sys.executable, str(here / s)] + (extra if s.startswith("07") else []))
    if r.returncode:
        sys.exit(f"{s} failed")
print("\nAll steps finished. Outputs: results/tables, results/figures, results/optimisation")
