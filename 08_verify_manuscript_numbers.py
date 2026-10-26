#!/usr/bin/env python3
"""Recompute every derived number quoted in the manuscript from the tabulated data and write a report.

Output: results/tables/manuscript_consistency_report.md (+ .csv).  Exit code 1 if any check fails.
"""
import _bootstrap  # noqa: F401
import sys
from adc12pem import verify
from adc12pem.constants import RES_TABLES

df = verify.run_checks()
RES_TABLES.mkdir(parents=True, exist_ok=True)
df.to_csv(RES_TABLES / "manuscript_consistency_report.csv", index=False)
(RES_TABLES / "manuscript_consistency_report.md").write_text(verify.to_markdown(df))
print(f"{int(df.passed.sum())}/{len(df)} checks pass")
print(df[~df.passed].to_string(index=False) if (~df.passed).any() else "all checks pass")
sys.exit(0 if df.passed.all() else 1)
