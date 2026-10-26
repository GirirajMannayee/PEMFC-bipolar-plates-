"""Shared constants, paths and plotting identifiers."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
PROCESSED = DATA / "processed"
RAW = DATA / "raw"
RESULTS = REPO / "results"
RES_TABLES = RESULTS / "tables"
RES_FIGS = RESULTS / "figures"
RES_OPT = RESULTS / "optimisation"

SAMPLES = ["S1", "S2", "S3", "S4"]
LABELS = {"S1": "ADC12", "S2": "ADC12\u2013Cu", "S3": "ADC12\u2013CNT", "S4": "ADC12\u2013Cu\u2013CNT"}
# Okabe-Ito colour-blind-safe palette
COLOURS = {"S1": "#000000", "S2": "#0072B2", "S3": "#009E73", "S4": "#D55E00"}
MARKERS = {"S1": "o", "S2": "s", "S3": "^", "S4": "D"}

R_GAS = 8.314462618  # J mol-1 K-1
IACS_MS_PER_M = 58.0

# DOE 2020 bipolar-plate targets used in the manuscript
DOE_ICORR_UA_CM2 = 1.0
DOE_ICR_MOHM_CM2 = 10.0
DOE_SIGMA_S_PER_CM = 100.0
DOE_PRESSURE_MPA = 1.4

# Optimisation design space (manuscript Sec. 2.10 / Table 5)
DESIGN_BOUNDS = {"cu_wt": (0.0, 8.0), "cnt_wt": (0.0, 1.5), "rpm": (300.0, 700.0)}
RANDOM_SEEDS = [42, 43, 44, 45, 46, 47]
