"""Small helpers to read the processed data tables."""
import pandas as pd
from .constants import PROCESSED, RES_TABLES, RES_FIGS


def load(name: str) -> pd.DataFrame:
    """Load data/processed/<name>.csv indexed by sample (when a sample column exists)."""
    df = pd.read_csv(PROCESSED / f"{name}.csv")
    return df.set_index("sample") if "sample" in df.columns else df


def save_table(df: pd.DataFrame, name: str, index: bool = True) -> None:
    RES_TABLES.mkdir(parents=True, exist_ok=True)
    df.to_csv(RES_TABLES / f"{name}.csv", index=index)


def save_fig(fig, name: str, dpi: int = 600) -> None:
    """Save a figure as PNG (raster) and PDF (vector) into results/figures."""
    RES_FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(RES_FIGS / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(RES_FIGS / f"{name}.pdf", bbox_inches="tight")
