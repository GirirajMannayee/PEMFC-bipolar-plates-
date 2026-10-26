"""Reusable plotting helpers (journal-style, colour-blind-safe)."""
import numpy as np
from .constants import COLOURS, LABELS, SAMPLES


def sample_bars(ax, means, sds, letters=None, ylabel="", ylim=None, fmt="{:.1f}"):
    """Bar chart of S1-S4 (mean +/- SD) with optional Tukey letters."""
    x = np.arange(len(SAMPLES))
    ax.bar(x, [means[s] for s in SAMPLES], yerr=[sds[s] for s in SAMPLES],
           color=[COLOURS[s] for s in SAMPLES], alpha=0.85, edgecolor="black", linewidth=0.6,
           capsize=3, error_kw=dict(elinewidth=0.8, capthick=0.8))
    ax.set_xticks(x)
    short = {"S1": "ADC12", "S2": "+Cu", "S3": "+CNT", "S4": "+Cu+CNT"}
    ax.set_xticklabels([f"{s}\n{short[s]}" for s in SAMPLES], fontsize=6.5)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    top = ax.get_ylim()[1]
    for i, s in enumerate(SAMPLES):
        y = means[s] + sds[s] + 0.015 * top
        ax.text(i, y, (letters[s] if letters else "") , ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.margins(x=0.08)


def watermark(ax, text="SYNTHETIC DEMO DATA \u2013 not manuscript results", y=0.97):
    ax.text(0.5, y, text, transform=ax.transAxes, ha="center", va="top", fontsize=6.5,
            color="#B00020", alpha=0.85, fontweight="bold")
