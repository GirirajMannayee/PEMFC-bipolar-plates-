"""Journal-style matplotlib defaults (single/double column widths in inches)."""
import matplotlib as mpl

COL1, COL2 = 3.5, 7.2


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.8, "lines.linewidth": 1.2, "lines.markersize": 4.5,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "axes.spines.top": True, "axes.spines.right": True,
        "legend.frameon": False, "figure.dpi": 150,
        "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.default": "regular",
    })


def panel_label(ax, text, x=-0.18, y=1.04):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")
