"""
style.py -- one visual identity for every figure in the project.

Serif text, thin spines, no chart-junk. The figures are meant to be dropped
straight into the 4-page report, so they are sized in inches to match a
two-column layout (single column 3.4 in, double column 7.0 in) and exported
as both PDF (vector, for LaTeX) and PNG (for slides and GitHub).

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- palette
COL = {
    "ink":      "#1b1b1b",
    "grey":     "#6b6b6b",
    "light":    "#b8b8b8",
    "blue":     "#1f4e79",
    "red":      "#a8322d",
    "green":    "#2f6f4e",
    "orange":   "#c47f17",
    "purple":   "#5b3a70",
    "teal":     "#2a7f7f",
}
CYCLE = [COL["blue"], COL["red"], COL["green"], COL["orange"], COL["purple"], COL["teal"]]

STYLE = {
    "figure.dpi": 160,
    "savefig.dpi": 300,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "STIXGeneral"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 9,
    "axes.labelsize": 9.5,
    "axes.titlesize": 9.5,
    "axes.linewidth": 0.8,
    "axes.edgecolor": COL["ink"],
    "axes.labelcolor": COL["ink"],
    "axes.prop_cycle": plt.cycler(color=CYCLE),
    "axes.grid": True,
    "grid.color": "#dcdcdc",
    "grid.linewidth": 0.5,
    "grid.alpha": 0.8,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.fontsize": 8,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
    "lines.markersize": 3.2,
    "figure.autolayout": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
}


def apply():
    plt.rcParams.update(STYLE)


def panel_label(ax, text, x=0.035, y=0.96):
    """(a), (b), ... inside the top-left corner of the axes (journal style)."""
    ax.text(x, y, text, transform=ax.transAxes, ha="left", va="top",
            fontsize=9.5, fontweight="bold", color=COL["ink"],
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))


def log_ticks(ax, axis="x", values=None):
    """Clean, human-readable ticks on a log axis with a narrow range."""
    values = values or []
    if axis == "x":
        ax.set_xticks(values)
        ax.set_xticklabels([str(v) for v in values])
        ax.minorticks_off()
    else:
        ax.set_yticks(values)
        ax.set_yticklabels([str(v) for v in values])
        ax.minorticks_off()


def finish(fig, path_stem, formats=("pdf", "png")):
    """Save the figure under both formats and close it."""
    out = []
    for ext in formats:
        f = f"{path_stem}.{ext}"
        fig.savefig(f)
        out.append(f)
    plt.close(fig)
    return out
