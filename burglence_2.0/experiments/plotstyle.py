# Shared matplotlib style and PNG saver

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = {"blue": "#1f4e79", "red": "#a8322d", "green": "#2e6b4f",
       "orange": "#c47f1a", "grey": "#6b6b6b", "purple": "#5b4a76"}

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
RESDIR = os.path.join(HERE, "..", "results")


def setup():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 8.5,
        "axes.linewidth": 0.7,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "legend.fontsize": 7.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "lines.linewidth": 1.2,
    })


def save(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(FIGDIR, name + ".png"), dpi=220, bbox_inches="tight")
    print("wrote figures/" + name + ".png")
