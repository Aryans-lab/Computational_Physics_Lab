"""
common.py -- shared paths, initial conditions and small utilities.

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "figures")
RESDIR = os.path.join(ROOT, "results")
DATADIR = os.path.join(ROOT, "data")

for d in (FIGDIR, RESDIR, DATADIR):
    os.makedirs(d, exist_ok=True)

sys.path.insert(0, os.path.join(ROOT, "src"))


def sine_ic(n, length=2.0 * np.pi):
    """u0 = sin x on a periodic grid: the canonical shock-formation problem."""
    x = np.arange(n) * length / n
    return x, np.sin(x)


def random_ic(n, seed=1, n_modes=12, rms=1.0, length=2.0 * np.pi):
    """Band-limited random initial condition with random phases.

    u0(x) = sum_{k=1}^{n_modes} a_k cos(k x + phi_k),   a_k = k^-1,
    normalised to a prescribed rms. Phases are drawn from the project's own
    linear congruential generator (mylib.myrand), so the initial ensemble is
    reproducible and consistent with the course's random-number work.
    """
    from mylib import myrand

    myrand(seed)
    x = np.arange(n) * length / n
    u = np.zeros(n)
    for k in range(1, n_modes + 1):
        phi = 2.0 * np.pi * myrand()
        u += np.cos(k * x * 2.0 * np.pi / length + phi) / k
    u -= u.mean()
    u *= rms / np.sqrt(np.mean(u * u))
    return x, u


def wave_numbers(n):
    return np.arange(n // 2 + 1)


def save_json(name, payload):
    path = os.path.join(RESDIR, name)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    return path


def save_csv(name, header, rows):
    path = os.path.join(RESDIR, name)
    with open(path, "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(repr(float(v)) if isinstance(v, (int, float, np.floating))
                             else str(v) for v in r) + "\n")
    return path


class Tee:
    """Duplicate stdout into a results file, so every number in the report
    can be traced back to a run."""

    def __init__(self, path):
        self.path = path
        self.file = open(path, "w")

    def write(self, text):
        sys.__stdout__.write(text)
        self.file.write(text)

    def flush(self):
        sys.__stdout__.flush()
        self.file.flush()

    def close(self):
        self.file.close()
