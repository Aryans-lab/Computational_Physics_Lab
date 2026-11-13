"""
diagnostics.py -- the measurements.

Everything the report claims about the flow is produced here:

  * energy spectrum  E(k)  -- the k^-2 law and the dissipative cut-off
  * structure functions S_p(r) -- the bifractal exponents zeta_p
  * shock width      delta -- measured as (drop in u) / (max |u_x|)
  * local logarithmic slopes and least-squares fits with error estimates
  * block averaging for autocorrelated samples

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import numpy as np

from fftlib import fft


# ----------------------------------------------------------------------
# spectra
# ----------------------------------------------------------------------
def energy_spectrum(u, dx=None):
    """One-dimensional energy spectrum E(k), k = 0 ... N/2.

    With u_j = sum_k uhat_k exp(i k x_j) and the normalisation used in
    fftlib.fft (uhat = (1/N) sum_j u_j exp(-i k x_j)), Parseval gives

        <u^2> = sum_k |uhat_k|^2 = E(0) + sum_{k>0} E(k),
        E(0) = |uhat_0|^2,   E(k) = 2 |uhat_k|^2  (k > 0).

    The factor 2 counts the k and -k modes of the real field.
    """
    n = np.asarray(u).size
    uh = fft(u)
    kmax = n // 2
    E = np.zeros(kmax + 1)
    E[0] = np.abs(uh[0]) ** 2
    E[1:kmax] = 2.0 * np.abs(uh[1:kmax]) ** 2
    E[kmax] = np.abs(uh[kmax]) ** 2          # Nyquist mode is not paired
    return E


def compensated_spectrum(E, power=2.0):
    """E(k) * k^power -- a flat curve means E ~ k^-power."""
    k = np.arange(E.size, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(k > 0, E * k ** power, np.nan)


# ----------------------------------------------------------------------
# shock geometry
# ----------------------------------------------------------------------
def shock_width(u, dx, n_shocks=1, periodic=True):
    """Measure shock thickness as delta = (jump in u)/(max |u_x|).

    For the exact profile u = -U tanh(U(x-x0)/(2 nu)) the jump is 2U and the
    maximum slope is U^2/(2 nu), so this definition returns delta = 4 nu / U.

    The jump is measured between the two extrema that bracket the shock, so a
    narrow shock cannot be mistaken for a small one: starting at the steepest
    grid point we walk left and right while the field keeps moving in the
    direction of the shock, and stop at the turning points.

    Returns a list of (x_shock, delta) tuples, strongest shock first.
    """
    n = len(u)
    if periodic:
        du = (np.roll(u, -1) - np.roll(u, 1)) / (2.0 * dx)
    else:
        du = np.gradient(u, dx)
    order = np.argsort(np.abs(du))[::-1]
    picked = []
    for i in order:
        if all(min(abs(i - j), n - abs(i - j)) > max(8, n // 20) for j in picked):
            picked.append(int(i))
        if len(picked) == n_shocks:
            break

    out = []
    for i in picked:
        slope = abs(du[i])
        if slope < 1e-12:
            continue
        sgn = 1.0 if du[i] > 0 else -1.0      # +1: u rises through the shock

        jl = 0
        while jl < n // 2 and sgn * (u[(i - jl - 1) % n] - u[(i - jl) % n]) < 0:
            jl += 1
        jr = 0
        while jr < n // 2 and sgn * (u[(i + jr + 1) % n] - u[(i + jr) % n]) > 0:
            jr += 1

        jump = abs(u[(i - jl) % n] - u[(i + jr) % n])
        out.append((float(i * dx), jump / slope))
    return out


def velocity_increments(u, r_index):
    """Delta u(x) = u(x + r) - u(x) for r = r_index * dx (periodic)."""
    return np.roll(u, -r_index) - u


# ----------------------------------------------------------------------
# structure functions
# ----------------------------------------------------------------------
def structure_functions(fields, dx, orders=(1, 2, 3, 4), r_min=2, r_max=None):
    """Structure functions S_p(r) = < |u(x+r) - u(x)|^p >_x,realisations.

    Parameters
    ----------
    fields : list of 1-D arrays (same length), e.g. one realisation per member
    orders : the p values
    r_min, r_max : separations in grid points

    Returns
    -------
    r (array of physical separations), dict p -> S_p(r)
    """
    n = len(fields[0])
    r_max = r_max if r_max is not None else n // 8
    r_idx = np.arange(r_min, r_max + 1)
    S = {p: np.zeros(len(r_idx)) for p in orders}

    for u in fields:
        for i, m in enumerate(r_idx):
            du = np.abs(velocity_increments(u, int(m)))
            for p in orders:
                S[p][i] += np.mean(du ** p)
    for p in orders:
        S[p] /= len(fields)

    return r_idx * dx, S


def local_slope(r, s, i0=None, i1=None):
    """Slope of log S vs log r between indices i0 and i1 (default: all)."""
    i0 = 1 if i0 is None else i0
    i1 = len(r) - 1 if i1 is None else i1
    lr = np.log(r[i0:i1])
    ls = np.log(s[i0:i1])
    good = np.isfinite(ls) & (s[i0:i1] > 0)
    if good.sum() < 2:
        return np.nan, np.nan
    slope, intercept = np.polyfit(lr[good], ls[good], 1)
    resid = ls[good] - (slope * lr[good] + intercept)
    n = good.sum()
    err = np.sqrt(np.sum(resid ** 2) / max(n - 2, 1) / np.sum((lr[good] - lr[good].mean()) ** 2))
    return slope, err


def fit_line(x, y, mask=None):
    """Least-squares straight line y = m x + c with standard error on m."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if mask is not None:
        x, y = x[mask], y[mask]
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    m, c = np.polyfit(x, y, 1)
    resid = y - (m * x + c)
    n = len(x)
    sxx = np.sum((x - x.mean()) ** 2)
    err = np.sqrt(np.sum(resid ** 2) / max(n - 2, 1) / sxx) if n > 2 else np.nan
    return m, c, err


# ----------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------
def block_average(series, n_blocks=10):
    """Mean and standard error estimated with n_blocks contiguous blocks.

    Naive sigma/sqrt(N) underestimates the error for autocorrelated data;
    block averaging exposes the correlation by letting the variance grow
    with block length until it plateaus.
    """
    series = np.asarray(series, dtype=float)
    n = series.size
    size = max(n // n_blocks, 1)
    blocks = [series[i * size:(i + 1) * size] for i in range(n_blocks) if series[i * size:(i + 1) * size].size]
    means = np.array([b.mean() for b in blocks])
    return means.mean(), means.std(ddof=1) / np.sqrt(len(means))
