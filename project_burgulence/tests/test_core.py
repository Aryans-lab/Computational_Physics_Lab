"""
test_core.py -- the checks that must pass before any number in the report is
believed.

Runs as a plain script (`python tests/test_core.py`) and is also compatible
with pytest, since every check is a function named test_*.

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

from fftlib import fft, ifft, dft, PeriodicGrid, fourier_interpolate, dealias_mask
from solvers import spectral, cn_imex, upwind, lax_wendroff, _cyclic_tridiagonal
from exact import cole_hopf, tanh_shock, tanh_width, characteristic_shock_time
from diagnostics import energy_spectrum, shock_width, structure_functions


# ---------------------------------------------------------------- FFT
def test_fft_matches_dft():
    rng = np.random.default_rng(0)
    a = rng.normal(size=128) + 1j * rng.normal(size=128)
    assert np.max(np.abs(fft(a) - dft(a))) < 1e-12


def test_fft_roundtrip():
    rng = np.random.default_rng(1)
    a = rng.normal(size=256)
    assert np.max(np.abs(ifft(fft(a)) - a)) < 1e-13


def test_parseval():
    rng = np.random.default_rng(2)
    u = rng.normal(size=128)
    assert abs(np.mean(u * u) - np.sum(np.abs(fft(u)) ** 2)) < 1e-13


def test_spectral_derivative():
    g = PeriodicGrid(64)
    u = np.sin(g.x) + 0.3 * np.cos(3 * g.x)
    du = g.derivative(u)
    exact = np.cos(g.x) - 0.9 * np.sin(3 * g.x)
    assert np.max(np.abs(du - exact)) < 1e-11


def test_fourier_interpolation_is_exact():
    x = np.arange(32) * 2 * np.pi / 32
    u = np.sin(x)
    fine = fourier_interpolate(u, 256)
    xf = np.arange(256) * 2 * np.pi / 256
    assert np.max(np.abs(fine - np.sin(xf))) < 1e-12


def test_dealias_mask():
    m = dealias_mask(64)
    k = np.fft.fftfreq(64, d=1.0 / 64).astype(int)
    assert m[np.argmax(k)] == False and m[np.argmin(k)] == False   # noqa: E712
    assert m[np.argmax(k == 10)] == True                          # noqa: E712


# ---------------------------------------------------------- exact solution
def test_cole_hopf_reduces_to_ic():
    nu = 0.1
    for n in (32, 128):
        x = np.arange(n) * 2 * np.pi / n
        u = cole_hopf(np.sin(x), x, 1e-4, nu, n_quad=4096)
        assert np.max(np.abs(u - np.sin(x))) < 1e-3


def test_cole_hopf_satisfies_pde():
    """Finite-difference residual of the exact solution, spectrally computed."""
    n, nu, T = 256, 0.1, 0.5
    x = np.arange(n) * 2 * np.pi / n
    u0 = np.sin(x)
    u = cole_hopf(u0, x, T, nu, n_quad=8192)
    dt = 1e-4
    up = cole_hopf(u0, x, T + dt, nu, n_quad=8192)
    um = cole_hopf(u0, x, T - dt, nu, n_quad=8192)
    ut = (up - um) / (2 * dt)
    k = np.fft.fftfreq(n, d=1.0 / n)
    ux = np.real(np.fft.ifft(1j * k * np.fft.fft(u)))
    uxx = np.real(np.fft.ifft(-k * k * np.fft.fft(u)))
    residual = np.max(np.abs(ut + u * ux - nu * uxx))
    assert residual < 1e-6, residual


# ------------------------------------------------------------- solvers
def test_solvers_agree_with_exact():
    nu, T, n = 0.1, 0.5, 128
    x = np.arange(n) * 2 * np.pi / n
    u0 = np.sin(x)
    exact = cole_hopf(u0, x, T, nu, n_quad=8192)
    # the spectral error is 1e-9 in space but ~1e-6 from the Strang splitting
    for fn, tol in [(spectral, 2e-5), (cn_imex, 1e-3),
                    (upwind, 3e-2), (lax_wendroff, 3e-2)]:
        err = np.sqrt(np.mean((fn(u0, nu, T, n=n)["u_final"] - exact) ** 2))
        assert err < tol, (fn.__name__, err)


def test_mass_conservation():
    n, nu = 256, 0.01
    x = np.arange(n) * 2 * np.pi / n
    out = spectral(np.sin(x) + 0.4, nu, 1.0, n=n)
    assert abs(out["mass_drift"]) < 1e-12


def test_cyclic_tridiagonal_matches_dense():
    rng = np.random.default_rng(3)
    n = 40
    a = -0.3 * np.ones(n)
    b = 2.0 * np.ones(n)
    c = -0.3 * np.ones(n)
    alpha = -0.3
    d = rng.normal(size=n)
    x = _cyclic_tridiagonal(a, b, c, d, alpha)
    A = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
    A[0, n - 1] = alpha
    A[n - 1, 0] = alpha
    assert np.max(np.abs(A @ x - d)) < 1e-10


# --------------------------------------------------------- diagnostics
def test_shock_width_of_exact_tanh():
    x = np.linspace(-np.pi, np.pi, 4096, endpoint=False)
    dx = x[1] - x[0]
    for nu in (0.1, 0.02, 0.005):
        w = shock_width(tanh_shock(x, U=1.0, nu=nu), dx, periodic=False)[0][1]
        assert abs(w / tanh_width(nu) - 1.0) < 0.03


def test_energy_spectrum_normalisation():
    n = 256
    x = np.arange(n) * 2 * np.pi / n
    u = np.sin(x) + 0.5 * np.cos(2 * x)
    E = energy_spectrum(u)
    assert abs(np.sum(E) - np.mean(u ** 2)) < 1e-12


def test_structure_function_of_single_mode():
    """For u = A sin x, S_2(r) = A^2 (1 - cos r) exactly."""
    n = 1024
    x = np.arange(n) * 2 * np.pi / n
    A = 1.3
    u = A * np.sin(x)
    r, S = structure_functions([u], 2 * np.pi / n, orders=(2,), r_min=1, r_max=40)
    exact = A ** 2 * (1 - np.cos(r))
    assert np.max(np.abs(S[2] / exact - 1)) < 1e-3


def test_shock_time():
    assert abs(characteristic_shock_time(np.sin, np.linspace(0, 2 * np.pi, 512)) - 1.0) < 1e-3


if __name__ == "__main__":
    tests = [(name, obj) for name, obj in sorted(globals().items())
             if name.startswith("test_") and callable(obj)]
    failures = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as exc:
            failures += 1
            print(f"  FAIL  {name}   {exc}")
    print(f"\n{len(tests) - failures}/{len(tests)} checks passed")
    sys.exit(1 if failures else 0)
