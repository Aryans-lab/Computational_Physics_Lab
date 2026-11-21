"""
Eleven quick checks of the code.  Runs in about three seconds:

    python tests/test_burgers.py        (plain python, prints PASS/FAIL)
    pytest tests/test_burgers.py        (also works)

Nothing here needs a reference solution from a paper: every check is either a
definition (the grid), an exact property of the equation (the tanh shock, the
sawtooth spectrum) or an obvious limit (advection of a constant is zero).
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import burgers as b   # noqa: E402

NU = 0.01


# ---------------------------------------------------------------- basics ---
def test_grid():
    x, dx = b.make_grid(128)
    assert len(x) == 129
    assert abs(dx - 2 * math.pi / 128) < 1e-14
    assert abs(x[-1] - 2 * math.pi) < 1e-12


def test_initial_sine_boundaries():
    x, _ = b.make_grid(64)
    u = b.initial_sine(x)
    assert abs(u[0]) < 1e-15 and abs(u[-1]) < 1e-12
    assert abs(np.max(u) - 1.0) < 1e-12


def test_advection_of_a_constant_is_zero():
    u = np.full(101, 0.5)
    assert np.max(np.abs(b.advection(u, 0.01))) < 1e-15


def test_diffusion_of_a_straight_line_is_zero():
    x = np.linspace(0, 1, 101)
    u = 2.0 * x - 1.0
    assert np.max(np.abs(b.diffusion(u, x[1] - x[0], NU))) < 1e-10


# ------------------------------------------------------------- the solver ---
def test_pure_diffusion_matches_exact():
    """With the nonlinear term switched off, u = exp(-nu t) sin(x)."""
    nu, T, N = NU, 1.0, 128
    x, dx = b.make_grid(N)
    u = b.initial_sine(x)
    dt = 1e-3
    nsteps = int(T / dt)
    for _ in range(nsteps):
        u = b._implicit_diffusion(u, nu * dt / dx ** 2)
    exact = math.exp(-nu * T) * np.sin(x)
    assert np.max(np.abs(u - exact)) < 2e-3


def test_mean_stays_zero():
    """u0 = sin x is odd about x = pi, so the mean is zero for all time."""
    x, dx = b.make_grid(256)
    run = b.solve("cn", N=256, nu=NU, T=1.0, save_times=[1.0])
    u = run["u"][-1]
    assert abs(np.mean((u[:-1] + u[1:]) / 2.0)) < 1e-12


def test_cn_survives_a_huge_time_step():
    """The explicit scheme would need dt < dx^2/2nu; CN does not care."""
    x, dx = b.make_grid(256)
    dt_limit = dx ** 2 / (2 * NU)
    run = b.solve("cn", N=256, nu=NU, T=0.2, dt=50 * dt_limit, save_times=[0.2])
    assert np.all(np.isfinite(run["u"][-1]))
    assert np.max(np.abs(run["u"][-1])) < 2.0


# ---------------------------------------------------- the exact tanh shock ---
def _tanh_profile(N=4000, U=0.9, nu=NU):
    x, dx = b.make_grid(N)
    return x, dx, -U * np.tanh(U * (x - math.pi) / (2.0 * nu))


def test_tanh_width_is_4nu_over_U():
    x, dx, u = _tanh_profile()
    assert abs(b.shock_width(u, dx) - 4 * NU / 0.9) / (4 * NU / 0.9) < 0.03


def test_tanh_peak_slope_is_U2_over_2nu():
    x, dx, u = _tanh_profile()
    assert abs(b.max_slope(u, dx) - 0.9 ** 2 / (2 * NU)) / (0.9 ** 2 / (2 * NU)) < 0.05


# ------------------------------------------------------------- the spectrum ---
def test_sawtooth_spectrum_slope_is_minus_two():
    """A pure sawtooth has c_k ~ 1/k, so E(k) ~ k^-2.  No PDE involved."""
    n = 4096
    x = np.linspace(0, 2 * math.pi, n, endpoint=False)
    u = np.pi - x                      # straight line with a jump at x = 0
    k, E = b.energy_spectrum(np.append(u, u[0]), 2 * math.pi / n)
    slope, _ = b.fit_power_law(k, E, 4, 40)
    assert abs(slope + 2.0) < 0.15


def test_fit_power_law_recovers_a_known_slope():
    k = np.arange(1, 51, dtype=float)
    E = 3.0 * k ** -1.7
    slope, intercept = b.fit_power_law(k, E, 2, 40)
    assert abs(slope + 1.7) < 1e-6
    assert abs(math.exp(intercept) - 3.0) < 1e-6


# ------------------------------------------------------------------ runner ---
if __name__ == "__main__":
    tests = [(name, fn) for name, fn in sorted(globals().items())
             if name.startswith("test_") and callable(fn)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name[5:]}")
        except AssertionError as exc:
            failed += 1
            print(f"  FAIL  {name[5:]}  ({exc})")
        except Exception as exc:                       # noqa: BLE001
            failed += 1
            print(f"  ERROR {name[5:]}  ({exc})")
    print(f"\n  {len(tests) - failed}/{len(tests)} checks passed")
    sys.exit(1 if failed else 0)
