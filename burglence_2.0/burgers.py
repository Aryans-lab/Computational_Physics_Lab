"""
Viscous Burgers equation -- the one-dimensional Navier-Stokes equation.

    u_t + u u_x = nu u_xx,      0 <= x <= 2*pi,      u(0,t) = u(2*pi,t) = 0

Everything in the project goes through this file.  There are two solvers and
three small measuring tools, and nothing else: no spectral method, no Strang
splitting, no embedded Runge-Kutta, no Cole-Hopf integral.

Why Dirichlet boundary conditions u = 0 at both ends?
    The initial condition u0(x) = sin(x) is odd about x = 0.  The Burgers
    equation keeps that symmetry (u u_x and u_xx are both odd if u is odd), so
    the exact solution really is zero at x = 0 and x = 2*pi for all time.
    Setting u = 0 there is therefore exact, not an approximation.

House rules of the course
    * the algorithm lives here and leans on mylib.py (Thomas algorithm,
      poly_fit, myrand);
    * numpy is used only for array arithmetic and for the FFT in the spectrum
      experiment;
    * matplotlib is used only in experiments/.

Author: Aryan Bandyopadhyay (2411014)
"""

import math
import os
import sys

import numpy as np

# mylib.py sits one directory above this project folder
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mylib import poly_fit, thomas_algorithm   # noqa: E402


# ---------------------------------------------------------------------------
#  Grid and initial conditions
# ---------------------------------------------------------------------------
def make_grid(N):
    """N intervals on [0, 2 pi] -> N+1 grid points, including both ends."""
    x = np.linspace(0.0, 2.0 * math.pi, N + 1)
    return x, x[1] - x[0]


def initial_sine(x, amplitude=1.0):
    """The classic test: one smooth hump.  u0 = sin(x)."""
    return amplitude * np.sin(x)


# ---------------------------------------------------------------------------
#  Right-hand side:  du/dt = -u u_x + nu u_xx
# ---------------------------------------------------------------------------
def advection(u, dx):
    """-u * du/dx with centred differences.  Zero on the two boundary points."""
    f = np.zeros_like(u)
    f[1:-1] = -u[1:-1] * (u[2:] - u[:-2]) / (2.0 * dx)
    return f


def diffusion(u, dx, nu):
    """nu * d2u/dx2 with centred second differences."""
    f = np.zeros_like(u)
    f[1:-1] = nu * (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dx ** 2
    return f


# ---------------------------------------------------------------------------
#  Solver 1: fully explicit (forward Euler in time)
# ---------------------------------------------------------------------------
def step_explicit(u, dx, dt, nu):
    """
    u_new = u + dt * (-u u_x + nu u_xx)

    One line of physics, one line of code.  Stable only if

        dt  <  dx / max|u|        (advection, CFL)
        dt  <  dx^2 / (2 nu)      (diffusion)

    The second condition is the painful one: halving dx costs a factor 4 in
    the number of time steps.
    """
    return u + dt * (advection(u, dx) + diffusion(u, dx, nu))


# ---------------------------------------------------------------------------
#  Solver 2: Crank-Nicolson for the diffusion term
# ---------------------------------------------------------------------------
def _implicit_diffusion(u, r):
    """
    One Crank-Nicolson step for the diffusion alone:

        (I - r/2 D2) u_new = (I + r/2 D2) u

    D2 is the second-difference operator and r = nu dt / dx^2.  The matrix is
    tridiagonal: sub-diagonal a, diagonal b, super-diagonal c.  One call to
    mylib.thomas_algorithm solves it in O(N) operations.
    """
    n = len(u) - 2                      # interior points only
    d = [0.0] * n
    for i in range(n):
        j = i + 1
        d[i] = u[j] + 0.5 * r * (u[j + 1] - 2.0 * u[j] + u[j - 1])

    a = [-0.5 * r] * (n - 1)
    b = [1.0 + r] * n
    c = [-0.5 * r] * (n - 1)

    u_new = u.copy()
    u_new[1:-1] = np.asarray(thomas_algorithm(a, b, c, d), dtype=float)
    return u_new


def step_cn(u, dx, dt, nu):
    """
    One full step:

        half a step of diffusion  ->  full step of advection  ->  half a step
        of diffusion

    Splitting the step symmetrically like this ("Strang splitting") keeps the
    scheme second order in time.  A plain advection-then-diffusion split would
    only be first order, which is easy to see in the convergence plot.

    Because the diffusion is implicit, the step is stable however small nu*dt
    is -- that is the whole point of the Crank-Nicolson half.
    """
    r = nu * dt / dx ** 2

    u = _implicit_diffusion(u, 0.5 * r)          # first half of the diffusion

    k1 = advection(u, dx)                        # advection, predictor-corrector
    k2 = advection(u + dt * k1, dx)
    u = u + 0.5 * dt * (k1 + k2)

    return _implicit_diffusion(u, 0.5 * r)       # second half of the diffusion


# ---------------------------------------------------------------------------
#  Driver
# ---------------------------------------------------------------------------
def stable_dt(method, u0, dx, nu, safety=0.4):
    """
    Largest time step we are willing to use.

    explicit : the two conditions above, with a safety factor of 0.4
    cn       : only the CFL condition, because diffusion is implicit
    """
    umax = float(np.max(np.abs(u0))) + 1e-12
    dt_cfl = dx / umax
    if method == "explicit":
        dt_diff = dx ** 2 / (2.0 * nu)
        return safety * min(dt_cfl, dt_diff)
    if method == "cn":
        return safety * dt_cfl
    raise ValueError("method must be 'explicit' or 'cn'")


def solve(method="cn", N=512, nu=0.01, T=2.0, u0=None, save_times=(),
          dt=None, verbose=False):
    """
    March the Burgers equation from t = 0 to t = T.

    Returns a dictionary with the grid, the saved profiles and some counters.
    `save_times` is a list of times at which the profile is stored; the nearest
    time step is used, so small overshoots do not matter.
    """
    x, dx = make_grid(N)
    if u0 is None:
        u = initial_sine(x)
    elif callable(u0):
        u = np.asarray(u0(x), dtype=float)
    else:
        u = np.asarray(u0, dtype=float).copy()

    step = step_explicit if method == "explicit" else step_cn
    if dt is None:
        dt = stable_dt(method, u, dx, nu)
    nsteps = int(math.ceil(T / dt))
    dt = T / nsteps                       # land exactly on T

    times = [0.0]
    profiles = [u.copy()]
    targets = sorted(float(t) for t in save_times if t > 0)
    target_i = 0

    for n in range(1, nsteps + 1):
        u = step(u, dx, dt, nu)
        if not np.all(np.isfinite(u)):
            raise FloatingPointError("scheme blew up -- reduce dt")
        t = n * dt
        while target_i < len(targets) and t >= targets[target_i] - 0.5 * dt:
            times.append(t)
            profiles.append(u.copy())
            target_i += 1
        if verbose and n % max(1, nsteps // 10) == 0:
            print(f"    t = {t:5.2f}   max|u| = {np.max(np.abs(u)):.3e}")

    return {"x": x, "dx": dx, "nu": nu, "N": N, "method": method,
            "dt": dt, "nsteps": nsteps, "t": times, "u": profiles}


# ---------------------------------------------------------------------------
#  Measurements
# ---------------------------------------------------------------------------
def max_slope(u, dx):
    """max |du/dx| -- this is what diverges at the breaking time."""
    return float(np.max(np.abs(np.diff(u)) / dx))


def shock_width(u, dx):
    """
    Thickness of the shock, in the usual fluid-dynamics sense:

        delta = (jump in u across the shock) / (max |du/dx|)

    For a tanh profile the jump is 2U and the maximum slope is U^2/(2 nu), so
    delta = 4 nu / U exactly.

    This definition assumes there is ONE shock in the domain, which is true for
    the experiments in this project (u0 = sin x).  We read the jump as the full
    range of u on the grid.
    """
    return float(np.max(u) - np.min(u)) / max_slope(u, dx)


def energy_spectrum(u, dx):
    """
    Fourier spectrum of the profile.

        c_k = (1/N) sum_j u_j exp(-i k x_j),      E(k) = |c_k|^2

    For a train of sharp shocks E(k) falls off like k^-2.  This is the only
    place where the FFT is used, and numpy's FFT is allowed for it.
    """
    u = np.asarray(u, dtype=float)
    n = len(u) - 1                             # drop the repeated endpoint
    v = u[:n] - np.mean(u[:n])
    c = np.fft.rfft(v) / n
    k = 2.0 * math.pi * np.fft.rfftfreq(n, d=dx)
    return k, np.abs(c) ** 2


def fit_power_law(k, E, kmin, kmax):
    """
    Least-squares fit of log E against log k in [kmin, kmax], using
    mylib.poly_fit (degree 1).  Returns (slope, intercept).
    """
    sel = (k >= kmin) & (k <= kmax) & (E > 0)
    data = list(zip(np.log(k[sel]), np.log(E[sel])))
    coeffs, _, _ = poly_fit(data, 1)
    return coeffs[1], coeffs[0]


def l2_difference(u1, u2):
    """Root-mean-square difference of two profiles on the same grid."""
    return float(np.sqrt(np.mean((np.asarray(u1) - np.asarray(u2)) ** 2)))
