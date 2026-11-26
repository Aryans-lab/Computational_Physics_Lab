# 1D viscous Burgers equation: u_t + u * u_x = nu * u_xx on [0, 2*pi]
# Boundary conditions u(0, t) = u(2*pi, t) = 0 hold exactly because sin(x) stays odd
# Author: Aryan Bandyopadhyay (2411014)

import math
import os
import sys

import numpy as np

# mylib.py lives one folder up
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mylib import poly_fit, thomas_algorithm  # noqa: E402

def make_grid(N):
    # N intervals on [0, 2*pi] gives N+1 points
    x = np.linspace(0.0, 2.0 * math.pi, N + 1)
    return x, x[1] - x[0]

def initial_sine(x, amplitude=1.0):
    return amplitude * np.sin(x)

def advection(u, dx):
    # Central difference for nonlinear term, keeping boundaries at 0
    f = np.zeros_like(u)
    f[1:-1] = -u[1:-1] * (u[2:] - u[:-2]) / (2.0 * dx)
    return f

def diffusion(u, dx, nu):
    # Central second difference for viscosity term, keeping boundaries at 0
    f = np.zeros_like(u)
    f[1:-1] = nu * (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dx ** 2
    return f

def step_explicit(u, dx, dt, nu):
    # Forward Euler step (needs small dt for both advection and diffusion)
    return u + dt * (advection(u, dx) + diffusion(u, dx, nu))

def _implicit_diffusion(u, r):
    # Crank Nicolson diffusion step solved as a tridiagonal system via Thomas algorithm
    # Directly making the matrix
    n = len(u) - 2
    d = [u[j] + 0.5 * r * (u[j + 1] - 2.0 * u[j] + u[j - 1]) for j in range(1, n + 1)]
    a = [-0.5 * r] * (n - 1)
    b = [1.0 + r] * n
    c = [-0.5 * r] * (n - 1)

    u_new = u.copy()
    u_new[1:-1] = np.asarray(thomas_algorithm(a, b, c, d), dtype=float)
    return u_new

def step_cn(u, dx, dt, nu):
    # Symmetric Strang splitting: half diffusion, full Heun advection, half diffusion
    r = nu * dt / dx ** 2
    u = _implicit_diffusion(u, 0.5 * r)
    k1 = advection(u, dx)
    k2 = advection(u + dt * k1, dx)
    u = u + 0.5 * dt * (k1 + k2)
    return _implicit_diffusion(u, 0.5 * r)

def stable_dt(method, u0, dx, nu, safety=0.4):
    umax = float(np.max(np.abs(u0))) + 1e-12
    dt_cfl = dx / umax
    if method == "explicit":
        return safety * min(dt_cfl, dx ** 2 / (2.0 * nu))
    if method == "cn":
        # Diffusion is implicit, so only the CFL limit applies
        return safety * dt_cfl
    raise ValueError("method must be 'explicit' or 'cn'")

def solve(method="cn", N=512, nu=0.01, T=2.0, u0=None, save_times=(), dt=None, verbose=False):
    x, dx = make_grid(N)
    if u0 is None:
        u = initial_sine(x)
    elif callable(u0):
        u = np.asarray(u0(x), dtype=float)
    else:
        u = np.asarray(u0, dtype=float).copy()

    step = step_explicit if method == "explicit" else step_cn
    dt = stable_dt(method, u, dx, nu) if dt is None else dt

    # Adjust dt so an integer number of steps lands exactly at T
    nsteps = int(math.ceil(T / dt))
    dt = T / nsteps

    times, profiles = [0.0], [u.copy()]
    targets = sorted(float(t) for t in save_times if t > 0)
    ti = 0

    for n in range(1, nsteps + 1):
        u = step(u, dx, dt, nu)
        if not np.all(np.isfinite(u)):
            raise FloatingPointError("scheme blew up, reduce dt")

        t = n * dt
        while ti < len(targets) and t >= targets[ti] - 0.5 * dt:
            times.append(t)
            profiles.append(u.copy())
            ti += 1

        if verbose and n % max(1, nsteps // 10) == 0:
            print(f"    t = {t:5.2f}   max|u| = {np.max(np.abs(u)):.3e}")

    return {"x": x, "dx": dx, "nu": nu, "N": N, "method": method,
            "dt": dt, "nsteps": nsteps, "t": times, "u": profiles}

def max_slope(u, dx):
    return float(np.max(np.abs(np.diff(u)) / dx))

def shock_width(u, dx):
    # Fluid thickness definition: total jump divided by steepest slope (4*nu/U for tanh)
    return float(np.max(u) - np.min(u)) / max_slope(u, dx)

def energy_spectrum(u, dx):
    # Real FFT after dropping the repeated endpoint and subtracting the mean
    u = np.asarray(u, dtype=float)
    n = len(u) - 1
    v = u[:n] - np.mean(u[:n])
    c = np.fft.rfft(v) / n
    k = 2.0 * math.pi * np.fft.rfftfreq(n, d=dx)
    return k, np.abs(c) ** 2

def fit_power_law(k, E, kmin, kmax):
    # Linear fit on log(E) vs log(k) using mylib.poly_fit; returns (slope, intercept)
    sel = (k >= kmin) & (k <= kmax) & (E > 0)
    data = list(zip(np.log(k[sel]), np.log(E[sel])))
    coeffs, _, _ = poly_fit(data, 1)
    return coeffs[1], coeffs[0]

def l2_difference(u1, u2):
    return float(np.sqrt(np.mean((np.asarray(u1) - np.asarray(u2)) ** 2)))
