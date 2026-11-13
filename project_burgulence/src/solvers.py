"""
solvers.py -- five ways to solve the 1-D Navier-Stokes (Burgers) equation

    u_t + u u_x = nu u_xx ,     0 <= x < 2*pi,   periodic

The point of having five solvers is not redundancy: the project is about the
difference between them. A shock is a discontinuity, and each discretisation
disagrees with the continuum in its own way (smearing, ringing, or blow-up).
Those disagreements are the data we analyse.

Solvers
-------
  spectral_rk4 : Fourier pseudo-spectral, 2/3-rule de-aliasing, RK4 in time.
                 Spectral accuracy in space; used as the trusted reference.
  upwind       : first-order upwind advection + explicit diffusion.
                 Adds numerical viscosity -- deliberately.
  lax_wendroff : second-order explicit scheme; dispersive oscillations.
  ftcs         : centred explicit; second-order in space and slightly unstable.
  cn_imex      : Crank-Nicolson diffusion (periodic tridiagonal solve with the
                 Thomas algorithm from mylib.py) + second-order Heun advection.

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fftlib import PeriodicGrid, fft, ifft, dealias_mask  # noqa: E402

try:
    from mylib import thomas_algorithm  # central course library
except Exception:  # pragma: no cover - only for standalone use
    thomas_algorithm = None


# ======================================================================
# helpers
# ======================================================================
def mass(u, dx):
    """Periodic integral of u (should be conserved exactly by all schemes)."""
    return np.sum(u) * dx


def energy(u, dx):
    """Periodic L2 norm squared, (1/2) * integral u^2."""
    return 0.5 * np.sum(u * u) * dx


def _check_periodic(u):
    assert abs(u[0] - u[-1]) < 1e-12, "field is not periodic"


# ======================================================================
# 1. spectral solver  (reference solution)
# ======================================================================
def _spectral_rhs(uh, k, mask, nu, include_diffusion):
    """Right-hand side of the Burgers equation in Fourier space."""
    u = np.real(ifft(uh * mask))
    flux_hat = fft(0.5 * u * u)
    out = -1j * k * flux_hat
    if include_diffusion:
        out = out - nu * k * k * uh
    return out


def spectral(u0, nu, t_end, n=512, cfl=0.25, length=2.0 * np.pi,
             snapshot_times=None, dealias=True, time_method="strang",
             max_steps=400000, dt_fixed=None):
    """Fourier pseudo-spectral solution of the viscous Burgers equation.

    Space: the field is a truncated Fourier series; the nonlinear term is
    formed in physical space as the flux -d/dx(u^2/2) and transformed back.
    The top third of the modes is removed before the product is formed
    (Orszag 2/3 rule), which prevents aliasing from folding high-wavenumber
    energy back into the resolved range.

    Time: two options.
      "strang" (default) -- exact diffusion in Fourier space, half step, then
          a fourth-order Runge-Kutta step of the advection, then the second
          half diffusion step. Second order overall in dt, but the stiff
          diffusion term is handled exactly, so dt is set by the advective
          CFL condition alone. This is what makes long runs affordable.
      "rk4"  -- plain explicit RK4 on the full right-hand side. Fourth order,
          but dt is then limited by nu*k_max^2, i.e. by diffusion.

    Returns dict with the same keys as the finite-difference solvers.
    """
    grid = PeriodicGrid(n, length)
    u = np.array(u0, dtype=float)
    k = grid.k
    mask = dealias_mask(n) if dealias else np.ones(n, dtype=bool)
    snapshot_times = [] if snapshot_times is None else sorted(snapshot_times)
    fields, times = [], []
    t, steps = 0.0, 0
    m0 = mass(u, grid.dx)
    dt = np.nan

    while t < t_end - 1e-12 and steps < max_steps:
        umax = max(np.max(np.abs(u)), 1e-12)
        dt = cfl * grid.dx / umax
        if dt_fixed is not None:
            dt = dt_fixed
        if time_method == "rk4":
            kmax = (n // 2) * (2.0 * np.pi / length)
            dt = min(dt, 0.5 / (nu * kmax * kmax + 1e-30))
        dt = min(dt, t_end - t)

        uh = fft(u)

        if time_method == "strang":
            # exact diffusion for dt/2
            E = np.exp(-nu * k * k * dt * 0.5)
            uh = uh * E
            k1 = _spectral_rhs(uh, k, mask, nu, False)
            k2 = _spectral_rhs(uh + 0.5 * dt * k1, k, mask, nu, False)
            k3 = _spectral_rhs(uh + 0.5 * dt * k2, k, mask, nu, False)
            k4 = _spectral_rhs(uh + dt * k3, k, mask, nu, False)
            uh = (uh + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)) * E

        elif time_method == "rk4":
            k1 = _spectral_rhs(uh, k, mask, nu, True)
            k2 = _spectral_rhs(uh + 0.5 * dt * k1, k, mask, nu, True)
            k3 = _spectral_rhs(uh + 0.5 * dt * k2, k, mask, nu, True)
            k4 = _spectral_rhs(uh + dt * k3, k, mask, nu, True)
            uh = uh + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        else:
            raise ValueError(f"unknown time_method {time_method!r}")

        u = np.real(ifft(uh))
        t += dt
        steps += 1

        while snapshot_times and t >= snapshot_times[0] - 1e-12:
            fields.append(u.copy())
            times.append(snapshot_times.pop(0))

    return {
        "x": grid.x,
        "times": np.array(times),
        "fields": fields,
        "u_final": u,
        "dt": dt,
        "steps": steps,
        "mass_drift": mass(u, grid.dx) - m0,
        "n": n,
        "nu": nu,
    }


# backwards-compatible alias used in the validation script
def spectral_rk4(*a, **kw):
    kw.setdefault("time_method", "rk4")
    return spectral(*a, **kw)


# ======================================================================
# 2-4. explicit finite-difference schemes
# ======================================================================
def _explicit_solve(u0, nu, t_end, n, cfl, scheme, length=2.0 * np.pi,
                    snapshot_times=None, max_steps=2_000_000, dt_fixed=None):
    """Driver shared by the three explicit finite-difference schemes."""
    dx = length / n
    x = np.arange(n) * dx
    u = np.array(u0, dtype=float)
    snapshot_times = [] if snapshot_times is None else sorted(snapshot_times)
    fields, times = [], []
    t, steps = 0.0, 0

    while t < t_end - 1e-12 and steps < max_steps:
        umax = max(np.max(np.abs(u)), 1e-12)
        # advective CFL limit AND the explicit diffusion limit nu*dt/dx^2 < 1/2
        h = min(cfl * dx / umax, 0.45 * dx * dx / max(nu, 1e-30))
        if dt_fixed is not None:
            h = dt_fixed

        if scheme == "upwind":
            # first-order upwind for the advection, centred diffusion
            up = np.roll(u, -1)
            um = np.roll(u, 1)
            adv = np.where(u >= 0.0, u * (u - um) / dx, u * (up - u) / dx)
            diff = nu * (up - 2.0 * u + um) / dx ** 2
            u_new = u + h * (-adv + diff)

        elif scheme == "lax_wendroff":
            # Lax-Wendroff for u_t + f(u)_x = 0, f = u^2/2, plus explicit diffusion
            up = np.roll(u, -1)
            um = np.roll(u, 1)
            r = h / dx
            f = 0.5 * u * u
            fp, fm = 0.5 * up * up, 0.5 * um * um
            u_new = (u - 0.5 * r * (fp - fm)
                     + 0.5 * r ** 2 * (0.5 * (u + up) * (fp - f)
                                       - 0.5 * (u + um) * (f - fm))
                     + nu * h * (up - 2.0 * u + um) / dx ** 2)

        elif scheme == "ftcs":
            up = np.roll(u, -1)
            um = np.roll(u, 1)
            adv = u * (up - um) / (2.0 * dx)
            diff = nu * (up - 2.0 * u + um) / dx ** 2
            u_new = u + h * (-adv + diff)

        else:
            raise ValueError(f"unknown scheme {scheme!r}")

        h = min(h, t_end - t)
        u = u_new
        t += h
        steps += 1

        if not np.all(np.isfinite(u)) or np.max(np.abs(u)) > 1e6:
            return {"x": x, "times": np.array(times), "fields": fields,
                    "u_final": u, "dt": h, "steps": steps,
                    "blew_up": True, "t_blowup": t}

        while snapshot_times and t >= snapshot_times[0] - 1e-12:
            fields.append(u.copy())
            times.append(snapshot_times.pop(0))

    return {"x": x, "times": np.array(times), "fields": fields, "u_final": u,
            "dt": h, "steps": steps, "blew_up": False, "t_blowup": None,
            "n": n, "nu": nu}


def upwind(u0, nu, t_end, n, cfl=0.4, **kw):
    return _explicit_solve(u0, nu, t_end, n, cfl, "upwind", **kw)


def lax_wendroff(u0, nu, t_end, n, cfl=0.4, **kw):
    return _explicit_solve(u0, nu, t_end, n, cfl, "lax_wendroff", **kw)


def ftcs(u0, nu, t_end, n, cfl=0.4, **kw):
    return _explicit_solve(u0, nu, t_end, n, cfl, "ftcs", **kw)


# ======================================================================
# 5. Crank-Nicolson (diffusion) + Heun (advection), periodic
# ======================================================================
def _cyclic_tridiagonal(a, b, c, d, alpha):
    """Solve a periodic tridiagonal system.

        a[i] x[i-1] + b[i] x[i] + c[i] x[i+1] = d[i],  indices modulo n

    with A[0, n-1] = A[n-1, 0] = alpha.  Sherman-Morrison: the cyclic
    system is the tridiagonal system plus a rank-one update, so it needs
    two Thomas sweeps.  The sweeps call mylib.thomas_algorithm.

    Returns None if mylib is unavailable (then numpy is used as a fallback).
    """
    n = len(b)
    if thomas_algorithm is None:
        A = np.zeros((n, n))
        A += np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
        A[0, n - 1] = alpha
        A[n - 1, 0] = alpha
        return np.linalg.solve(A, d)

    gamma = -b[0]
    bb = np.array(b, dtype=float)
    bb[0] = b[0] - gamma
    bb[-1] = b[-1] - alpha * alpha / gamma

    a_sub = list(np.asarray(a[1:], dtype=float))
    c_sup = list(np.asarray(c[:-1], dtype=float))
    bb_l = list(bb)

    x = np.array(thomas_algorithm(a_sub, bb_l, c_sup, list(d)))

    u_rhs = np.zeros(n)
    u_rhs[0] = gamma
    u_rhs[-1] = alpha
    z = np.array(thomas_algorithm(a_sub, bb_l, c_sup, list(u_rhs)))

    fact = (x[0] + alpha * x[-1] / gamma) / (1.0 + z[0] + alpha * z[-1] / gamma)
    return x - fact * z


def cn_imex(u0, nu, t_end, n, cfl=0.5, length=2.0 * np.pi,
            snapshot_times=None, max_steps=200000, dt_fixed=None):
    """Crank-Nicolson in the diffusion, Heun in the advection.

        (1 + 2*alpha) u_i^{n+1} - alpha (u_{i-1} + u_{i+1})^{n+1}
            = u_i^n + alpha (u_{i-1} + u_{i+1})^n - dt * N_i^{n+1/2}

    with alpha = nu*dt/(2*dx^2) and N = u u_x evaluated at the midpoint
    (predictor: explicit Euler; corrector: average of the two slopes).
    Unconditionally stable for the diffusion; the advection still obeys a
    CFL limit, which sets dt.

    The matrix is tridiagonal and periodic, hence one Sherman-Morrison
    correction around two Thomas sweeps.
    """
    dx = length / n
    x = np.arange(n) * dx
    u = np.array(u0, dtype=float)
    snapshot_times = [] if snapshot_times is None else sorted(snapshot_times)
    fields, times = [], []
    t, steps = 0.0, 0

    def advection(f):
        fp = np.roll(f, -1)
        fm = np.roll(f, 1)
        return f * (fp - fm) / (2.0 * dx)

    def diffusion(f):
        return (np.roll(f, -1) - 2.0 * f + np.roll(f, 1)) / dx ** 2

    while t < t_end - 1e-12 and steps < max_steps:
        umax = max(np.max(np.abs(u)), 1e-12)
        dt = cfl * dx / umax if dt_fixed is None else dt_fixed
        dt = min(dt, t_end - t)
        alpha = nu * dt / (2.0 * dx ** 2)

        n_old = advection(u)
        d_old = diffusion(u)

        # predictor (explicit Euler) for N at the half step
        u_pred = u + dt * (-n_old + nu * d_old)
        n_mid = 0.5 * (n_old + advection(u_pred))

        rhs = u + alpha * (np.roll(u, -1) - 2.0 * u + np.roll(u, 1)) - dt * n_mid

        a = np.full(n, -alpha)
        b = np.full(n, 1.0 + 2.0 * alpha)
        c = np.full(n, -alpha)
        u = _cyclic_tridiagonal(a, b, c, rhs, -alpha)

        t += dt
        steps += 1

        while snapshot_times and t >= snapshot_times[0] - 1e-12:
            fields.append(u.copy())
            times.append(snapshot_times.pop(0))

    return {"x": x, "times": np.array(times), "fields": fields, "u_final": u,
            "dt": dt, "steps": steps, "n": n, "nu": nu}
