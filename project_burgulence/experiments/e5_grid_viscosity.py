"""
E5 -- How much of the viscosity is the grid?

Every finite-difference scheme adds dissipation of its own. For first-order
upwind the truncation error of the advection term is a diffusion term with
coefficient

        nu_num = (U dx / 2) (1 - C),        C = U dt/dx   (Courant number),

so on a coarse grid the scheme can be far more viscous than the fluid. The
same is not true of a centred scheme, which is non-dissipative but dispersive.

This experiment measures nu_num directly -- from the damping of a single
Fourier mode of a pure advection problem, where the exact answer is known --
and then compares it with the physical viscosity the user asked for.

Run:  python experiments/e5_grid_viscosity.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import style
from common import FIGDIR, RESDIR, save_json, save_csv, Tee
from diagnostics import fit_line

style.apply()
import matplotlib.pyplot as plt


def mode_amplitude(u, k=1):
    """Amplitude of the k-th Fourier mode: u(x) ~ 2|a_k| cos(kx + phi)."""
    x = np.arange(len(u)) * 2.0 * np.pi / len(u)
    return 2.0 * abs(np.mean(u * np.exp(-1j * k * x)))


def linear_advection(n, u0, U=1.0, cfl=0.4, t_end=1.0, scheme="upwind"):
    """Solve u_t + U u_x = 0 with constant U (a pure test of the scheme).

    upwind : first-order, dissipative, |G| < 1  ->  numerical viscosity
    central: Crank-Nicolson in time, centred in space, |G| = 1 -> no damping,
             only a phase (dispersion) error.
    """
    dx = 2.0 * np.pi / n
    u = np.array(u0, dtype=float)
    t = 0.0
    while t < t_end - 1e-12:
        dt = min(cfl * dx / U, t_end - t)
        if scheme == "upwind":
            up, um = np.roll(u, -1), np.roll(u, 1)
            # u_i^{n+1} = u_i - C (u_i - u_{i-1}),  U > 0
            u = u - (U * dt / dx) * (u - um)
        else:
            up, um = np.roll(u, -1), np.roll(u, 1)
            c = U * dt / (4.0 * dx)
            # (1 + c D0) u^{n+1} = (1 - c D0) u^n  with D0 = centred difference
            rhs = u - c * (up - um)
            # cyclic tridiagonal solve, alpha = c
            from solvers import _cyclic_tridiagonal
            a = np.full(n, c)
            b = np.ones(n)
            u = _cyclic_tridiagonal(a, b, a, rhs, c)
        t += dt
    return u, dx


def main():
    tee = Tee(os.path.join(RESDIR, "e5_grid_viscosity.txt"))
    sys.stdout = tee
    print("=" * 78)
    print("E5  Numerical viscosity: how much of the dissipation is the grid?")
    print("=" * 78)

    # ------------------------------------------------ (a) mode damping
    print("\nPure advection u_t + u_x = 0, u0 = sin x, one period, C = 0.4:")
    print("  N      nu_num (measured)  (U dx/2)(1-C)   ratio")
    rows = []
    for n in [32, 64, 128, 256, 512, 1024]:
        x = np.arange(n) * 2.0 * np.pi / n
        u0 = np.sin(x)
        u_up, dx = linear_advection(n, u0, scheme="upwind")
        amp = mode_amplitude(u_up)
        nu_num = -np.log(max(amp, 1e-16)) / 1.0
        theory = 0.5 * dx * (1.0 - 0.4)
        rows.append((n, dx, 0.4, nu_num, theory, 0.0))
        print(f"  {n:5d}   {nu_num:10.4e}    {theory:10.4e}   {nu_num/theory:5.2f}")

    m, c, e = fit_line(np.log([r[1] for r in rows]), np.log([r[3] for r in rows]))
    print(f"  fitted nu_num ~ dx^{m:.3f}   (theory: 1)")

    # Courant-number dependence at fixed dx
    print("\nSame grid (N=256), varying Courant number (upwind):")
    cfl_rows = []
    n = 256
    dx = 2.0 * np.pi / n
    for C in [0.1, 0.2, 0.4, 0.6, 0.8]:
        x = np.arange(n) * 2.0 * np.pi / n
        u_up, _ = linear_advection(n, np.sin(x), cfl=C, scheme="upwind")
        amp = mode_amplitude(u_up)
        nu_num = -np.log(max(amp, 1e-16))
        cfl_rows.append((C, nu_num, 0.5 * dx * (1.0 - C)))
        print(f"  C={C:4.2f}   nu_num = {nu_num:.4e}   theory = {0.5*dx*(1-C):.4e}   "
              f"ratio {nu_num/(0.5*dx*(1-C)):.2f}")

    # ------------------------------------------------ (b) grid vs fluid
    nu_phys = 2e-3
    print(f"\nGrid viscosity relative to a physical nu = {nu_phys}:")
    ratio_rows = []
    for n in [128, 256, 512, 1024, 2048, 4096]:
        dx = 2.0 * np.pi / n
        nu_num = 0.5 * dx * (1.0 - 0.4)
        delta_phys = 4.0 * nu_phys
        ratio_rows.append((n, dx, nu_num, nu_num / nu_phys, delta_phys / dx))
        print(f"  N={n:5d}  dx={dx:.4e}  nu_num={nu_num:.4e}  "
              f"nu_num/nu = {nu_num/nu_phys:6.1f}   points inside 4nu/U = {delta_phys/dx:5.1f}")

    # ------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.55))

    ax = axes[0]
    dxs = np.array([r[1] for r in rows])
    nus = np.array([r[3] for r in rows])
    ax.loglog(dxs, nus, "o", color=style.COL["blue"], ms=4, label="upwind (measured)")
    ax.loglog(dxs, 0.5 * dxs * 0.6, "--", color=style.COL["grey"], lw=0.9,
              label=r"$\frac{1}{2}U\Delta x(1-C)$")
    ax.set_xlabel(r"$\Delta x$")
    ax.set_ylabel(r"$\nu_{\rm num}$")
    ax.set_ylim(1e-4, 1)
    ax.set_title(r"(a)  numerical viscosity", loc="left", fontsize=9)
    ax.legend(fontsize=7, loc="upper right")

    ax = axes[1]
    ns = np.array([r[0] for r in ratio_rows], dtype=float)
    rr = np.array([r[3] for r in ratio_rows])
    ax.loglog(ns, rr, "o-", color=style.COL["red"], ms=4, lw=1.0)
    ax.axhline(1.0, color=style.COL["grey"], lw=0.9, ls="--")
    ax.axhline(0.25, color=style.COL["blue"], lw=0.9, ls=":")
    ax.text(140, 1.15, r"grid = fluid", fontsize=7, color=style.COL["grey"])
    ax.text(140, 0.28, r"25% criterion", fontsize=7, color=style.COL["blue"])
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\nu_{\rm num}/\nu$")
    ax.set_xlim(110, 5e3)
    ax.set_ylim(0.05, 1e2)
    ax.set_title(r"(b)  with $\nu=2\times10^{-3}$", loc="left", fontsize=9)
    style.log_ticks(ax, "x", [128, 512, 2048])
    ax.set_xticklabels(["128", "512", "2048"])

    fig.tight_layout(w_pad=1.5)
    style.finish(fig, os.path.join(FIGDIR, "fig06_grid_viscosity"))
    print("\nwrote figures/fig06_grid_viscosity.{pdf,png}")

    save_json("e5_grid_viscosity.json", {
        "nu_vs_dx": [{"N": r[0], "dx": r[1], "C": r[2], "nu_num_upwind": r[3],
                      "theory": r[4], "nu_num_central": r[5]} for r in rows],
        "courant_scan": [{"C": r[0], "nu_num": r[1], "theory": r[2]} for r in cfl_rows],
        "grid_to_fluid": [{"N": r[0], "dx": r[1], "nu_num": r[2], "ratio": r[3],
                           "points_in_shock": r[4]} for r in ratio_rows],
        "nu_num_slope": float(m), "nu_phys": nu_phys,
    })
    save_csv("e5_grid_viscosity.csv", ["N", "dx", "nu_num", "nu_num_over_nu", "points_in_shock"],
             ratio_rows)

    tee.close()
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
