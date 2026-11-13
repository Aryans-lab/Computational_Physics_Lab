"""
E2 -- The shock: when does the wave break, and how thin is the shock?

Start from u0 = sin x. In the inviscid limit every point of the wave travels at
its own speed u(x), so the crest catches up with the trough and the profile
steepens until characteristics cross at

        t_s = -1 / min_x u0'(x) = 1        (at x = pi).

Before t_s the solution is the implicit characteristic solution, which we can
plot exactly with no time stepping at all. After t_s a shock exists whose
thickness is set by the only length scale left in the problem, nu/U:

        delta = 4 nu / U        (amplitude to maximum-slope width).

Two figures:
  fig02 -- space-time diagram and the steepening profiles (the physics)
  fig03 -- divergence of the strain rate at t_s, and the measured delta(nu)
           scaling (the quantitative test)

Run:  python experiments/e2_shock_formation.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import style
from common import FIGDIR, RESDIR, save_json, save_csv, Tee
from exact import characteristic_solution, tanh_shock, tanh_width
from fftlib import PeriodicGrid
from solvers import spectral
from diagnostics import shock_width

style.apply()
import matplotlib.pyplot as plt

NU_MAIN = 2e-3
N_MAIN = 4096
T_SNAPS = [0.0, 0.4, 0.8, 1.0, 1.4, 2.2]
NU_LIST = [0.02, 0.01, 0.005, 0.002, 0.0015]


def strain(u, grid):
    """max |u_x| measured spectrally."""
    return float(np.max(np.abs(grid.derivative(u))))


def main():
    tee = Tee(os.path.join(RESDIR, "e2_shock_formation.txt"))
    sys.stdout = tee
    print("=" * 78)
    print("E2  Shock formation from u0 = sin x")
    print("=" * 78)

    n, nu = N_MAIN, NU_MAIN
    x = np.arange(n) * 2.0 * np.pi / n
    u0 = np.sin(x)
    grid = PeriodicGrid(n, 2.0 * np.pi)

    dense = np.linspace(0.0, 2.2, 221)
    snaps = sorted(set([round(t, 6) for t in dense] + T_SNAPS))
    print(f"spectral run: N={n}, nu={nu}, snapshots={len(snaps)}")
    out = spectral(u0, nu, 2.2, n=n, cfl=0.4, snapshot_times=snaps)
    print(f"steps = {out['steps']}, dt = {out['dt']:.2e}, mass drift = {out['mass_drift']:.2e}")

    times = np.array(out["times"])
    fields = out["fields"]

    # ---------------------------------------------------- steepening curves
    strain_t = np.array([strain(f, grid) for f in fields])
    print("\n  t      max|u_x|     inviscid 1/(t_s-t)")
    for target in [0.2, 0.5, 0.8, 0.95]:
        i = int(np.argmin(np.abs(times - target)))
        inv = 1.0 / max(1.0 - times[i], 1e-12)
        print(f"  {times[i]:.2f}   {strain_t[i]:9.3f}   {inv:10.3f}")

    # saturation value predicted by the tanh profile: U^2/(2 nu) with U ~ 1
    strain_max_pred = 1.0 / (2.0 * nu)
    print(f"\n  viscous saturation, predicted U^2/(2 nu) = {strain_max_pred:.1f}, "
          f"measured {strain_t.max():.1f}")

    # ------------------------------------------------------- width vs nu
    width_rows = []
    for nu_i in NU_LIST:
        n_i = 4096 if nu_i >= 2e-3 else 8192
        x_i = np.arange(n_i) * 2.0 * np.pi / n_i
        out_i = spectral(np.sin(x_i), nu_i, 2.0, n=n_i, cfl=0.4)
        w = shock_width(out_i["u_final"], 2.0 * np.pi / n_i, n_shocks=1)[0][1]
        pred = tanh_width(nu_i, U=1.0)
        width_rows.append((nu_i, w, pred, n_i))
        print(f"  nu={nu_i:7.4f}  N={n_i:5d}  delta={w:.4f}   4 nu/U={pred:.4f}   "
              f"ratio={w/pred:.2f}")

    nus = np.array([r[0] for r in width_rows])
    wid = np.array([r[1] for r in width_rows])
    slope, intercept = np.polyfit(np.log(nus), np.log(wid), 1)
    print(f"  fitted delta ~ nu^{slope:.3f}   (exact: 1)")

    # ------------------------------------------------------------- figures
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1.15, 1]})

    ax = axes[0]
    X, T = np.meshgrid(x, times)
    U = np.array(fields)
    pc = ax.pcolormesh(X, T, U, cmap="RdBu_r", shading="auto", vmin=-1, vmax=1,
                       rasterized=True, zorder=2)
    ax.grid(False)
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    ax.axhline(1.0, color="k", lw=0.7, ls="--")
    ax.text(0.12, 1.05, r"$t_s=1$", fontsize=7.5)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$t$")
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(0, 2.2)
    ax.set_xticks([0, np.pi, 2 * np.pi])
    ax.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax.set_title(r"(a)  steepening, $\nu=2\times10^{-3}$", loc="left", fontsize=9)
    cb = fig.colorbar(pc, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(r"$u(x,t)$", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    ax = axes[1]
    # the post-shock profile, zoomed, against the exact tanh solution
    i2 = int(np.argmin(np.abs(times - 2.2)))
    i08 = int(np.argmin(np.abs(times - 0.8)))
    win = 0.045
    m = np.abs(x - np.pi) < win
    ax.plot((x[m] - np.pi) * 1e2, fields[i08][m], color=style.COL["blue"], lw=1.1,
            label=r"$t=0.8$ (pre-shock)")
    ax.plot((x[m] - np.pi) * 1e2, fields[i2][m], "o", ms=2.0, color=style.COL["red"],
            label=r"$t=2.2$ (solver)", markevery=3)
    ax.plot((x[m] - np.pi) * 1e2, tanh_shock(x[m], U=0.935, nu=nu, x0=np.pi),
            color="k", lw=0.9, ls="--", label=r"exact $\tanh$")
    ax.set_xlabel(r"$(x-\pi)\times 10^{2}$")
    ax.set_ylabel(r"$u$")
    ax.set_title("(b)  the asymptotic profile", loc="left", fontsize=9)
    ax.legend(fontsize=6.8, loc="lower left", frameon=True, framealpha=0.92,
              edgecolor="none", handlelength=1.4)

    fig.tight_layout(w_pad=1.4)
    style.finish(fig, os.path.join(FIGDIR, "fig02_shock_formation"))

    # ------------------------------------------------------------ fig03
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.4))

    ax = axes[0]
    ax.semilogy(times, strain_t, "o", ms=1.8, color=style.COL["blue"], label="measured")
    tt = np.linspace(0.02, 0.985, 400)
    ax.semilogy(tt, 1.0 / (1.0 - tt), "-", color=style.COL["grey"], lw=0.9,
                label=r"inviscid $1/(t_s-t)$")
    ax.axhline(strain_max_pred, color=style.COL["red"], ls="--", lw=0.8)
    ax.text(0.15, strain_max_pred * 1.35, r"$U^2/2\nu$", color=style.COL["red"], fontsize=7.5)
    ax.set_ylim(0.6, 900)
    ax.axvline(1.0, color="k", lw=0.6, ls=":")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\max|u_x|$")
    ax.set_xlim(0, 2.3)
    ax.set_title(r"(a)  strain diverges, then saturates", loc="left", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")

    ax = axes[1]
    ax.loglog(nus, wid, "o", color=style.COL["green"], ms=4, label="measured")
    ax.loglog(nus, np.exp(intercept) * nus ** slope, "-", color=style.COL["green"],
              lw=0.8, label=rf"fit: $\delta\propto\nu^{{{slope:.2f}}}$")
    ax.loglog(nus, 4 * nus, "--", color=style.COL["grey"], lw=0.9, label=r"$4\nu/U$")
    ax.set_xlabel(r"$\nu$")
    ax.set_ylabel(r"shock width $\delta$")
    style.log_ticks(ax, "x", [0.002, 0.005, 0.01, 0.02])
    ax.set_xticklabels(["0.002", "0.005", "0.01", "0.02"])
    ax.set_yticks([0.006, 0.01, 0.02, 0.04, 0.08])
    ax.set_yticklabels(["0.006", "0.01", "0.02", "0.04", "0.08"])
    ax.minorticks_off()
    ax.set_title(r"(b)  shock thickness scales as $\nu$", loc="left", fontsize=9)
    ax.legend(fontsize=7, loc="upper left")

    fig.tight_layout(w_pad=1.4)
    style.finish(fig, os.path.join(FIGDIR, "fig03_shock_scaling"))
    print("\nwrote figures/fig02_shock_formation and fig03_shock_scaling")

    save_json("e2_shock_formation.json", {
        "nu": nu, "n": n, "t_s": 1.0,
        "strain_max_pred": float(strain_max_pred),
        "strain_max_measured": float(strain_t.max()),
        "width_vs_nu": [{"nu": r[0], "delta": r[1], "predicted": r[2], "N": r[3]}
                        for r in width_rows],
        "width_slope": float(slope),
    })
    save_csv("e2_strain_history.csv", ["t", "max_strain"],
             [(times[i], strain_t[i]) for i in range(len(times))])
    save_csv("e2_width_vs_nu.csv", ["nu", "delta", "predicted_4nu_over_U", "N"], width_rows)

    tee.close()
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
