"""
E1 -- Are the solvers right?  Convergence, cost and stability.

This is the experiment that has to succeed before any physics can be claimed.
A smooth initial condition (nu = 0.1, integrated to T = 0.5, well before the
shock forms at t = 1) is compared with the exact Cole-Hopf solution.

Three questions:
  (i)   does each scheme converge at its advertised order in dx?
  (ii)  what does the same accuracy cost in time steps?
  (iii) where exactly does the explicit scheme become unstable, and how fast?

Run:  python experiments/e1_validation.py
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import style
from common import FIGDIR, RESDIR, save_json, save_csv, Tee
from exact import cole_hopf
from solvers import spectral, cn_imex, upwind, lax_wendroff, ftcs

style.apply()
import matplotlib.pyplot as plt

NU = 0.1
T = 0.5
RESOLUTIONS = [16, 32, 64, 128, 256]

DT_FIXED = 2e-4          # isolates the *spatial* error of every scheme

SCHEMES = [
    ("upwind ($N^{-1}$)", upwind, {"dt_fixed": DT_FIXED}, style.COL["red"], "upwind", 1.0),
    ("Lax-Wendroff ($N^{-2}$)", lax_wendroff, {"dt_fixed": DT_FIXED}, style.COL["orange"], "Lax-Wendroff", 2.4),
    ("Crank-Nicolson ($N^{-2}$)", cn_imex, {"dt_fixed": DT_FIXED}, style.COL["green"], "Crank-Nicolson", 0.42),
    ("spectral (exponential)", spectral, {"dt_fixed": DT_FIXED}, style.COL["blue"], "spectral", 1.0),
]


def l2(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


def main():
    tee = Tee(os.path.join(RESDIR, "e1_validation.txt"))
    sys.stdout = tee
    print("=" * 78)
    print("E1  Validation of the Burgers solvers   nu = 0.1,  T = 0.5,  u0 = sin x")
    print("=" * 78)
    print("The reference is the exact Cole-Hopf solution, computed on a 8192-point")
    print("quadrature grid with spectral interpolation of u0 (grid-independent to 1e-12).\n")

    convergence = {}
    for name, fn, kwargs, _, _short, _off in SCHEMES:
        errs, steps = [], []
        for n in RESOLUTIONS:
            x = np.arange(n) * 2.0 * np.pi / n
            u0 = np.sin(x)
            exact = cole_hopf(u0, x, T, NU, n_quad=8192)
            out = fn(u0, NU, T, n=n, **kwargs)
            errs.append(l2(out["u_final"], exact))
            steps.append(out["steps"])
            print(f"{name:26s} N={n:4d}  steps={out['steps']:6d}  L2 error={errs[-1]:.3e}")
        convergence[name] = {"N": RESOLUTIONS, "l2": errs, "steps": steps}

        err = np.array(errs)
        n = np.array(RESOLUTIONS, dtype=float)
        good = err > 1e-8
        order = -np.polyfit(np.log(n[good]), np.log(err[good]), 1)[0]
        print(f"{'':26s} fitted order: {order:.2f}\n")

    # ---------------------------------------------------------------- cost
    print("Cost of the same problem at N = 256, T = 0.5 (nu = 0.1):")
    x256 = np.arange(256) * 2.0 * np.pi / 256
    u0 = np.sin(x256)
    cost = {}
    for name, fn, kw in [("spectral, RK4 explicit", spectral, {"time_method": "rk4"}),
                         ("spectral, Strang split", spectral, {"time_method": "strang"}),
                         ("Crank-Nicolson", cn_imex, {}),
                         ("upwind", upwind, {})]:
        t0 = time.time()
        out = fn(u0, NU, T, n=256, **kw)
        cost[name] = {"steps": out["steps"], "seconds": time.time() - t0}
        print(f"  {name:24s} steps={out['steps']:6d}   wall={time.time()-t0:5.2f} s")
    print("  An explicit treatment of diffusion needs dt < dx^2/(2 nu); Crank-Nicolson")
    print("  has no such restriction, which is why it needs ~40x fewer steps here.\n")

    # ------------------------------------------------------------ stability
    # von Neumann analysis predicts that a centred explicit scheme amplifies
    # the grid-scale mode by |G| = |1 - 4 alpha| per step, alpha = nu dt/dx^2.
    # We seed that mode and measure the growth directly. Because sin(x)
    # dominates max|u|, the measurement uses the Nyquist amplitude itself.
    n = 256
    x = np.arange(n) * 2.0 * np.pi / n
    u0 = np.sin(x)
    sign = (-1.0) ** np.arange(n)
    dx = 2.0 * np.pi / n
    dt_crit = 0.5 * dx * dx / NU
    print(f"\nFTCS diffusion stability limit at N={n}: dt_crit = dx^2/(2 nu) = {dt_crit:.3e}")

    series, growth = {}, {}
    runs = [("FTCS, $0.9\\,\\Delta t_{\\rm crit}$", "ftcs", 0.9 * dt_crit, 3.0, style.COL["blue"]),
            ("FTCS, $1.1\\,\\Delta t_{\\rm crit}$", "ftcs", 1.1 * dt_crit, 3.0, style.COL["red"]),
            ("CN, $20\\,\\Delta t_{\\rm crit}$", "cn", 20.0 * dt_crit, 3.0, style.COL["green"])]

    for label, kind, dt, t_end, colour in runs:
        u = u0 + 1e-8 * sign
        hist_t, hist_a = [0.0], [abs(np.mean(u * sign))]
        t, k = 0.0, 0
        while t < t_end and k < 50000:
            if kind == "ftcs":
                up, um = np.roll(u, -1), np.roll(u, 1)
                u = u + dt * (-u * (up - um) / (2 * dx) + NU * (up - 2 * u + um) / dx ** 2)
            else:
                alpha = NU * dt / (2 * dx ** 2)
                nadv = -0.5 * u * (np.roll(u, -1) - np.roll(u, 1)) / dx
                up, um = np.roll(u, -1), np.roll(u, 1)
                rhs = u + alpha * (up - 2 * u + um) - dt * nadv
                a = np.full(n, -alpha); b = np.full(n, 1.0 + 2.0 * alpha)
                from solvers import _cyclic_tridiagonal
                u = _cyclic_tridiagonal(a, b, a, rhs, -alpha)
            t += dt
            k += 1
            hist_t.append(t)
            hist_a.append(abs(np.mean(u * sign)))
            if not np.all(np.isfinite(u)) or hist_a[-1] > 1e6:
                break
        series[label] = (np.array(hist_t), np.array(hist_a), colour)

        loga = np.log(np.clip(np.array(hist_a), 1e-300, None))
        good = (np.array(hist_a) > 1e-16) & (np.array(hist_a) < 1e-3)
        if good.sum() > 10:
            rate = np.polyfit(np.array(hist_t)[good], loga[good], 1)[0]
            alpha = NU * dt / (dx * dx)
            if kind == "ftcs":
                g = 1.0 - 4.0 * alpha                       # von Neumann, theta = pi
            else:
                g = (1.0 - 2.0 * alpha) / (1.0 + 2.0 * alpha)
            theory = np.log(abs(g)) / dt
            growth[label] = {"measured_per_time": float(rate), "von_neumann": float(theory), "dt": float(dt)}
            print(f"  {label:28s} dt={dt:.2e}  measured {rate:+.1f}/time   "
                  f"von Neumann {theory:+.1f}/time  (ratio {rate/theory:.3f})")
        else:
            growth[label] = {"measured_per_time": float("nan"), "dt": float(dt)}
            print(f"  {label:28s} dt={dt:.2e}  mode decays (unconditionally stable)")

    # ------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.4))

    ax = axes[0]
    for name, _, _, colour, short, off in SCHEMES:
        n = np.array(convergence[name]["N"], dtype=float)
        e = np.array(convergence[name]["l2"])
        ax.loglog(n, e, "o-", ms=3.0, color=colour, label=name, lw=1.0)
        # inline label, so the curves are never hidden by a legend box
        ax.text(300, e[-1] * off, "  " + short, color=colour, fontsize=6.5,
                va="center", ha="left")
    nn = np.array([16, 256], dtype=float)
    e0 = convergence["upwind ($N^{-1}$)"]["l2"][0]
    ax.loglog(nn, e0 * (nn / 16) ** -1, "--", color=style.COL["light"], lw=0.8)
    ax.loglog(nn, e0 * (nn / 16) ** -2, ":", color=style.COL["grey"], lw=0.8)
    ax.set_xlabel("grid points $N$")
    ax.set_ylabel(r"$L_2$ error vs exact")
    ax.set_ylim(1e-11, 1e-1)
    ax.set_xlim(13, 1000)
    style.log_ticks(ax, "x", [16, 64, 256])
    ax.set_title("(a)  spatial convergence", loc="left", fontsize=9)


    ax = axes[1]
    for label, (tt, aa, colour) in series.items():
        yy = np.where(aa > 2e-14, aa, np.nan)
        ax.semilogy(tt, yy, color=colour, label=label, lw=1.0)
    ax.set_xlabel("time $t$")
    ax.set_ylabel(r"Nyquist amplitude $|a_{N/2}|$")
    ax.set_ylim(2e-14, 1e0)
    ax.set_title("(b)  grid-scale mode", loc="left", fontsize=9)
    ax.legend(loc="lower left", fontsize=6.6, handlelength=1.2, frameon=True,
              framealpha=0.92, edgecolor="none", borderpad=0.4)

    ax = axes[2]
    names = list(cost)
    vals = [cost[k]["steps"] for k in names]
    cols = [style.COL["blue"], style.COL["purple"], style.COL["green"], style.COL["red"]]
    ax.barh(range(len(names)), vals, color=cols, height=0.62)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlim(30, max(vals) * 4)
    ax.set_xlabel("steps to $T=0.5$")
    for i, v in enumerate(vals):
        ax.text(v * 1.1, i, str(v), va="center", fontsize=7.5)
    ax.grid(axis="y", visible=False)
    ax.set_title("(c)  cost", loc="left", fontsize=9)

    fig.tight_layout(w_pad=1.6)
    style.finish(fig, os.path.join(FIGDIR, "fig01_validation"))
    print("\nwrote figures/fig01_validation.{pdf,png}")

    save_json("e1_validation.json", {
        "nu": NU, "T": T,
        "convergence": convergence,
        "cost": cost,
        "stability": {k.replace("$", "").replace("\\", ""): v for k, v in growth.items()},
        "dt_crit": float(dt_crit),
    })
    rows = []
    for name, d in convergence.items():
        for i in range(len(d["N"])):
            rows.append((name.replace("$", "").replace("\\", ""), d["N"][i], d["l2"][i], d["steps"][i]))
    save_csv("e1_convergence.csv", ["scheme", "N", "l2_error", "steps"], rows)

    tee.close()
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
