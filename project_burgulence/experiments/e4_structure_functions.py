"""
E4 -- Burgulence: structure functions and the bifractal exponents.

With a random, band-limited initial condition the flow breaks into a train of
shocks of random strength and separation. There is no longer a single shock to
measure, so we ask a statistical question instead: how does the p-th moment of
the velocity increment grow with separation r,

        S_p(r) = < |u(x+r) - u(x)|^p >_x  ~  r^{zeta_p} ?

For Burgers turbulence the answer is known exactly in the limit nu -> 0:

        zeta_p = p     for p <= 1,        zeta_p = 1     for p >= 1.

The reason is that a shock contributes a jump of fixed size, so every moment
p >= 1 is dominated by the same set of discontinuities. A Kolmogorov-type
argument would instead give zeta_p = p/3; the two predictions differ by a
factor of three at p = 4, which is what the fit has to resolve.

Run:  python experiments/e4_structure_functions.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import style
from common import FIGDIR, RESDIR, save_json, save_csv, Tee, random_ic
from diagnostics import structure_functions, local_slope, energy_spectrum
from solvers import spectral

style.apply()
import matplotlib.pyplot as plt

NU = 0.002
N = 2048
T = 2.0
M = 24                      # realisations
ORDERS = (0.5, 1, 2, 3, 4)
R_FIT = (30, 120)           # inertial window, in grid points
R_FIT_ALT = [(20, 80), (30, 120), (40, 160)]   # for the systematic


def main():
    tee = Tee(os.path.join(RESDIR, "e4_structure_functions.txt"))
    sys.stdout = tee
    print("=" * 78)
    print("E4  Burgulence: structure functions")
    print(f"    nu={NU}, N={N}, T={T}, {M} realisations")
    print("=" * 78)

    fields = []
    for m in range(M):
        x, u0 = random_ic(N, seed=1000 + 7 * m, n_modes=8, rms=1.0)
        out = spectral(u0, NU, T, n=N, cfl=0.3)
        fields.append(out["u_final"])
        if (m + 1) % 6 == 0:
            print(f"  {m+1}/{M} realisations done")

    r, S = structure_functions(fields, 2.0 * np.pi / N, orders=ORDERS, r_min=2, r_max=N // 6)

    print("\n  p    zeta_p (fit)    stat err   window spread   min(p,1)   p/3")
    zetas, zeta_err, zeta_sys = {}, {}, {}
    dx = 2.0 * np.pi / N
    for p in ORDERS:
        z, e = local_slope(r, S[p], np.searchsorted(r, R_FIT[0] * dx),
                           np.searchsorted(r, R_FIT[1] * dx))
        alts = [local_slope(r, S[p], np.searchsorted(r, a * dx),
                            np.searchsorted(r, b * dx))[0]
                for a, b in R_FIT_ALT]
        zetas[p], zeta_err[p] = z, e
        zeta_sys[p] = float(np.std(alts))
        print(f"  {p:3.1f}   {z:6.3f} +- {e:.3f}      {zeta_sys[p]:.3f}        "
              f"{min(p,1.0):6.3f}     {p/3:6.3f}")

    # ------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.55))

    ax = axes[0]
    for i, p in enumerate(ORDERS):
        ax.loglog(r, S[p], "o", ms=2.0, color=style.CYCLE[i], label=rf"$p={p}$")
        seg = (r >= dx * R_FIT[0]) & (r <= dx * R_FIT[1])
        ax.loglog(r[seg], S[p][seg], "-", color=style.CYCLE[i], lw=0.7)
    ax.axvspan(dx * R_FIT[0], dx * R_FIT[1],
               color=style.COL["light"], alpha=0.18, lw=0)
    ax.set_xlabel(r"separation $r$")
    ax.set_ylabel(r"$S_p(r)$")
    ax.set_title("(a)  structure functions", loc="left", fontsize=9)
    ax.legend(fontsize=6.8, loc="lower right", ncol=2, handlelength=1.0)

    ax = axes[1]
    ps = np.array(sorted(zetas))
    zs = np.array([zetas[p] for p in ps])
    es = np.array([np.hypot(zeta_err[p], zeta_sys[p]) for p in ps])
    ax.errorbar(ps, zs, yerr=es, fmt="o", ms=4, color=style.COL["blue"],
                capsize=2, label="measured")
    pp = np.linspace(0, 4.2, 200)
    ax.plot(pp, np.minimum(pp, 1.0), "-", color=style.COL["red"], lw=1.0,
            label=r"Burgers: $\min(p,1)$")
    ax.plot(pp, pp / 3.0, "--", color=style.COL["grey"], lw=0.9,
            label=r"K41: $p/3$")
    ax.set_xlabel(r"order $p$")
    ax.set_ylabel(r"$\zeta_p$")
    ax.set_xlim(0, 4.3)
    ax.set_ylim(0, 2.2)
    ax.set_title("(b)  measured scaling exponents", loc="left", fontsize=9)
    ax.legend(fontsize=7, loc="upper left")

    fig.tight_layout(w_pad=1.4)
    style.finish(fig, os.path.join(FIGDIR, "fig05_structure_functions"))
    print("\nwrote figures/fig05_structure_functions.{pdf,png}")

    save_json("e4_structure_functions.json", {
        "nu": NU, "N": N, "T": T, "M": M,
        "zetas": {str(p): zetas[p] for p in zetas},
        "zeta_err": {str(p): zeta_err[p] for p in zeta_err},
        "zeta_systematic": {str(p): zeta_sys[p] for p in zeta_sys},
        "fit_window_grid_points": R_FIT,
    })
    save_csv("e4_zeta.csv", ["p", "zeta", "stat_err", "window_spread", "min_p_1", "p_over_3"],
             [(p, zetas[p], zeta_err[p], zeta_sys[p], min(p, 1.0), p / 3.0) for p in ps])

    tee.close()
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
