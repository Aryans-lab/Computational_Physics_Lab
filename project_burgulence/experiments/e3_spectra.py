"""
E3 -- Fourier analysis of a shock train: the k^-2 spectrum and its cut-off.

A periodic shock train is, in the language of Fourier analysis, a sawtooth.
Its spectrum is a power law,

        E(k) ~ k^-2 ,

over the whole range of scales between the domain size and the shock
thickness, and then falls off exponentially once k exceeds 1/delta. The
viscosity therefore does not change the shape of the inertial range -- it only
decides where the range ends. That is the statement this experiment measures.

Run:  python experiments/e3_spectra.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import style
from common import FIGDIR, RESDIR, save_json, save_csv, Tee
from diagnostics import energy_spectrum, shock_width, fit_line
from solvers import spectral
from exact import tanh_width

style.apply()
import matplotlib.pyplot as plt

NV = [0.01, 0.005, 0.002]
N = 4096
T = 2.0


def main():
    tee = Tee(os.path.join(RESDIR, "e3_spectra.txt"))
    sys.stdout = tee
    print("=" * 78)
    print("E3  Energy spectra of a shock train")
    print("=" * 78)

    spectra, rows = {}, []
    for nu in NV:
        x = np.arange(N) * 2.0 * np.pi / N
        out = spectral(np.sin(x), nu, T, n=N, cfl=0.3)
        u = out["u_final"]
        E = energy_spectrum(u)
        k = np.arange(E.size)
        delta = shock_width(u, 2.0 * np.pi / N, n_shocks=1)[0][1]
        k_cut = 1.0 / delta

        # inertial range: between the domain scale and a third of the cut-off
        k1, k2 = 4, max(int(3 * k_cut / 10), 8)
        slope, intercept, err = fit_line(np.log(k[k1:k2]), np.log(E[k1:k2]))
        # exponential tail: E k^2 ~ exp(-k delta)
        tail = (k > 0.5 * k_cut) & (k < 0.9 * k * 1.0 + 1e-9) & (E > 0)
        tail = (k > 0.6 * k_cut) & (k < 0.95 * k_cut) & (E * k ** 2 > 0)
        m_tail, c_tail, e_tail = fit_line(k[tail], np.log(E[tail] * k[tail] ** 2))

        spectra[nu] = (k, E, delta, k_cut)
        rows.append((nu, delta, k_cut, slope, err, -m_tail))
        print(f"nu={nu:6.4f}  delta={delta:.4f}  k_cut=1/delta={k_cut:7.1f}  "
              f"inertial slope={slope:.3f}+-{err:.3f}  tail decay length={-m_tail:.4f}  "
              f"(4nu/U={tanh_width(nu):.4f})")

    nus = np.array(NV)
    kcuts = np.array([spectra[nu][3] for nu in NV])
    m, c, e = fit_line(np.log(nus), np.log(kcuts))
    print(f"\n  k_cut ~ nu^{m:.3f}  (prediction: k_cut = U/(4 nu), slope -1)")
    print(f"  slope of the inertial range: {np.mean([r[3] for r in rows]):.3f} (prediction -2)")

    # ------------------------------------------------- resolution artefact
    print("\nResolution study at nu = 0.002 (physical shock width 0.008):")
    res_rows = []
    for n in [256, 512, 1024, 2048, 4096]:
        x = np.arange(n) * 2.0 * np.pi / n
        out = spectral(np.sin(x), 0.002, T, n=n, cfl=0.3)
        E = energy_spectrum(out["u_final"])
        k = np.arange(E.size)
        # fit the same decade of k in every run: 4 <= k <= 40
        m1 = (k >= 4) & (k <= 40)
        s1, _, _ = fit_line(np.log(k[m1]), np.log(E[m1]))
        d = shock_width(out["u_final"], 2.0 * np.pi / n, n_shocks=1)[0][1]
        res_rows.append((n, s1, d))
        print(f"  N={n:5d}  points per shock = {d / (2*np.pi/n):5.2f}   "
              f"apparent slope (4<k<40) = {s1:.3f}   measured delta = {d:.4f}")

    # ------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45))

    ax = axes[0]
    for nu in NV:
        k, E, delta, k_cut = spectra[nu]
        good = E > 1e-14
        ax.loglog(k[good], E[good], color=style.CYCLE[NV.index(nu)], lw=1.0,
                  label=rf"$\nu={nu}$")
    kk = np.array([5, 80], dtype=float)
    ax.loglog(kk, 2e-2 * kk ** -2, "--", color=style.COL["grey"], lw=0.9)
    ax.text(55, 2e-2 * 55 ** -2 * 2.2, r"$k^{-2}$", color=style.COL["grey"], fontsize=8)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$E(k)$")
    ax.set_xlim(1, 1e3)
    ax.set_ylim(1e-13, 1)
    ax.set_title("(a)  energy spectra", loc="left", fontsize=9)
    ax.set_ylim(1e-14, 3)
    ax.legend(fontsize=6.8, loc="lower left", handlelength=1.0)

    ax = axes[1]
    for nu in NV:
        k, E, delta, k_cut = spectra[nu]
        good = (k > 0) & (E * k ** 2 > 1e-12)
        ax.semilogy(k[good], (E * k ** 2)[good], color=style.CYCLE[NV.index(nu)], lw=1.0)
        ax.axvline(k_cut, color=style.CYCLE[NV.index(nu)], ls=":", lw=0.7)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$E(k)\,k^{2}$")
    ax.set_xlim(1, 1e3)
    ax.set_ylim(1e-8, 1e0)
    ax.set_title(r"(b)  compensated spectrum", loc="left", fontsize=9)
    style.log_ticks(ax, "x", [10, 100, 1000])

    ax = axes[2]
    ax.loglog(nus, kcuts, "o", color=style.COL["blue"], label="measured $1/\\delta$")
    ax.loglog(nus, 1.0 / (4 * nus), "--", color=style.COL["grey"], lw=0.9,
              label=r"$U/4\nu$")
    ax.set_xlabel(r"$\nu$")
    ax.set_ylabel(r"$k_{\rm cut}$")
    style.log_ticks(ax, "x", [0.002, 0.005, 0.01])
    ax.set_xticklabels(["0.002", "0.005", "0.01"])
    ax.set_title("(c)  cut-off scales as $1/\\nu$", loc="left", fontsize=9)
    ax.legend(fontsize=7, loc="lower left")

    fig.tight_layout(w_pad=1.4)
    style.finish(fig, os.path.join(FIGDIR, "fig04_spectra"))
    print("\nwrote figures/fig04_spectra.{pdf,png}")

    save_json("e3_spectra.json", {
        "nu_list": NV, "N": N, "T": T,
        "rows": [{"nu": r[0], "delta": r[1], "k_cut": r[2], "slope": r[3],
                  "slope_err": r[4], "tail_length": r[5]} for r in rows],
        "kcut_vs_nu_slope": float(m),
        "resolution_study": [{"N": r[0], "apparent_slope": r[1], "delta": r[2]} for r in res_rows],
    })
    save_csv("e3_spectra_summary.csv",
             ["nu", "delta", "k_cut", "inertial_slope", "slope_err", "tail_length"], rows)
    save_csv("e3_resolution.csv", ["N", "apparent_slope_k4_40", "measured_delta"], res_rows)

    tee.close()
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
