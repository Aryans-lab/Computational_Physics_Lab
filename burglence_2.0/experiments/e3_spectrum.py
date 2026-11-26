# Experiment 3: Fourier energy spectrum and dissipative cutoff of the shock

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import burgers as b        # noqa: E402
import plotstyle as ps     # noqa: E402

ps.setup()
results = {}

N, T = 2048, 2.0
KMIN, KMAX = 4, 12          # fitting window inside the power law range

spectra = {}
for nu in (0.01, 0.005):
    x, dx = b.make_grid(N)
    u = b.solve("cn", N=N, nu=nu, T=T, save_times=[T], dt=0.6 * dx)["u"][-1]
    k, E = b.energy_spectrum(u, dx)
    slope, intercept = b.fit_power_law(k, E, KMIN, KMAX)

    # Cutoff is the highest k where E * k^2 stays above 70% of its plateau
    plateau = float(np.median(E[(k >= 3) & (k <= 8)] * k[(k >= 3) & (k <= 8)] ** 2))
    k_cut = float(k[np.where(E * k ** 2 > 0.7 * plateau)[0][-1]])

    spectra[nu] = {"k": k, "E": E, "slope": slope, "intercept": intercept, "k_cut": k_cut, "plateau": plateau}
    print(f"  nu = {nu:5.3f}:  slope = {slope:+.3f}   k_cut = {k_cut:.1f}   (1/delta = {1.0/(4*nu/0.75):.1f})")

    results[f"nu_{nu}"] = {"slope": float(slope), "k_cut": k_cut, "N": N, "T": T, "kmin": KMIN, "kmax": KMAX}

# Plot E(k) and compensated E(k)*k^2
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.2))
colours = {0.01: ps.COL["red"], 0.005: ps.COL["blue"]}

ax = axes[0]
for nu, data in spectra.items():
    ax.loglog(data["k"][1:], data["E"][1:], ".", ms=3, color=colours[nu], label=rf"$\nu = {nu:g}$")
kg = np.array([2.5, 30.0])
ax.loglog(kg, spectra[0.005]["E"][4] * (kg / 4.0) ** -2.0, "--", color=ps.COL["grey"], lw=1.0, label=r"$k^{-2}$")
ax.set_xlim(2, 200)
ax.set_ylim(1e-7, 1e-1)
ax.set_xlabel(r"$k$")
ax.set_ylabel(r"$E(k) = |c_k|^2$")
ax.set_title(r"(a) both slopes are $-2$", fontsize=9)
ax.legend(loc="lower left")
ax.text(0.03, 0.30, "fitted slopes: " + ", ".join(f"{d['slope']:.2f}" for d in spectra.values()),
        transform=ax.transAxes, fontsize=7.5)

ax = axes[1]
for nu, data in spectra.items():
    k, E = data["k"], data["E"]
    ax.semilogx(k[1:], E[1:] * k[1:] ** 2, ".", ms=3, color=colours[nu], label=rf"$\nu = {nu:g}$")
    ax.axvline(data["k_cut"], color=colours[nu], ls=":", lw=0.9)
    ax.text(data["k_cut"] * 1.05, 3e-3, rf"$k_{{\rm cut}}={data['k_cut']:.0f}$", fontsize=7, color=colours[nu])
ax.set_xlim(2, 200)
ax.set_ylim(0, 0.12)
ax.set_xlabel(r"$k$")
ax.set_ylabel(r"$E(k)\,k^2$")
ax.set_title(r"(b) flat = $k^{-2}$; the knee is the shock", fontsize=9)
ax.legend(loc="lower left")

fig.tight_layout()
ps.save(fig, "fig03_spectrum")

os.makedirs(ps.RESDIR, exist_ok=True)
with open(os.path.join(ps.RESDIR, "e3_spectrum.json"), "w") as f:
    json.dump(results, f, indent=2)

print("\n  fitted window:", KMIN, "< k <", KMAX)
for nu, d in spectra.items():
    print(f" nu = {nu:5.3f}: slope {d['slope']:+.3f}, cut-off k ~ {d['k_cut']:.0f}")
