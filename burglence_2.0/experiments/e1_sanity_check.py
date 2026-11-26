# Experiment 1: sanity checks for diffusion, convergence order, and explicit stability

import json
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import burgers as b                    
import plotstyle as ps                         
from mylib import pde_heat_1d_crank_nicolson   

ps.setup()
results = {}

# (a) Pure diffusion against exact solution u = exp(-nu*t) * sin(x)
nu, T, N = 0.01, 1.0, 64
x_heat, snaps = pde_heat_1d_crank_nicolson(
    2.0 * math.pi, nu * T, N, 400, lambda x: math.sin(x), 0.0, 0.0, save_steps=[400]
)
u_num = np.array(snaps[400])
u_exact = math.exp(-nu * T) * np.sin(x_heat)
err_diff = float(np.max(np.abs(u_num - u_exact)))
results["diffusion_max_error"] = err_diff

# (b) Grid convergence before shock formation (t = 0.5) using N = 4096 as reference
T_conv, nu_conv = 0.5, 0.01
x_ref, dx_ref = b.make_grid(4096)
dt_ref = b.stable_dt("cn", b.initial_sine(x_ref), dx_ref, nu_conv) / 4.0
u_ref = b.solve("cn", N=4096, nu=nu_conv, T=T_conv, save_times=[T_conv], dt=dt_ref)["u"][-1]

Ns = [64, 128, 256, 512]
errors = [
    b.l2_difference(b.solve("cn", N=Nc, nu=nu_conv, T=T_conv, save_times=[T_conv])["u"][-1],
                    u_ref[::4096 // Nc])
    for Nc in Ns
]
orders = [math.log(errors[i] / errors[i + 1]) / math.log(2.0) for i in range(len(errors) - 1)]
results["convergence"] = {"N": Ns, "error": errors, "order": orders}

# (c) Explicit solver stability right below and above dt = dx^2 / (2*nu)
N_stab, nu_stab, n_steps = 512, 0.01, 400
x_s, dx_s = b.make_grid(N_stab)
dt_limit = dx_s ** 2 / (2.0 * nu_stab)

histories = {}
for factor in (0.9, 1.1):
    u = b.initial_sine(x_s)
    dt = factor * dt_limit
    track = [float(np.max(np.abs(u)))]
    for _ in range(n_steps):
        u = b.step_explicit(u, dx_s, dt, nu_stab)
        track.append(float(np.max(np.abs(u))) if np.all(np.isfinite(u)) else np.nan)
    histories[factor] = np.array(track)
    print(f"  dt = {factor} dt_limit : final max|u| = {track[-1]:.3e}")

results["stability"] = {
    "dt_limit": dt_limit,
    "final_0.9": float(np.nanmax(histories[0.9][-1:])),
    "final_1.1": float(np.nanmax(histories[1.1][-1:])),
}

# Plot all three panels
fig, axes = plt.subplots(1, 3, figsize=(10.0, 2.9))

ax = axes[0]
ax.plot(x_heat, u_exact, color=ps.COL["blue"], lw=2.4, label="exact")
ax.plot(x_heat[::2], u_num[::2], "o", ms=3.2, color=ps.COL["red"], label="Crank–Nicolson (mylib)")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$u$")
ax.set_title(r"(a) diffusion only, $\nu = 0.01$, $t = 1$", fontsize=9)
ax.legend(loc="lower left")
ax.text(0.03, 0.93, f"max error = {err_diff:.1e}", transform=ax.transAxes, fontsize=7.5, va="top")

ax = axes[1]
ax.loglog(Ns, errors, "o-", color=ps.COL["blue"], label="measured")
ax.loglog(Ns, errors[0] * (np.array(Ns) / Ns[0]) ** -2.0, "--", color=ps.COL["grey"], lw=1.0, label=r"slope $-2$")
ax.set_xlabel(r"$N$")
ax.set_xticks(Ns)
ax.set_xticklabels([str(n) for n in Ns])
ax.minorticks_off()
ax.set_ylabel(r"rms error at $t = 0.5$")
ax.set_title("(b) convergence before breaking", fontsize=9)
ax.legend()
ax.text(0.04, 0.10, "orders: " + ", ".join(f"{o:.2f}" for o in orders), transform=ax.transAxes, fontsize=7.5)

ax = axes[2]
for factor, colour in ((0.9, ps.COL["green"]), (1.1, ps.COL["red"])):
    ax.semilogy(histories[factor], color=colour, lw=1.2, label=rf"$dt = {factor}\,\Delta t_{{\rm limit}}$")
ax.set_xlabel("time step number")
ax.set_ylabel(r"$\max|u|$")
ax.set_title(r"(c) explicit scheme near $dt = dx^2/2\nu$", fontsize=9)
ax.set_ylim(1e-6, 1e6)
ax.legend(loc="center left")

fig.tight_layout()
ps.save(fig, "fig01_sanity")

os.makedirs(ps.RESDIR, exist_ok=True)
with open(os.path.join(ps.RESDIR, "e1_sanity.json"), "w") as f:
    json.dump(results, f, indent=2)

print("\ndiffusion error:", f"{err_diff:.2e}")
print("convergence orders:", ", ".join(f"{o:.2f}" for o in orders))
print("dt limit (N=512):", f"{dt_limit:.4e}")
print("max|u| at 0.9 limit:", f"{results['stability']['final_0.9']:.3f}")
print("max|u| at 1.1 limit:", f"{results['stability']['final_1.1']:.3e}")