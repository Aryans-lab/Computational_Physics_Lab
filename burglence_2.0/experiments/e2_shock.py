# Experiment 2: shock formation, slope saturation, tanh profile, and thickness scaling

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

ps.setup()
results = {}

nu, N, T, t_s = 0.01, 2048, 3.0, 1.0
x, dx = b.make_grid(N)
times = [0.5, 1.0, 2.0, 3.0]
hist_t = np.linspace(0.0, T, 60)

# Run once and store snapshots at requested times
run = b.solve("cn", N=N, nu=nu, T=T, save_times=list(times) + list(hist_t), dt=0.6 * dx)
saved = dict(zip(run["t"], run["u"]))
u_final = run["u"][-1]
results["run"] = {"N": N, "nu": nu, "T": T, "dt": run["dt"], "nsteps": run["nsteps"]}

# Compare final shock near x = pi with exact tanh formula
U = 0.5 * (np.max(u_final) - np.min(u_final))
delta_theory = 4.0 * nu / U
delta_measured = b.shock_width(u_final, dx)
u_tanh = -U * np.tanh(U * (x - math.pi) / (2.0 * nu))
err_tanh = float(np.max(np.abs(u_final - u_tanh)[np.abs(x - math.pi) < 0.1]))

results["shock"] = {
    "U": float(U),
    "delta_measured": float(delta_measured),
    "delta_theory": float(delta_theory),
    "ratio": float(delta_measured / delta_theory),
    "tanh_max_error": err_tanh,
    "max_slope_final": float(b.max_slope(u_final, dx)),
    "plateau_theory": float(U ** 2 / (2 * nu)),
}

# Sweep viscosity nu to check delta = 4*nu/U (finer grid for smaller nu)
nu_list = [0.04, 0.02, 0.01, 0.005]
N_list = [512, 512, 1024, 2048]
widths, U_list, res_list = [], [], []
for nui, Ni in zip(nu_list, N_list):
    xi, dxi = b.make_grid(Ni)
    ui = b.solve("cn", N=Ni, nu=nui, T=2.0, save_times=[2.0], dt=0.6 * dxi)["u"][-1]
    widths.append(float(b.shock_width(ui, dxi)))
    Ui = float(0.5 * (np.max(ui) - np.min(ui)))
    U_list.append(Ui)
    res_list.append(4.0 * nui / Ui / dxi)
    print(f"  nu = {nui:5.3f}  N = {Ni:5d}  delta = {widths[-1]:.4f}"
          f"  4nu/U = {4*nui/Ui:.4f}  ({res_list[-1]:.1f} points)")

slope_w = math.log(widths[-1] / widths[0]) / math.log(nu_list[-1] / nu_list[0])
results["width_sweep"] = {
    "nu": nu_list, "N": N_list, "delta": widths,
    "U": U_list, "points_per_shock": res_list, "slope": slope_w,
}

# Plot 2x2 figure
fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.6))

ax = axes[0][0]
for t, colour in zip([0.0, 0.5, 1.0, 2.0, 3.0],
                     [ps.COL["grey"], ps.COL["green"], ps.COL["orange"], ps.COL["purple"], ps.COL["blue"]]):
    u = saved[min(saved, key=lambda s: abs(s - t))]
    ax.plot(x, u, color=colour, lw=1.2, label="initial" if t == 0 else f"$t = {t:g}$")
ax.axvline(math.pi, color=ps.COL["red"], lw=0.8, ls=":")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$u(x,t)$")
ax.set_title(r"(a) steepening and breaking, $\nu = 10^{-2}$", fontsize=9)
ax.legend(loc="lower left", ncol=1, fontsize=7)

ax = axes[0][1]
t_all = np.array(run["t"][1:])
slope_all = np.array([b.max_slope(u, dx) for u in run["u"][1:]])
U_all = np.array([0.5 * (np.max(u) - np.min(u)) for u in run["u"][1:]])
ax.semilogy(t_all, slope_all, "o", ms=2.5, color=ps.COL["blue"], label="measured")
tt = np.linspace(0.02, 0.95, 80)
ax.semilogy(tt, 1.0 / (t_s - tt), "-", color=ps.COL["grey"], lw=1.0, label=r"$1/(t_s-t)$ (inviscid)")
ax.semilogy(t_all, U_all ** 2 / (2 * nu), "-", color=ps.COL["red"], lw=1.1, label=r"$U(t)^2/2\nu$ (viscous)")
ax.axvline(t_s, color=ps.COL["red"], lw=0.8, ls=":")
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$\max|u_x|$")
ax.set_ylim(0.8, 200)
ax.set_title("(b) strain grows, then saturates", fontsize=9)
ax.legend(loc="lower right")

ax = axes[1][0]
sel = np.abs(x - math.pi) < 0.09
ax.plot((x[sel] - math.pi) * 1e3, u_final[sel], "o", ms=3.0, color=ps.COL["red"], label="numerical, $t = 3$")
ax.plot((x[sel] - math.pi) * 1e3, u_tanh[sel], "--", color=ps.COL["blue"], lw=1.2, label=r"$-U\tanh[U(x-\pi)/2\nu]$")
ax.set_xlabel(r"$(x - \pi)\times 10^{3}$")
ax.set_ylabel(r"$u$")
ax.set_title("(c) the profile IS a tanh", fontsize=9)
ax.legend(loc="upper right")

ax = axes[1][1]
nu_arr = np.array(nu_list)
ax.loglog(nu_arr, widths, "o", ms=5, color=ps.COL["blue"], label="measured")
ax.loglog(nu_arr, 4 * nu_arr / np.array(U_list), "s--", ms=4, lw=1.0, color=ps.COL["red"], label=r"$4\nu/U$")
ax.set_xlabel(r"$\nu$")
ax.set_xticks(nu_list)
ax.set_xticklabels([f"{v:g}" for v in nu_list])
ax.minorticks_off()
ax.set_ylabel(r"shock thickness  $\delta$")
ax.set_title(rf"(d) $\delta \propto \nu^{{{slope_w:.2f}}}$", fontsize=9)
ax.legend(loc="upper left")

fig.tight_layout()
ps.save(fig, "fig02_shock")

os.makedirs(ps.RESDIR, exist_ok=True)
with open(os.path.join(ps.RESDIR, "e2_shock.json"), "w") as f:
    json.dump(results, f, indent=2)

print("\nU:", f"{U:.3f}")
print("delta measured:", f"{delta_measured:.4f}")
print("delta = 4 nu / U:", f"{delta_theory:.4f}", f"(ratio {delta_measured/delta_theory:.3f})")
print("tanh max error:", f"{err_tanh:.3e}")
print("max|u_x| at t=3:", f"{b.max_slope(u_final, dx):.1f}", "  U^2/2nu =", f"{U**2/(2*nu):.1f}")
print("delta ~ nu^slope:", f"{slope_w:.3f}")
