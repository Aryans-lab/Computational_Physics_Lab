"""
question4.py — week01_random
----------------------------
Problem : Simulate the radioactive decay chain A -> B -> C with decay
          constants lambda_A = 0.0347, lambda_B = 0.0198, starting from
          (N_A, N_B, N_C) = (500, 0, 0).  At each time step every nucleus
          independently "decays" with probability lambda*dt (decided by the
          in-house LCG), so the ensemble follows the Bateman equations in
          ensemble average without any ODE solver.  Plot N_A, N_B, N_C vs t.
Usage   : python question4.py [output_file]
          Defaults: output/q4_output.txt, figure figures/q4_radioactive_decay.png
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe: never open a display
import matplotlib.pyplot as plt  # noqa: E402

from mylib import myrand, myrand_reset  # noqa: E402

BASE_DIR = Path(__file__).resolve().parent
N_A0, N_B0, N_C0 = 500, 0, 0
LAMBDA_A, LAMBDA_B = 0.0347, 0.0198
T_MAX, DT = 400.0, 0.1
SEED = 1


def radioactive_decay_simulation(
    na0: int, nb0: int, nc0: int,
    lambda_a: float, lambda_b: float, t_max: float, dt: float, seed: int,
) -> tuple[list[float], list[int], list[int], list[int]]:
    """Kinetic Monte-Carlo of A -> B -> C with per-step decay probabilities."""
    na, nb, nc = na0, nb0, nc0
    time, n_a, n_b, n_c = [0.0], [na0], [nb0], [nc0]
    n_steps = int(round(t_max / dt))

    # Pre-generate the random pool once: at most na0 + nb0 <= na0 fresh draws
    # per step are needed, since A + B never exceeds the initial inventory.
    myrand_reset(seed)
    pool = [myrand() for _ in range(na0 * n_steps)]
    idx = 0

    for _ in range(n_steps):
        draws_a = [pool[idx + j] for j in range(na)]
        idx += na
        draws_b = [pool[idx + j] for j in range(nb)]
        idx += nb

        decays_a = sum(1 for r in draws_a if r < lambda_a * dt)
        decays_b = sum(1 for r in draws_b if r < lambda_b * dt)

        na -= decays_a
        nb += decays_a - decays_b
        nc += decays_b
        time.append(time[-1] + dt)
        n_a.append(na)
        n_b.append(nb)
        n_c.append(nc)
    return time, n_a, n_b, n_c


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q4_output.txt"
    fig_file = BASE_DIR / "figures" / "q4_radioactive_decay.png"

    time, n_a, n_b, n_c = radioactive_decay_simulation(
        N_A0, N_B0, N_C0, LAMBDA_A, LAMBDA_B, T_MAX, DT, SEED
    )

    # Exact Bateman curves for overlay (deterministic reference).
    exact_a = [N_A0 * math.exp(-LAMBDA_A * t) for t in time]
    exact_b = [
        N_A0 * LAMBDA_A / (LAMBDA_B - LAMBDA_A)
        * (math.exp(-LAMBDA_A * t) - math.exp(-LAMBDA_B * t))
        for t in time
    ]

    i_pk = max(range(len(n_b)), key=lambda i: n_b[i])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(time, n_a, label="N_A (MC)", marker="o", ms=2)
    ax.plot(time, n_b, label="N_B (MC)", marker="o", ms=2)
    ax.plot(time, n_c, label="N_C (MC)", marker="o", ms=2)
    ax.plot(time, exact_a, "--", color="gray", alpha=0.6, label="N_A exact")
    ax.plot(time, exact_b, "--", color="gray", alpha=0.6, label="N_B exact")
    ax.set_xlabel("time")
    ax.set_ylabel("number of nuclei")
    ax.set_title("Radioactive decay chain A -> B -> C (kinetic Monte Carlo)")
    ax.legend(loc="center right", fontsize=8)
    fig.tight_layout()
    fig_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_file, dpi=120)
    plt.close(fig)

    lines = [
        "Kinetic Monte-Carlo decay chain A -> B -> C",
        f"lambda_A = {LAMBDA_A},  lambda_B = {LAMBDA_B},  dt = {DT},  t_max = {T_MAX}",
        f"initial conditions: N_A = {N_A0}, N_B = {N_B0}, N_C = {N_C0}  (seed = {SEED})",
        "",
        f"peak of N_B: {n_b[i_pk]} nuclei at t = {time[i_pk]:.1f}",
        f"  (Bateman maximum t = ln(la/lnb)/(lnb-lna) "
        f"= {math.log(LAMBDA_A / LAMBDA_B) / (LAMBDA_B - LAMBDA_A):.1f})",
        f"final counts: N_A = {n_a[-1]}, N_B = {n_b[-1]}, N_C = {n_c[-1]}",
        f"conservation check: N_A + N_B + N_C = {n_a[-1] + n_b[-1] + n_c[-1]} (must be {N_A0})",
        "",
        "The MC trace (points) follows the exact Bateman solution (dashed)",
        "with statistical fluctuations ~ sqrt(N) — no ODE was solved.",
        "",
        f"Figure: {fig_file.relative_to(BASE_DIR)}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")
    print(f"Figure written to: {fig_file}")


if __name__ == "__main__":
    main()
