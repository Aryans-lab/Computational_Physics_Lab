"""
question3.py — week01_random
----------------------------
Problem : Estimate pi by the "throwing darts" method: throw N uniform points
          in the unit square and count the fraction that fall inside the
          quarter circle x^2 + y^2 <= 1, so that  pi ~ 4 * (hits / N).
          The estimates for N = 20 ... 5000 are plotted (with residuals) to
          visualise the O(N^{-1/2}) statistical error scaling.
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt, figure figures/q3_pi_estimate.png
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
N_MIN, N_MAX = 20, 5000
SEED = 1


def estimate_pi(n: int, seed: int) -> float:
    """Monte-Carlo pi estimate with n throws.

    Each N reseeds the LCG with the same seed, so the estimate for N uses the
    first N samples of one fixed sequence — a *nested* design, which is what
    makes the convergence plot meaningful.
    """
    myrand_reset(seed)
    hits = 0
    for _ in range(n):
        x = myrand()
        y = myrand()
        hits += x * x + y * y <= 1.0
    return 4.0 * hits / n


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"
    fig_file = BASE_DIR / "figures" / "q3_pi_estimate.png"

    n_values = list(range(N_MIN, N_MAX + 1))
    estimates = [estimate_pi(n, SEED) for n in n_values]
    residuals = [p - math.pi for p in estimates]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.scatter(n_values, estimates, s=1)
    ax1.axhline(math.pi, color="red", ls="--", label="true pi")
    ax1.set_ylabel("estimated pi")
    ax1.set_title("Monte-Carlo pi estimation (quarter-circle throws)")
    ax1.legend()
    ax2.scatter(n_values, residuals, s=1, color="tab:blue")
    ax2.set_xlabel("number of throws N")
    ax2.set_ylabel("residual (estimate - true pi)")
    ax2.set_title("Residual plot — error envelope ~ 1/sqrt(N)")
    fig.tight_layout()
    fig_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_file, dpi=120)
    plt.close(fig)

    report_N = [100, 250, 500, 1000, 2500, 5000]
    lines = [
        f"Monte-Carlo pi: quarter circle of unit radius, seed = {SEED}",
        "pi_est(N) = 4 * (points inside x^2 + y^2 <= 1) / N",
        "",
        f"{'N':>6}  {'pi estimate':>14}  {'|error|':>12}  {'~1/sqrt(N)':>12}",
    ]
    for n in report_N:
        p = estimates[n - N_MIN]
        lines.append(f"{n:>6d}  {p:>14.8f}  {abs(p - math.pi):>12.2e}  {1/math.sqrt(n):>12.2e}")
    lines += [
        "",
        f"Final estimate at N = {N_MAX}: {estimates[-1]:.6f}   (true pi = {math.pi:.6f})",
        "The scatter of the residuals traces the Central-Limit-Theorem bound",
        "sigma(pi_est) = 2/sqrt(N): halving the error requires 4x the throws.",
        "",
        f"Figure: {fig_file.relative_to(BASE_DIR)}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")
    print(f"Figure written to: {fig_file}")


if __name__ == "__main__":
    main()
