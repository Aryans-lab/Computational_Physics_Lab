"""
question3.py — week06b_simpson_montecarlo
-----------------------------------------
Problem : Use Monte Carlo integration (with the in-house LCG pRNG) to
          estimate  int_{-1}^{1} sin^2(x) dx  accurate to 4 decimal places.
          Plot the integral value vs N and the absolute deviation vs N, and
          report the first N for which the estimate is within 1e-4 of the
          exact value.
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt, figures/q3_integral_vs_N.png,
                    figures/q3_deviation_vs_N.png
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe: never open a display
import matplotlib.pyplot as plt  # noqa: E402

from mylib import monte_carlo  # noqa: E402

BASE_DIR = Path(__file__).resolve().parent
ACCURACY = 1e-4
N_VALUES = [100 * i for i in range(11, 501)]  # 1100 ... 50000


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"
    fig1 = BASE_DIR / "figures" / "q3_integral_vs_N.png"
    fig2 = BASE_DIR / "figures" / "q3_deviation_vs_N.png"

    # Exact value: int sin^2 x dx = [x/2 - sin(2x)/4]_{-1}^{1} = 1 - sin(2)/2
    exact = 1.0 - math.sin(2.0) / 2.0

    integrals: list[float] = []
    deviations: list[float] = []
    seeds: list[int] = []

    for i, N in enumerate(N_VALUES):
        seed = i + 1
        value, _, _ = monte_carlo(lambda x: math.sin(x) ** 2, -1.0, 1.0, N, seed)
        integrals.append(value)
        deviations.append(abs(value - exact))
        seeds.append(seed)

    first_idx = next((i for i, d in enumerate(deviations) if d < ACCURACY), None)

    # -- figure 1: integral value vs N
    plt.figure(figsize=(8, 6))
    plt.plot(N_VALUES, integrals, marker="o", ms=4, label="Monte Carlo estimate")
    plt.axhline(exact, color="red", ls="--", label=f"exact = {exact:.8f}")
    plt.xlabel("N")
    plt.ylabel("integral value")
    plt.title("Monte Carlo integration: estimate vs N")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    fig1.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig1, dpi=150)
    plt.close()

    # -- figure 2: absolute deviation vs N
    plt.figure(figsize=(8, 6))
    plt.plot(N_VALUES, deviations, marker="o", ms=4, color="tab:blue",
             label="absolute deviation")
    plt.axhline(ACCURACY, color="red", ls="--", label=r"required accuracy $10^{-4}$")
    plt.xlabel("N")
    plt.ylabel("absolute deviation")
    plt.title("Monte Carlo integration: deviation vs N")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig2, dpi=150)
    plt.close()

    lines = [
        "Monte Carlo estimate of int_{-1}^1 sin^2(x) dx (in-house LCG)",
        f"exact value = 1 - sin(2)/2 = {exact:.10f}",
        "N runs: 1100 to 50000 in steps of 100, seeds 1..490",
        "",
    ]
    if first_idx is not None:
        lines += [
            f"first N with |estimate - exact| < {ACCURACY:g} :",
            f"  N = {N_VALUES[first_idx]}   (seed {seeds[first_idx]})",
            f"  estimate = {integrals[first_idx]:.8f}",
            f"  deviation  = {deviations[first_idx]:.3e}",
        ]
    else:
        lines.append(f"accuracy of {ACCURACY:g} was not reached up to N = 50000")
    lines += [
        "",
        "The deviation wanders inside the ~1/sqrt(N) envelope (central limit",
        "theorem): the error does not shrink smoothly, it is a statistical",
        "fluctuation — hence the wiggly lines, in contrast to the smooth",
        "error decay of the deterministic quadrature rules in week06a.",
        "",
        f"Figures: {fig1.name}, {fig2.name} (in figures/)",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")
    print(f"Figures written to: {fig1}, {fig2}")


if __name__ == "__main__":
    main()
