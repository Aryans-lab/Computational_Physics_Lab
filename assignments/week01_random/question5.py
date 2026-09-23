"""
question5.py — week01_random
----------------------------
Problem : Generate an exponentially distributed pRNG with PDF exp(-x) from a
          uniform pRNG on [0, 1) using the inverse-transform method
          y = -ln(u), and histogram at least 5,000 samples.
Usage   : python question5.py [output_file]
          Defaults: output/q5_output.txt, figure figures/q5_exponential_distribution.png
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
N_SAMPLES = 5000
LAM = 1.0
SEED = 1


def exponential_sample(n: int, lam: float, seed: int) -> list[float]:
    """n exponential deviates via y = -(1/lam) ln(u), u in (0, 1]."""
    myrand_reset(seed)
    samples: list[float] = []
    for _ in range(n):
        u = myrand()
        while u <= 0.0:  # guard the measure-zero u = 0
            u = myrand()
        samples.append(-(1.0 / lam) * math.log(u))
    return samples


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q5_output.txt"
    fig_file = BASE_DIR / "figures" / "q5_exponential_distribution.png"

    samples = exponential_sample(N_SAMPLES, LAM, SEED)
    mean = sum(samples) / N_SAMPLES
    var = sum((s - mean) ** 2 for s in samples) / N_SAMPLES

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(samples, bins=50, color="skyblue", edgecolor="black", density=True)
    ax.plot(
        [x * 0.02 for x in range(1, 301)],
        [LAM * math.exp(-LAM * x) for x in [xx * 0.02 for xx in range(1, 301)]],
        color="red", lw=2, label="exact e^{-x}",
    )
    ax.set_xlim(0, 6)
    ax.set_xlabel("y")
    ax.set_ylabel("probability density")
    ax.set_title(f"Exponential distribution from the LCG (N = {N_SAMPLES})")
    ax.legend()
    fig.tight_layout()
    fig_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_file, dpi=120)
    plt.close(fig)

    lines = [
        f"Inverse-transform sampling: y = -(1/lambda) ln(u),  lambda = {LAM}",
        f"N = {N_SAMPLES} samples (seed = {SEED})",
        "",
        f"sample mean   = {mean:.5f}   (exact 1/lambda = {1.0 / LAM:.5f})",
        f"sample stddev = {math.sqrt(var):.5f}   (exact 1/lambda = {1.0 / LAM:.5f})",
        "For the exponential distribution the mean equals the standard deviation,",
        "a clean one-number check that the transformation worked.",
        "",
        f"Figure: {fig_file.relative_to(BASE_DIR)}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")
    print(f"Figure written to: {fig_file}")


if __name__ == "__main__":
    main()
