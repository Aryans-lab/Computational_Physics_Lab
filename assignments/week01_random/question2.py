"""
question2.py — week01_random
----------------------------
Problem : Write a Linear Congruential Generator (LCG) with the GCC parameters
          a = 1103515245, c = 12345, m = 32768, and check its correlation by
          plotting x_i vs x_{i+k} for k = 5.
          The generator is mylib.myrand (full-period: the parameters satisfy
          the Hull-Dobell conditions for m = 2^15).
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt, figure figures/q2_correlation_lcg.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe: never open a display
import matplotlib.pyplot as plt  # noqa: E402

from mylib import LCG_A, LCG_C, LCG_M, myrand, myrand_reset  # noqa: E402

BASE_DIR = Path(__file__).resolve().parent
N_SAMPLES = 2000
K = 5
X0 = 1.2  # seed = int(x0) = 1 (the LCG state is an integer modulo m)


def lcg_numbers(n: int, seed: int) -> list[float]:
    """Draw n uniform deviates from the LCG, seeded reproducibly."""
    myrand_reset(seed)
    return [myrand() for _ in range(n)]


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"
    fig_file = BASE_DIR / "figures" / "q2_correlation_lcg.png"

    deviates = lcg_numbers(N_SAMPLES, seed=int(X0))

    x_vals = deviates[:-K]
    y_vals = deviates[K:]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(x_vals, y_vals, s=2)
    ax.set_xlabel("x_i")
    ax.set_ylabel("x_{i+k}")
    ax.set_title(f"LCG correlation, k = {K}, N = {N_SAMPLES}")
    fig.tight_layout()
    fig_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_file, dpi=120)
    plt.close(fig)

    # Simple statistical diagnostics of the sequence.
    mean = sum(deviates) / N_SAMPLES
    var = sum((d - mean) ** 2 for d in deviates) / N_SAMPLES

    lines = [
        f"LCG: x_{{i+1}} = ({LCG_A} * x_i + {LCG_C}) mod {LCG_M}  (GCC parameters)",
        f"Seed: {int(X0)} (x_0 = {X0} truncated), N = {N_SAMPLES}, lag k = {K}",
        "",
        "Full-period check (Hull-Dobell, m = 2^15):",
        f"  a == 1 (mod 4)      : {LCG_A % 4 == 1}",
        f"  c odd               : {LCG_C % 2 == 1}",
        f"  m a power of two    : {LCG_M & (LCG_M - 1) == 0}",
        "  => period = m = 32768 (verified in the test suite)",
        "",
        f"Sample mean = {mean:.6f}   (exact 0.5)",
        f"Sample variance = {var:.6f}   (exact 1/12 = {1/12:.6f})",
        "",
        "Correlation plot (x_i vs x_{i+5}) shows the known lattice structure",
        "of low-bit LCGs — acceptable for these physics simulations, but a",
        "modern PRNG (xorshift/PCG) would scatter the points uniformly.",
        "",
        f"Figure: {fig_file.relative_to(BASE_DIR)}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")
    print(f"Figure written to: {fig_file}")


if __name__ == "__main__":
    main()
