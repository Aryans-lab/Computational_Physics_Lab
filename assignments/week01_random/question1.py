"""
question1.py — week01_random
----------------------------
Problem : Use the logistic map  x_{i+1} = c x_i (1 - x_i)  to generate 1,000
          "random" numbers from x_0 = 0.1, and expose the correlation between
          successive values by plotting x_i vs x_{i+k}.  Three chaos
          parameters c are used, with lags k = 3, 5, 10.
          (The same routine is the workhorse of the correlation study in
          question2.py, where a full-period LCG is compared.)
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt, figure figures/q1_correlation_chaos.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe: never open a display
import matplotlib.pyplot as plt  # noqa: E402

BASE_DIR = Path(__file__).resolve().parent

C_VALUES = [3.85, 3.95, 4.0]   # chaotic regime (c = 4 is fully chaotic)
K_VALUES = [3, 5, 10]          # lags
N_SAMPLES = 1000
X0 = 0.1


def random_number_iterative(x0: float, c: float, n: int) -> list[float]:
    """Generate n iterates of the logistic map x_{i+1} = c x_i (1 - x_i)."""
    random_numbers: list[float] = []
    xi = x0
    for _ in range(n):
        xi = c * xi * (1.0 - xi)
        random_numbers.append(xi)
    return random_numbers


def plot_panel(numbers: list[float], k: int, c: float, ax: plt.Axes) -> None:
    """Scatter x_i vs x_{i+k} on the given axes."""
    x_vals = numbers[:-k]
    y_vals = numbers[k:]
    ax.scatter(x_vals, y_vals, s=2)
    ax.set_xlabel("x_i")
    ax.set_ylabel("x_{i+k}")
    ax.set_title(f"c={c}, k={k}")


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"
    fig_file = BASE_DIR / "figures" / "q1_correlation_chaos.png"

    fig, axes = plt.subplots(3, 3, figsize=(10, 9))
    fig.suptitle(f"Logistic-map correlation study (N={N_SAMPLES})")

    for i, c in enumerate(C_VALUES):
        numbers = random_number_iterative(X0, c, N_SAMPLES)
        for j, k in enumerate(K_VALUES):
            plot_panel(numbers, k, c, axes[i][j])

    fig.tight_layout()
    fig_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_file, dpi=120)
    plt.close(fig)

    lines = [
        "Logistic map x_{i+1} = c x_i (1 - x_i), x_0 = 0.1, N = 1000",
        "Correlation plots x_i vs x_{i+k}, k = 3, 5, 10 (see figures/).",
        "",
        "Parameters studied: c = " + ", ".join(str(c) for c in C_VALUES),
        "  c = 3.85, 3.95: chaotic but with visible lag structure at small k.",
        "  c = 4.0       : fully chaotic; successive iterates decorrelate fast.",
        "",
        f"Figure: {fig_file.relative_to(BASE_DIR)}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")
    print(f"Figure written to: {fig_file}")


if __name__ == "__main__":
    main()
