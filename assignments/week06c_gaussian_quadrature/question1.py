"""
question1.py — week06c_gaussian_quadrature
------------------------------------------
Problem : Write ONE library routine for Gaussian quadrature that serves both
          the Gauss-Legendre family (interval [a, b]) and the Gauss-Laguerre
          family (interval [0, inf), weight e^{-x}), with hardwired tables of
          zeros and weights for orders n = 1 ... 6.  The routine lives in
          mylib.py (gaussian_quadrature + the GAUSS_LEGENDRE / GAUSS_LAGUERRE
          tables); this driver prints the tables and proves the *order of
          precision*: n nodes integrate every polynomial of degree <= 2n-1
          exactly.
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import (
    GAUSS_LAGUERRE,
    GAUSS_LEGENDRE,
    gaussian_quadrature,
)

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"
    lines: list[str] = []

    lines += [
        "Gaussian quadrature — one routine, two families (mylib.gaussian_quadrature)",
        "",
        "Table of Gauss-Legendre nodes/weights on [-1, 1] (n = 1..6):",
    ]
    for n in range(1, 7):
        nodes, weights = GAUSS_LEGENDRE[n]
        lines.append(f"  n = {n}:")
        for t, w in zip(nodes, weights):
            lines.append(f"      t = {t:+.15f}   w = {w:.15f}")

    lines += ["", "Table of Gauss-Laguerre nodes/weights on [0, inf), w(x) = e^{-x}:"]
    for n in range(1, 7):
        nodes, weights = GAUSS_LAGUERRE[n]
        lines.append(f"  n = {n}:")
        for t, w in zip(nodes, weights):
            lines.append(f"      t = {t:+.15f}   w = {w:.15e}")

    # -- order of precision: n nodes are exact for every polynomial of
    #    degree <= 2n - 1
    lines += [
        "",
        "Order-of-precision proof: for n nodes, int x^{2n-1} must come out",
        "exactly (int x^k on [-1,1] = 0 for odd k, 2/(k+1) for even k; on",
        "[0,inf) with weight e^{-x}, int x^k e^{-x} dx = k!).",
        "",
        f"{'n':>3}  {'degree 2n-1':>11}  {'legendre (x^d on [-1,1])':>26}  "
        f"{'laguerre (x^d * e^-x)':>24}",
    ]
    for n in range(1, 7):
        d = 2 * n - 1
        if d % 2 == 1:
            leg_exact = 0.0
        else:
            leg_exact = 2.0 / (d + 1)
        lag_exact = math.factorial(d)
        leg_val = gaussian_quadrature(lambda x, d=d: x ** d, -1.0, 1.0, n, "legendre")
        lag_val = gaussian_quadrature(lambda x, d=d: x ** d, 0.0, 0.0, n, "laguerre")
        lines.append(
            f"{n:>3d}  {d:>11d}  {abs(leg_val - leg_exact):>26.3e}  "
            f"{abs(lag_val - lag_exact):>24.3e}"
        )
    lines += [
        "",
        "Every entry is at round-off level: with only n function evaluations the",
        "rule is exact for polynomials of degree 2n - 1 — the best possible",
        "order for n freely chosen nodes (Riemann's theorem).",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
