"""
question2.py — week06b_simpson_montecarlo
-----------------------------------------
Problem : Evaluate
              int_1^2 (1/x) dx   and   int_0^{pi/2} x cos x dx
          with the Midpoint rule and with Simpson's 1/3 rule, accurate to
          6 decimal places.  The number of subintervals N is determined from
          the classical error-bound formulae (mylib.integration_error_bound_N)
          rather than by trial and error.
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import (
    integration_error_bound_N,
    midpoint,
    simpson,
)

BASE_DIR = Path(__file__).resolve().parent
TARGET = 1e-6

# (label, integrand, a, b, M2 = max|f''|, M4 = max|f''''|, exact value)
PROBLEMS = [
    ("1/x on [1, 2]", lambda x: 1.0 / x, 1.0, 2.0,
     2.0,          # max|f''|  of 1/x on [1,2] is 2 (at x = 1)
     24.0,         # max|f''''| of 1/x on [1,2] is 24 (at x = 1)
     math.log(2.0)),
    ("x cos x on [0, pi/2]", lambda x: x * math.cos(x), 0.0, math.pi / 2.0,
     2.3,          # max |f''| = max |2 sin x + x cos x| ~= 2.3 on [0, pi/2]
     4.2,          # max |f''''| = max |4 sin x + x cos x| ~= 4.2 on [0, pi/2]
     math.pi / 2.0 - 1.0),
]


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"
    lines: list[str] = [
        f"Midpoint + Simpson 1/3, target error {TARGET:g} (from the error bounds)",
        "",
    ]

    for label, f, a, b, M2, M4, exact in PROBLEMS:
        n_mid = integration_error_bound_N("midpoint", a, b, M2, TARGET)
        n_simp = integration_error_bound_N("simpson", a, b, M4, TARGET)

        mid = midpoint(f, a, b, n_mid)
        simp, evals = simpson(f, a, b, n_simp)

        lines += [
            f"Integral: {label}   (exact = {exact:.10f})",
            f"  derivative bounds: max|f''| = {M2:g},  max|f''''| = {M4:g}",
            f"  N from error bound (midpoint,    M2)  = {n_mid}",
            f"  N from error bound (Simpson 1/3, M4)  = {n_simp}",
            "",
            f"  Midpoint   (N = {n_mid:>4d}): {mid:.10f}   |error| = {abs(mid - exact):.3e}",
            f"  Simpson 1/3 (N = {n_simp:>4d}): {simp:.10f}   |error| = {abs(simp - exact):.3e}"
            f"   ({evals} function evaluations)",
            "",
        ]

    lines += [
        "As expected for the same error target, Simpson needs far fewer",
        "subintervals (O(eps^{-1/4}) vs O(eps^{-1/2}) in N) — the O(h^4)",
        "rule pays off as soon as the integrand is C^4.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
