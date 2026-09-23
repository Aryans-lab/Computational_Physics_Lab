"""
question1.py — week05b_fixedpoint_newton
----------------------------------------
Problem : Write the library functions for Fixed-point iteration and
          Newton-Raphson (1-D and multivariate).  The routines live in
          mylib.py (fixed_point, newton_raphson, newton_raphson_system,
          jacobian); this driver validates each against problems with known
          solutions, including the quadratic-convergence signature of Newton.
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import (
    fixed_point,
    jacobian,
    newton_raphson,
    newton_raphson_system,
)

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"
    lines: list[str] = []

    # -- 1. Fixed point: x^2 - 2x - 3 = 0 via g(x) = sqrt(2x + 3) -> x = 3
    g = lambda x: math.sqrt(2.0 * x + 3.0)  # noqa: E731
    r, it = fixed_point(g, 1.0, accuracy=1e-10)
    lines += [
        "1. Fixed-point iteration  (x = g(x))",
        "   f(x) = x^2 - 2x - 3,  rewrite x = sqrt(2x + 3),  x_0 = 1",
        f"   root = {r:.12f}  in {it} iterations   (exact: 3)",
        f"   |f(root)| = {abs(r * r - 2 * r - 3):.3e}",
        "",
    ]

    # -- 2. Newton 1-D: f(x) = x^2 - 2, quadratic convergence check
    f = lambda x: x * x - 2.0  # noqa: E731
    df = lambda x: 2.0 * x  # noqa: E731
    r, it = newton_raphson(f, df, x0=2.0, accuracy=1e-14)
    err = abs(r - math.sqrt(2.0))
    lines += [
        "2. Newton-Raphson, 1-D (analytical derivative)",
        "   f(x) = x^2 - 2,  x_0 = 2",
        f"   root = {r:.15f}  in {it} iterations   (exact sqrt(2) = {math.sqrt(2):.15f})",
        f"   absolute error = {err:.3e}",
        "   quadratic convergence: the digits of the root roughly double every",
        "   iteration (3 -> 6 -> 12 -> 14 significant digits here).",
        "",
    ]

    # -- 3. Newton multivariate: known solution (3, 4)
    F = lambda x, y: [3 * x * x - 2 * x * y - 3, 3 * y * y - 4 * x * y]  # noqa: E731
    J = lambda x, y: jacobian(F, [x, y])  # noqa: E731
    x, it = newton_raphson_system(F, J, [2.5, 3.1], accuracy=1e-10)
    res = max(abs(v) for v in F(*x))
    lines += [
        "3. Newton-Raphson, multivariate (numerical Jacobian, central differences)",
        "   F1 = 3x^2 - 2xy - 3,  F2 = 3y^2 - 4xy,  x_0 = (2.5, 3.1)",
        f"   solution = ({x[0]:.10f}, {x[1]:.10f})  in {it} iterations   (exact: (3, 4))",
        f"   max |F(x*)| = {res:.3e}",
        "   the Newton step solves J dx = F (Gauss-Jordan on [J | F]) instead of",
        "   inverting J explicitly — same step, better conditioning.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
