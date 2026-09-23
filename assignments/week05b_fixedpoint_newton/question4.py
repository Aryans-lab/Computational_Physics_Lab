"""
question4.py — week05b_fixedpoint_newton
----------------------------------------
Problem : Solve the coupled non-linear system
              F1(x, y) = 3x^2 - 2xy - 3 = 0
              F2(x, y) = 3y^2 - 4xy = 0
          with multivariate Newton-Raphson, using a numerical Jacobian
          (central finite differences), starting from (2.5, 3.1).
Usage   : python question4.py [output_file]
          Defaults: output/q4_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import jacobian, newton_raphson_system

BASE_DIR = Path(__file__).resolve().parent


def F(x: float, y: float) -> list[float]:
    return [3.0 * x * x - 2.0 * x * y - 3.0, 3.0 * y * y - 4.0 * x * y]


def J(x: float, y: float) -> list[list[float]]:
    """Numerical Jacobian via central differences (mylib.jacobian)."""
    return jacobian(F, [x, y])


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q4_output.txt"

    x, it = newton_raphson_system(F, J, [2.5, 3.1], accuracy=1e-10)
    residual = max(abs(v) for v in F(*x))

    lines = [
        "Multivariate Newton-Raphson (numerical Jacobian, central differences)",
        "",
        "System:",
        "  F1 = 3x^2 - 2xy - 3 = 0",
        "  F2 = 3y^2 - 4xy     = 0",
        "initial guess: (2.5, 3.1),  accuracy = 1e-10",
        "",
        f"Final root approximation: ({x[0]:.8f}, {x[1]:.8f}) after {it} iterations.",
        "exact solution:           (3.00000000, 4.00000000)",
        f"residual max |F(x*)| = {residual:.3e}",
        "",
        "Jacobian at the solution (numerical vs analytical):",
        f"  numerical:  [[{J(*x)[0][0]: .6f}, {J(*x)[0][1]: .6f}], [{J(*x)[1][0]: .6f}, {J(*x)[1][1]: .6f}]]",
        "  analytical: [[6x - 2y, -2x], [-4y, 6y - 4x]] at (3, 4) = [[12, -6], [-16, 12]]",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
