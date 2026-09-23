"""
question1.py — week05c_polynomial_roots
---------------------------------------
Problem : Write the library functions for Laguerre's method and synthetic
          deflation, finding real roots only, taking only the coefficients
          as input.  The routines live in mylib.py (laguerre,
          synthetic_division, laguerre_roots, polynomial_value); this driver
          validates them on a cubic with known roots and checks every root
          by substitution.
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import laguerre_roots, polynomial_value

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"

    # P(x) = x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
    coeffs = [1.0, -6.0, 11.0, -6.0]
    roots, iters = laguerre_roots(coeffs, b0=2.5)
    exact = [1.0, 2.0, 3.0]

    lines = [
        "Validation of Laguerre's method + synthetic deflation",
        "P(x) = x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3),  initial guess b0 = 2.5",
        "",
        f"{'root':>18}  {'iterations':>10}  |P(root)|    error vs exact",
    ]
    for r, it, e in sorted(zip(roots, iters, exact), key=lambda t: t[0]):
        lines.append(
            f"{r:>18.12f}  {it:>10d}  {abs(polynomial_value(coeffs, r)):>11.3e}  "
            f"{abs(r - e):.3e}"
        )
    lines += [
        "",
        "Deflation flow: each root is divided out with synthetic division",
        "(Horner's scheme), and the next root search warm-starts from the root",
        "just found — the last (linear) factor is solved exactly.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
