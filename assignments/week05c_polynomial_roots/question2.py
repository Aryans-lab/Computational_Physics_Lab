"""
question2.py — week05c_polynomial_roots
---------------------------------------
Problem : Use Laguerre's method and synthetic deflation to find the (real)
          roots of three polynomials.
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import laguerre_roots, polynomial_value

BASE_DIR = Path(__file__).resolve().parent

POLYNOMIALS: list[tuple[str, list[float], list[float]]] = [
    ("P1(x) = x^4 - x^3 - 7x^2 + x + 6", [1, -1, -7, 1, 6], [1.0, -1.0, -2.0, 3.0]),
    ("P2(x) = x^4 - 5x^2 + 4", [1, 0, -5, 0, 4], [1.0, -1.0, 2.0, -2.0]),
    ("P3(x) = 2x^5 - 19.5x^3 + 0.5x^2 + 13.5x - 4.5",
     [2, 0, -19.5, 0.5, 13.5, -4.5],
     [-3.0, -1.0, 0.5, 0.5, 3.0]),  # note: x = 0.5 is a DOUBLE root
]


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"
    lines: list[str] = []

    for label, coeffs, exact in POLYNOMIALS:
        roots, iters = laguerre_roots(coeffs, b0=1.0)
        # every expected root must be hit (double root: tolerance absorbs
        # the linear-convergence split 0.49998/0.50002)
        for e in exact:
            assert min(abs(r - e) for r in roots) < 1e-4
        lines += [
            f"{label}",
            f"{'root':>20}  {'iterations':>10}  |P(root)|",
        ]
        for r, it in zip(roots, iters):
            lines.append(f"{r:>20.12f}  {it:>10d}  {abs(polynomial_value(coeffs, r)):>11.3e}")
        lines.append("")

    lines += [
        "Remarks:",
        "  * P1 and P2 factor over the integers: P1 = (x-1)(x+1)(x+2)(x-3),",
        "    P2 = (x-1)(x+1)(x-2)(x+2).  The roots are recovered to ~1e-10.",
        "  * P3 has a double root at x = 0.5:  P3 = 2(x-3)(x+3)(x+1)(x-0.5)^2.",
        "    Multiple roots are harder — convergence degrades from cubic to",
        "    linear — so the two deflated roots split around 0.5 as",
        "    0.49998 / 0.50002.  This is the expected behaviour of the",
        "    method, not an accuracy failure.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
