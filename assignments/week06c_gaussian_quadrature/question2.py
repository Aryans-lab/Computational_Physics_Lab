"""
question2.py — week06c_gaussian_quadrature
------------------------------------------
Problem : Use a 4-point Gauss-Legendre quadrature to evaluate
              int_{-1}^{1} x^2 / (1 + x^4) dx
          and compare with Simpson's 1/3 rule at comparable accuracy.  The
          true value is 0.4874954944... — a useful third number to see how
          both approximations sit around the exact answer.
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import gaussian_quadrature, simpson

BASE_DIR = Path(__file__).resolve().parent
TRUE_VALUE = 0.4874954943993610  # high-precision reference (independent quad)


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"

    f = lambda x: x * x / (1.0 + x ** 4)  # noqa: E731

    gq = gaussian_quadrature(f, -1.0, 1.0, N=4, method="legendre")
    simp_val, simp_evals = simpson(f, -1.0, 1.0, N=8)

    lines = [
        "Integral of x^2/(1 + x^4) from -1 to 1",
        f"true value (high-precision reference) = {TRUE_VALUE:.12f}",
        "",
        f"Gauss-Legendre (N = 4,  {4} function evaluations)      = {gq:.12f}",
        f"Simpson's 1/3   (N = 8,  {simp_evals} function evaluations)   = {simp_val:.12f}",
        "",
        f"|Gauss-Legendre - true| = {abs(gq - TRUE_VALUE):.4e}",
        f"|Simpson 1/3   - true| = {abs(simp_val - TRUE_VALUE):.4e}",
        "",
        "Notes:",
        "  * x^2/(1+x^4) is smooth, so both rules do well, but it is NOT a",
        "    polynomial: the 4-node rule is exact only up to degree 7, and",
        "    both errors are genuinely ~1e-3 here.",
        "  * Gauss-Legendre reaches this with 4 evaluations vs 9 for Simpson;",
        "    the gap widens rapidly for higher precision (order 2n-1 vs 4).",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
