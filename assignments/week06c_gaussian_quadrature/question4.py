"""
question4.py — week06c_gaussian_quadrature
------------------------------------------
Problem : Integrate  int_0^inf e^{-x} x^4 dx  with an appropriate n-point
          Gauss-Laguerre rule and compare with the analytical value 4! = 24.
          Since x^4 has degree 4, the 3-point rule (exact up to degree 5)
          must return the answer exactly.
Usage   : python question4.py [output_file]
          Defaults: output/q4_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import gaussian_quadrature

BASE_DIR = Path(__file__).resolve().parent
ANALYTICAL = 24.0  # 4!


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q4_output.txt"

    f = lambda x: x ** 4  # noqa: E731
    value = gaussian_quadrature(f, 0.0, 0.0, N=3, method="laguerre")

    lines = [
        "Integral of e^{-x} x^4 from 0 to infinity  (3-point Gauss-Laguerre)",
        "",
        "x^4 is a polynomial of degree 4, and an n-node Gauss rule is exact",
        "for degrees <= 2n - 1.  With n = 3 that ceiling is 5 >= 4, so the",
        "rule must reproduce the exact answer 4! = 24 with zero error.",
        "",
        f"Gauss-Laguerre (N = 3) = {value:.15f}",
        f"analytical value (4!)  = {ANALYTICAL:.1f}",
        f"difference             = {abs(value - ANALYTICAL):.3e}",
        f"relative error         = {abs(value - ANALYTICAL) / ANALYTICAL:.3e}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
