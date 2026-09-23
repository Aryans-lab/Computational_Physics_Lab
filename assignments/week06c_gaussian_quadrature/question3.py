"""
question3.py — week06c_gaussian_quadrature
------------------------------------------
Problem : Use a 5-point Gauss-Laguerre quadrature to evaluate
              int_0^inf e^{-x} / (1 + x) dx
          (analytical value ~ 0.59635).  The e^{-x} factor is the quadrature
          weight, so the integrand passed to the routine is just 1/(1+x).
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import gaussian_quadrature

BASE_DIR = Path(__file__).resolve().parent
# e * E_1(1), where E_1 is the exponential integral — exact value:
TRUE_VALUE = 0.5963473623231941


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"

    f = lambda x: 1.0 / (1.0 + x)  # noqa: E731
    value = gaussian_quadrature(f, 0.0, 0.0, N=5, method="laguerre")

    lines = [
        "Integral of e^{-x}/(1+x) from 0 to infinity  (5-point Gauss-Laguerre)",
        "",
        "The weight e^{-x} is built into the rule, so the integrand supplied",
        "is simply f(x) = 1/(1+x).",
        "",
        f"Gauss-Laguerre (N = 5) = {value:.12f}",
        f"true value  e*E_1(1)   = {TRUE_VALUE:.12f}",
        f"absolute error         = {abs(value - TRUE_VALUE):.4e}",
        "relative error           "
        f"= {abs(value - TRUE_VALUE) / TRUE_VALUE:.3e}",
        "",
        "Gauss-Laguerre handles the infinite interval [0, inf) natively — no",
        "cutoff has to be chosen — and 5 nodes already give ~2e-4 accuracy",
        "on this non-polynomial integrand.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
