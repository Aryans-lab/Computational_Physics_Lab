"""
question2.py — week06a_midpoint_trapezoid
-----------------------------------------
Problem : With both the Midpoint and the Trapezoidal rule, evaluate
              int_1^2 (1/x) dx,   int_0^{pi/2} x cos x dx,   int_0^1 x arctan(x) dx
          for N = 4, 8, 15, 20 and compare with the analytical results in
          tabular form.
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import midpoint, trapezoidal

BASE_DIR = Path(__file__).resolve().parent

INTEGRALS = [
    ("1/x  on [1, 2]", lambda x: 1.0 / x, 1.0, 2.0, math.log(2.0)),
    ("x cos x  on [0, pi/2]", lambda x: x * math.cos(x), 0.0, math.pi / 2.0,
     math.pi / 2.0 - 1.0),
    ("x arctan x  on [0, 1]", lambda x: x * math.atan(x), 0.0, 1.0,
     math.pi / 4.0 - 0.5),
]
NS = [4, 8, 15, 20]


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"
    lines: list[str] = []

    for label, f, a, b, exact in INTEGRALS:
        lines += [
            f"Integral: {label}",
            f"exact value = {exact:.10f}",
            "",
            f"{'N':>5}  {'midpoint':>16}  {'|error|':>11}  {'trapezoidal':>16}  {'|error|':>11}",
        ]
        for N in NS:
            m = midpoint(f, a, b, N)
            t = trapezoidal(f, a, b, N)
            lines.append(
                f"{N:>5d}  {m:>16.10f}  {abs(m - exact):>11.3e}  "
                f"{t:>16.10f}  {abs(t - exact):>11.3e}"
            )
        lines.append("")

    lines += [
        "Analytical values used:  ln 2 = 0.6931471806,",
        "  int x cos x = [x sin x + cos x]_0^{pi/2} = pi/2 - 1 = 0.5707963268,",
        "  int x arctan x = [(x^2/2) arctan x - x/2 + arctan x/2]_0^1 = pi/4 - 1/2 = 0.2853981634.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
