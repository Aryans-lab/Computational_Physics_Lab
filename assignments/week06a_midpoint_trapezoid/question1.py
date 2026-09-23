"""
question1.py — week06a_midpoint_trapezoid
-----------------------------------------
Problem : Write library functions for the Midpoint (rectangle) and
          Trapezoidal integration rules.  The routines live in mylib.py
          (midpoint, trapezoidal); this driver validates both on
          int_1^2 dx/x = ln 2 and demonstrates the O(h^2) convergence:
          doubling N should shrink the error by ~4.
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import midpoint, trapezoidal

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"

    f = lambda x: 1.0 / x  # noqa: E731
    exact = math.log(2.0)
    Ns = [1, 2, 4, 8, 16]

    lines = [
        "Validation of midpoint & trapezoidal rules on int_1^2 (1/x) dx = ln 2",
        f"exact value = {exact:.12f}",
        "",
        f"{'N':>5}  {'midpoint':>16}  {'error':>11}  {'trapezoid':>16}  {'error':>11}  {'mid ratio':>10}  {'trap ratio':>11}",
    ]
    prev = None
    for N in Ns:
        m = midpoint(f, 1.0, 2.0, N)
        t = trapezoidal(f, 1.0, 2.0, N)
        em, et = abs(m - exact), abs(t - exact)
        rm = f"{em / prev[0]:.3f}" if prev else "   -"
        rt = f"{et / prev[1]:.3f}" if prev else "   -"
        lines.append(f"{N:>5d}  {m:>16.10f}  {em:>11.3e}  {t:>16.10f}  {et:>11.3e}  {rm:>10}  {rt:>11}")
        prev = (em, et)
    lines += [
        "",
        "error-ratio column: ~0.25 (= 1/4) on doubling N — the O(h^2) = O(1/N^2)",
        "order of both rules, exactly as the error bounds promise:",
        "  midpoint:    |E| <= (b-a)^3 / (24 N^2) * max|f''|",
        "  trapezoidal: |E| <= (b-a)^3 / (12 N^2) * max|f''|",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
