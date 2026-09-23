"""
question3.py — week05b_fixedpoint_newton
----------------------------------------
Problem : Find the root of  f(x) = x^2 - 2x - 3  using the Fixed-point
          method.  For the root x = 3, use  g(x) = sqrt(2x + 3);  for the
          root x = -1, use  g(x) = 3 / (x - 2)  (both satisfy |g'| < 1 near
          their respective fixed points).
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import fixed_point

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"

    g_pos = lambda x: math.sqrt(2.0 * x + 3.0)    # x^2 = 2x + 3 -> x = sqrt(2x+3)  (root 3)
    g_neg = lambda x: 3.0 / (x - 2.0)             # x(x-2) = 3  -> x = 3/(x-2)     (root -1)
    f = lambda x: x * x - 2.0 * x - 3.0  # noqa: E731

    r_pos, it_pos = fixed_point(g_pos, 1.0, accuracy=1e-10)
    r_neg, it_neg = fixed_point(g_neg, 0.0, accuracy=1e-10)

    lines = [
        "Fixed-point method for f(x) = x^2 - 2x - 3 = 0   (exact roots: -1 and 3)",
        "",
        "Root x = 3,  with g(x) = sqrt(2x + 3),  x_0 = 1:",
        f"   the root is {r_pos:.12f}",
        f"   the number of iterations taken is {it_pos}",
        f"   |f(root)| = {abs(f(r_pos)):.3e}",
        "   (converges because |g'(3)| = 1/3 < 1)",
        "",
        "Root x = -1,  with g(x) = 3/(x - 2),  x_0 = 0:",
        f"   the root is {r_neg:.12f}",
        f"   the number of iterations taken is {it_neg}",
        f"   |f(root)| = {abs(f(r_neg)):.3e}",
        "   (converges because |g'(-1)| = 1/9 < 1)",
        "",
        "The choice of the rewrite x = g(x) is the whole art of the method:",
        "g(x) = x^2 - 2x - 3 + x would diverge from both roots (|g'| > 1 there).",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
