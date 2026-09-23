"""
question2.py — week05a_bisection_regf
-------------------------------------
Problem : Find the root of  f(x) = log(x/2) - sin(5x/2)  to an accuracy of
          1e-6 in the interval [1.0, 1.5] (f changes sign there), using both
          Bisection and Regula Falsi.
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import bisection, regula_falsi

BASE_DIR = Path(__file__).resolve().parent
ACCURACY = 1e-6


def f(x: float) -> float:
    import math

    return math.log(x / 2.0) - math.sin(5.0 * x / 2.0)


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"

    a, b = 1.0, 1.5
    assert f(a) * f(b) < 0.0, "bracket must contain a sign change"

    r_b, it_b = bisection(f, a, b, accuracy=ACCURACY)
    r_r, it_r = regula_falsi(f, a, b, accuracy=ACCURACY)
    # high-accuracy reference for error reporting
    r_ref, _ = bisection(f, a, b, accuracy=1e-13, max_iterations=200)

    lines = [
        f"Root of f(x) = log(x/2) - sin(5x/2) on [{a}, {b}],  accuracy = {ACCURACY:g}",
        f"(f({a}) = {f(a):+.4f}, f({b}) = {f(b):+.4f}  — sign change present)",
        "",
        "1. Bisection method:",
        f"   the root is approximately {r_b:.12f}",
        f"   the number of iterations taken is {it_b}",
        f"   |f(root)| = {abs(f(r_b)):.3e}",
        "",
        "2. Regula Falsi method:",
        f"   the root is approximately {r_r:.12f}",
        f"   the number of iterations taken is {it_r}",
        f"   |f(root)| = {abs(f(r_r)):.3e}",
        "",
        f"reference (bisection to 1e-13): {r_ref:.12f}",
        f"  bisection error  = {abs(r_b - r_ref):.3e}",
        f"  regula-falsi error = {abs(r_r - r_ref):.3e}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
