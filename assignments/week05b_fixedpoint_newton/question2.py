"""
question2.py — week05b_fixedpoint_newton
----------------------------------------
Problem : Find the root of  f(x) = 3x + sin(x) - e^x  in [-1.5, 1.5] to an
          accuracy of 1e-6 using Bisection, Regula Falsi and Newton-Raphson
          (with analytical *and* numerical derivative), and compare the
          convergence rates.
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import bisection, newton_raphson, regula_falsi

BASE_DIR = Path(__file__).resolve().parent
ACCURACY = 1e-6


def f(x: float) -> float:
    return 3.0 * x + math.sin(x) - math.exp(x)


def df(x: float) -> float:
    return 3.0 + math.cos(x) - math.exp(x)


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"

    a, b = -1.5, 1.5
    assert f(a) * f(b) < 0.0, "bracket must contain a sign change"

    r_b, it_b = bisection(f, a, b, accuracy=ACCURACY)
    r_r, it_r = regula_falsi(f, a, b, accuracy=ACCURACY)
    r_n, it_n = newton_raphson(f, df, x0=0.0, accuracy=ACCURACY)
    r_nn, it_nn = newton_raphson(f, df=None, x0=0.0, accuracy=ACCURACY)
    r_ref, _ = bisection(f, a, b, accuracy=1e-13, max_iterations=200)

    lines = [
        f"Root of f(x) = 3x + sin(x) - e^x on [{a}, {b}],  accuracy = {ACCURACY:g}",
        "",
        "1. Bisection method:",
        f"   the root is approximately {r_b:.12f}",
        f"   the number of iterations taken is {it_b}",
        "",
        "2. Regula Falsi method:",
        f"   the root is approximately {r_r:.12f}",
        f"   the number of iterations taken is {it_r}",
        "",
        "3. Newton-Raphson method (analytical derivative f' = 3 + cos x - e^x):",
        f"   the root is {r_n:.12f}",
        f"   the number of iterations taken is {it_n}",
        "",
        "4. Newton-Raphson (numerical derivative, central difference, h = 1e-5):",
        f"   the root is {r_nn:.12f}",
        f"   the number of iterations taken is {it_nn}",
        "",
        f"reference (bisection to 1e-13): {r_ref:.12f}",
        "",
        "Convergence-rate comparison (iterations to 1e-6):",
        f"   bisection      {it_b:>4d}   (linear,  ~log2(3/1e-6) = {math.log2(3e6):.1f} needed)",
        f"   regula falsi   {it_r:>4d}   (superlinear, but one endpoint sticks)",
        f"   Newton         {it_n:>4d}   (quadratic — the error squares each step)",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
