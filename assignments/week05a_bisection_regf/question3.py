"""
question3.py — week05a_bisection_regf
-------------------------------------
Problem : For  f(x) = -x - cos(x)  find an appropriate interval bracket in
          which a root lies, *starting from the interval [2, 3]* (which
          contains no sign change), and then refine with Regula Falsi.
          This exercises the automatic bracket-expansion step of the
          library (mylib.bracket_root).
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import bracket_root, regula_falsi

BASE_DIR = Path(__file__).resolve().parent


def f(x: float) -> float:
    return -x - math.cos(x)


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"

    a, b = 2.0, 3.0
    assert f(a) * f(b) >= 0.0, "the starting bracket must have NO sign change"
    a, b = bracket_root(f, a, b)
    r, it = regula_falsi(f, a, b, accuracy=1e-8)

    # reference: f is monotone (f' = -1 + sin x <= 0), so the bracket contains
    # exactly one root; find it to high precision for the error report.
    r_ref = bisection_reference(f, a, b)

    lines = [
        f"Bracketing f(x) = -x - cos(x), starting interval [{2.0}, {3.0}]",
        f"(f(2) = {f(2.0):+.4f}, f(3) = {f(3.0):+.4f} — no sign change)",
        "",
        f"after automatic expansion: the root is in the interval [{a:.6f}, {b:.6f}]",
        "(expansion pushes out the side with the smaller |f|, growing the step",
        " by 10% each push, until f(a) and f(b) have opposite signs)",
        "",
        "Regula Falsi on the bracketed interval (accuracy = 1e-8):",
        f"   the root is approximately {r:.12f}",
        f"   the number of iterations taken is {it}",
        f"   |f(root)| = {abs(f(r)):.3e}",
        f"   error vs high-accuracy reference = {abs(r - r_ref):.3e}",
        "",
        "f is monotone (f' = -1 + sin x <= 0), so the bracket contains exactly",
        "one root — the negative one, near -0.739 (the other root of",
        "-x = cos(x) is at +0.739 and lies outside the expanded interval).",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


def bisection_reference(f, a: float, b: float) -> float:
    """Plain bisection to 1e-14 for a reference root (test-local helper)."""
    for _ in range(100):
        c = 0.5 * (a + b)
        if f(c) * f(a) < 0.0:
            b = c
        else:
            a = c
        if b - a < 1e-14:
            break
    return 0.5 * (a + b)


if __name__ == "__main__":
    main()
