"""
question1.py — week05a_bisection_regf
-------------------------------------
Problem : Code up the library functions Bisection and Regula Falsi.  The
          routines live in mylib.py (bisection, regula_falsi, with
          bracket_root for automatic bracket expansion); this driver
          validates both on f(x) = x^2 - 2 against the known root sqrt(2),
          reporting iteration counts and achieved accuracy.
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from mylib import bisection, bracket_root, regula_falsi

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"

    f = lambda x: x * x - 2.0  # noqa: E731
    a, b = 1.0, 2.0
    true_root = math.sqrt(2.0)
    acc = 1e-10

    r_b, it_b = bisection(f, a, b, accuracy=acc)
    r_r, it_r = regula_falsi(f, a, b, accuracy=acc)

    # Auto-bracketing demo: start from a bracket that does NOT contain the
    # root and let bracket_root expand it.
    a2, b2 = bracket_root(f, 3.0, 4.0)

    lines = [
        "Validation of the bracketing routines on f(x) = x^2 - 2,  root = sqrt(2)",
        f"(bracket [{a}, {b}], accuracy = {acc:g})",
        "",
        f"bisection     : root = {r_b:.12f}  in {it_b:3d} iterations",
        f"                |f(root)| = {abs(f(r_b)):.3e},   error = {abs(r_b - true_root):.3e}",
        f"regula falsi  : root = {r_r:.12f}  in {it_r:3d} iterations",
        f"                |f(root)| = {abs(f(r_r)):.3e},   error = {abs(r_r - true_root):.3e}",
        "",
        "Convergence-rate check (bisection halves the interval each step,):",
        f"  iters for 1e-{10} from a 1-wide bracket ~ log2(1/1e-10) = {math.log2(1e10):.1f} -> {it_b}",
        "",
        f"Auto-bracketing demo: start from [{3.0}, {4.0}] (no sign change)",
        f"  bracket_root expanded to [{a2:.6f}, {b2:.6f}]",
        f"  bisection on the expanded bracket: root = {bisection(f, a2, b2, accuracy=acc)[0]:.12f}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
