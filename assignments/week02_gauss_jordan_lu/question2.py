"""
question2.py — week02_gauss_jordan_lu
-------------------------------------
Problem : Check Gauss-Jordan on the three-variable system
              2y + 5z = 1,   3x - y + 2z = -2,   x - y + 3z = 3
          (augmented matrix in data/asgn2_mat1).
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    gauss_jordan_elimination_augmented,
    matrix_residual,
    read_matrix_from_file,
)

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"

    augmented = read_matrix_from_file(BASE_DIR / "data" / "asgn2_mat1")
    n = len(augmented)
    A = [row[:n] for row in augmented]
    b = [row[n] for row in augmented]

    x = gauss_jordan_elimination_augmented(augmented)
    res = matrix_residual(A, x, b)

    lines = [
        "Gauss-Jordan elimination (partial pivoting)",
        "System (from data/asgn2_mat1):",
        "  2y + 5z = 1",
        "  3x - y + 2z = -2",
        "  x - y + 3z = 3",
        "",
        "The values of x, y, z are:",
        f"  x = {x[0]:.10f}",
        f"  y = {x[1]:.10f}",
        f"  z = {x[2]:.10f}",
        f"  residual ||Ax - b||_2 = {res:.3e}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
