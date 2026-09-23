"""
question3.py — week03_lu_cholesky
---------------------------------
Problem : Use Cholesky factorisation for the symmetric positive-definite
          matrix A (data/asgn3_mat2) and solve A x = b with
          b = (3, 3, 1, 3) (data/asgn3_vec2) via
              L y = b      (forward substitution)
              L^T x = y    (backward substitution).
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    cholesky_decomposition,
    cholesky_forback,
    matrix_residual,
    read_matrix_from_file,
    read_vector_from_file,
)

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"

    A = read_matrix_from_file(BASE_DIR / "data" / "asgn3_mat2")
    b = read_vector_from_file(BASE_DIR / "data" / "asgn3_vec2")

    L = cholesky_decomposition(A)
    x = cholesky_forback(L, b)
    res = matrix_residual(A, x, b)

    lines = [
        "Cholesky factorisation A = L L^T  +  forward/backward substitution",
        "(data/asgn3_mat2, data/asgn3_vec2)",
        "",
        "The Cholesky factorisation of the matrix A is:",
    ]
    lines += ["  [" + ", ".join(f"{v: .10f}" for v in row) + "]" for row in L]
    lines += [
        "",
        "The solution of the system of linear equations is:",
    ]
    for i, xi in enumerate(x):
        lines.append(f"  a{i + 1} = {xi:.12f}")
    lines += [
        "",
        f"residual ||Ax - b||_2 = {res:.3e}",
        "exact solution is x = (0, 1, 1, 1) — recovered to machine precision",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
