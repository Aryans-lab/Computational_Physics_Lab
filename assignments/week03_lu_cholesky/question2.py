"""
question2.py — week03_lu_cholesky
---------------------------------
Problem : Solve the 6x6 system in data/asgn3_mat1 with right-hand side
          data/asgn3_vec1 using LU decomposition with forward-backward
          substitution (Doolittle scheme, partial pivoting).
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    lu_forback,
    matrix_residual,
    read_matrix_from_file,
    read_vector_from_file,
)

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"

    A = read_matrix_from_file(BASE_DIR / "data" / "asgn3_mat1")
    b = read_vector_from_file(BASE_DIR / "data" / "asgn3_vec1")

    x = lu_forback(A, b)  # Doolittle LU + partial pivoting + forward/backward
    res = matrix_residual(A, x, b)

    lines = [
        "LU decomposition (Doolittle, partial pivoting) + forward-backward",
        "substitution for the 6x6 system (data/asgn3_mat1, data/asgn3_vec1).",
        "",
        "The solution of the system of linear equations is:",
    ]
    for i, xi in enumerate(x):
        lines.append(f"  a{i + 1} = {xi:.10f}")
    lines += [
        "",
        f"residual ||Ax - b||_2 = {res:.3e}",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
