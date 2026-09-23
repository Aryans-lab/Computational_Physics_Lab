"""
question2.py — week04_jacobi_gs
-------------------------------
Problem : Solve the 6x6 system in data/asgn4_mat1 (b in data/asgn4_vec1)
          with the Jacobi and the Gauss-Seidel methods and compare.
Usage   : python question2.py [output_file]
          Defaults: output/q2_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    gauss_seidel,
    jacobi_it,
    matrix_residual,
    read_matrix_from_file,
    read_vector_from_file,
)

BASE_DIR = Path(__file__).resolve().parent
TOL = 1e-6


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q2_output.txt"

    A = read_matrix_from_file(BASE_DIR / "data" / "asgn4_mat1")
    b = read_vector_from_file(BASE_DIR / "data" / "asgn4_vec1")

    x_j, it_j = jacobi_it(A, b, tol=TOL)
    x_g, it_g = gauss_seidel(A, b, tol=TOL)

    lines = [
        f"Jacobi and Gauss-Seidel on the 6x6 system (tol = {TOL:g})",
        "",
        f"Jacobi converged after {it_j} iterations:",
        "  x = [" + ", ".join(f"{v:.10f}" for v in x_j) + "]",
        f"  residual ||Ax - b||_2 = {matrix_residual(A, x_j, b):.3e}",
        "",
        f"Gauss-Seidel converged after {it_g} iterations:",
        "  x = [" + ", ".join(f"{v:.10f}" for v in x_g) + "]",
        f"  residual ||Ax - b||_2 = {matrix_residual(A, x_g, b):.3e}",
        "",
        f"speed-up: Jacobi/Gauss-Seidel = {it_j / it_g:.2f}x iterations",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
