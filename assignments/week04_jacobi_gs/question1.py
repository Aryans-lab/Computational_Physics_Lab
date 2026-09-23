"""
question1.py — week04_jacobi_gs
-------------------------------
Problem : Prepare library routines for the Jacobi and Gauss-Seidel iterative
          solvers (and their SOR generalisation).  The routines live in
          mylib.py (jacobi_it, gauss_seidel, sor_gauss_seidel); this driver
          benchmarks all three on the diagonally-dominant 6x6 system in
          data/asgn4_mat1 and reports iteration counts and residuals.
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    gauss_seidel,
    is_diagonally_dominant,
    jacobi_it,
    matrix_residual,
    read_matrix_from_file,
    read_vector_from_file,
    sor_gauss_seidel,
)

BASE_DIR = Path(__file__).resolve().parent
TOL = 1e-8


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"

    A = read_matrix_from_file(BASE_DIR / "data" / "asgn4_mat1")
    b = read_vector_from_file(BASE_DIR / "data" / "asgn4_vec1")
    strict_dd = is_diagonally_dominant(A, strict=True)

    x_j, it_j = jacobi_it(A, b, tol=TOL)
    x_g, it_g = gauss_seidel(A, b, tol=TOL)
    x_s, it_s = sor_gauss_seidel(A, b, tol=TOL, omega=1.3)

    lines = [
        "Iterative solvers benchmark (tol = 1e-8, max |dx| stopping test)",
        "System: 6x6 from data/asgn4_mat1 / data/asgn4_vec1",
        f"strictly diagonally dominant: {strict_dd}  "
        "(guarantees convergence of Jacobi and Gauss-Seidel)",
        "",
        f"{'method':<22} {'iterations':>10}  {'residual ||Ax-b||_2':>20}  solution",
    ]
    for name, x, it in (("Jacobi", x_j, it_j),
                        ("Gauss-Seidel", x_g, it_g),
                        ("SOR (omega = 1.3)", x_s, it_s)):
        lines.append(
            f"{name:<22} {it:>10d}  {matrix_residual(A, x, b):>20.3e}  "
            + "[" + ", ".join(f"{v:.6f}" for v in x) + "]"
        )
    lines += [
        "",
        "Notes:",
        "  * Jacobi updates every component from the previous iterate only,",
        "    so the sweep is embarrassingly parallel.",
        "  * Gauss-Seidel reuses fresh values in-place: typically ~2x fewer",
        "    iterations on this system.",
        "  * SOR with 1 < omega < 2 accelerates further; omega is tuned per",
        "    problem (question3.py studies omega = 1.57 on a Poisson model).",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
