"""
question3.py — week04_jacobi_gs
-------------------------------
Problem : Solve the 10x10 tridiagonal (Poisson-type) model problem
              -x_{i-1} + 2 x_i - x_{i+1} = 1,  b = (1, ..., 1)
          with Gauss-Seidel and with SOR (omega = 1.57), and compare the
          convergence rates.  Exact solution: (5, 9, 12, 14, 15, 15, 14, 12, 9, 5).
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    gauss_seidel,
    matrix_residual,
    sor_gauss_seidel,
)

BASE_DIR = Path(__file__).resolve().parent
ORDER = 10
OMEGA = 1.57
TOL = 1e-12


def build_tridiagonal(n: int) -> tuple[list[list[float]], list[float]]:
    """Poisson model: 2 on the diagonal, -1 on the off-diagonals, b = ones."""
    A = [[2.0 if i == j else (-1.0 if abs(i - j) == 1 else 0.0)
          for j in range(n)] for i in range(n)]
    b = [1.0] * n
    return A, b


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"

    A, b = build_tridiagonal(ORDER)
    x_gs, it_gs = gauss_seidel(A, b, tol=TOL, max_iterations=10000)
    x_sor, it_sor = sor_gauss_seidel(A, b, tol=TOL, max_iterations=10000, omega=OMEGA)
    x_exact = [5.0, 9.0, 12.0, 14.0, 15.0, 15.0, 14.0, 12.0, 9.0, 5.0]

    lines = [
        f"Gauss-Seidel vs SOR on the {ORDER}x{ORDER} Poisson model problem",
        f"(tol = {TOL:g})",
        "",
        f"Gauss-Seidel (omega = 1)      : converged after {it_gs:>5d} iterations",
        "  x = [" + ", ".join(f"{v:.6f}" for v in x_gs) + "]",
        f"  residual ||Ax - b||_2 = {matrix_residual(A, x_gs, b):.3e}",
        "",
        f"SOR (omega = {OMEGA})         : converged after {it_sor:>5d} iterations",
        "  x = [" + ", ".join(f"{v:.6f}" for v in x_sor) + "]",
        f"  residual ||Ax - b||_2 = {matrix_residual(A, x_sor, b):.3e}",
        "",
        f"speed-up: Gauss-Seidel/SOR = {it_gs / it_sor:.2f}x iterations",
        f"exact solution  x = {x_exact}",
        "The tridiagonal matrix is only *weakly* diagonally dominant, so the",
        "spectral radius sits close to 1 and plain Gauss-Seidel crawls; the",
        "over-relaxation with omega = 1.57 (the optimal value for this",
        "geometry is 2/(1+sin(pi/11)) ~= 1.56) cuts the iteration count by a",
        "factor of about 6.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
