"""
question1.py — week02_gauss_jordan_lu
-------------------------------------
Problem : Prepare library routines for Gauss-Jordan elimination (with row
          pivoting) and LU decomposition.  Choice of LU scheme: **Doolittle**
          (unit lower-triangular L).  The routines live in mylib.py
          (gauss_jordan_elimination_augmented, lu_decomposition, lu_forback);
          this driver exercises them, including a system that forces a pivot
          swap in a *later* column (the case naive non-pivoted codes break on).
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    gauss_jordan_elimination_augmented,
    lu_decomposition,
    lu_forback,
    matrix_residual,
    print_matrix,
)

BASE_DIR = Path(__file__).resolve().parent

# A 3x3 system whose second column has a zero on the diagonal without a
# row exchange: A = [[1, 1, 0], [0, 0, 1], [0, 1, 0]],  b = [1, 2, 3].
# Exact solution: x = (-2, 3, 2).  A non-pivoted factorisation fails at step 2.
Pivot_TEST_A = [[1.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
Pivot_TEST_b = [1.0, 2.0, 3.0]


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"
    lines: list[str] = []

    # -- 1. Gauss-Jordan on the pivot-stress system
    aug = [Pivot_TEST_A[i][:] + [Pivot_TEST_b[i]] for i in range(3)]
    x = gauss_jordan_elimination_augmented(aug)
    lines += [
        "Gauss-Jordan elimination (partial pivoting) — pivot-stress system",
        "  A = [[1, 1, 0], [0, 0, 1], [0, 1, 0]],  b = [1, 2, 3]",
        "  (column 2 has a zero pivot unless rows 2 and 3 are exchanged)",
        f"  solution x = {x}",
        f"  residual ||Ax - b||_2 = {matrix_residual(Pivot_TEST_A, x, Pivot_TEST_b):.3e}",
        "",
    ]

    # -- 2. Doolittle LU with pivoting, verified by reassembly
    A = [[1.0, 2.0, 4.0], [3.0, 8.0, 14.0], [2.0, 6.0, 13.0]]
    L, U, perm = lu_decomposition(A, return_perm=True)
    lines += [
        "LU decomposition — Doolittle scheme, partial pivoting (P A = L U)",
        "  A = [[1, 2, 4], [3, 8, 14], [2, 6, 13]]",
        "  row permutation applied: perm = " + str(perm) + "  (first column: 3 > 1 > 2)",
        "  L =",
    ]
    lines += ["    " + "  ".join(f"{v:.6f}" for v in row) for row in L]
    lines += ["  U ="]
    lines += ["    " + "  ".join(f"{v:.6f}" for v in row) for row in U]

    # Reassembly check: L @ U must equal P @ A.
    LU = [[sum(L[i][k] * U[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    PA = [[A[perm[i]][j] for j in range(3)] for i in range(3)]
    max_err = max(abs(LU[i][j] - PA[i][j]) for i in range(3) for j in range(3))
    lines += [
        f"  reassembly check  ||L U - P A||_max = {max_err:.3e}",
        "",
    ]

    # -- 3. Forward-backward substitution on the same matrix
    b = [1.0, 2.0, 3.0]
    x2 = lu_forback(A, b)
    lines += [
        "Forward-backward substitution on  A x = b,  b = [1, 2, 3]:",
        f"  solution x = {x2}",
        f"  residual ||Ax - b||_2 = {matrix_residual(A, x2, b):.3e}",
        "  (cross-check: gauss_jordan on [A|b] gives the same vector)",
    ]
    x3 = gauss_jordan_elimination_augmented([A[i][:] + [b[i]] for i in range(3)])
    agree = max(abs(x2[i] - x3[i]) for i in range(3))
    lines.append(f"  max |x_LU - x_GJ| = {agree:.3e}")

    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")
    print_matrix(x2, label="solution x (LU)")


if __name__ == "__main__":
    main()
