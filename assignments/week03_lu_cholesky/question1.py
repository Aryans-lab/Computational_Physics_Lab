"""
question1.py — week03_lu_cholesky
---------------------------------
Problem : Prepare a library routine for Cholesky factorisation A = L L^T that
          includes an explicit check for symmetric matrices.  The routine
          lives in mylib.py (cholesky_decomposition); this driver verifies it
          on the symmetric positive-definite matrix in data/asgn3_mat2 and
          demonstrates both failure modes (non-symmetric, non-positive-
          definite input).
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import (
    cholesky_decomposition,
    is_symmetric,
    read_matrix_from_file,
)

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"

    A = read_matrix_from_file(BASE_DIR / "data" / "asgn3_mat2")
    n = len(A)

    lines = [
        "Cholesky factorisation A = L L^T (with symmetry check)",
        f"A (from data/asgn3_mat2), symmetric: {is_symmetric(A)}",
        "  [" + ", ".join(f"{v:g}" for v in A[0]) + "]",
        "  [" + ", ".join(f"{v:g}" for v in A[1]) + "]",
        "  [" + ", ".join(f"{v:g}" for v in A[2]) + "]",
        "  [" + ", ".join(f"{v:g}" for v in A[3]) + "]",
        "",
        "L (lower triangular, positive diagonal):",
    ]
    L = cholesky_decomposition(A)
    lines += ["  [" + ", ".join(f"{v: .10f}" for v in row) + "]" for row in L]

    # Verification: L L^T must reproduce A.
    reassembled = [[sum(L[i][k] * L[j][k] for k in range(n)) for j in range(n)]
                   for i in range(n)]
    max_err = max(abs(reassembled[i][j] - A[i][j]) for i in range(n) for j in range(n))
    lines += [
        "",
        "Verification  ||L L^T - A||_max = "
        f"{max_err:.3e}   ->  factorisation exact to round-off",
        "",
        "Guard rails (both must raise ValueError):",
    ]
    try:
        cholesky_decomposition([[1.0, 2.0], [3.0, 4.0]])
        lines.append("  non-symmetric input : NO ERROR (unexpected)")
    except ValueError as e:
        lines.append(f"  non-symmetric input : ValueError -> {e}")
    try:
        cholesky_decomposition([[-1.0, 0.0], [0.0, 1.0]])
        lines.append("  non-positive-definite input : NO ERROR (unexpected)")
    except ValueError as e:
        lines.append(f"  non-positive-definite input : ValueError -> {e}")

    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
