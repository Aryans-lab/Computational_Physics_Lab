"""
question3.py — week02_gauss_jordan_lu
-------------------------------------
Problem : Find the LU decomposition of the matrix
              A = [[1, 2, 4], [3, 8, 14], [2, 6, 13]]   (data/asgn2_mat2)
          and verify it by reassembling L * U.
Usage   : python question3.py [output_file]
          Defaults: output/q3_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import lu_decomposition, read_matrix_from_file

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q3_output.txt"

    A = read_matrix_from_file(BASE_DIR / "data" / "asgn2_mat2")
    n = len(A)
    L, U, perm = lu_decomposition(A, return_perm=True)

    # Reassembly: L U must equal P A (permuted A), not A itself, because a
    # pivot swap was needed (|3| > |1| in the first column).
    PA = [[A[perm[i]][j] for j in range(n)] for i in range(n)]
    reassembled = [[sum(L[i][k] * U[k][j] for k in range(n)) for j in range(n)]
                   for i in range(n)]
    max_err = max(abs(reassembled[i][j] - PA[i][j])
                  for i in range(n) for j in range(n))

    lines = [
        "LU decomposition (Doolittle, partial pivoting):  P A = L U",
        "A (from data/asgn2_mat2):",
        "  " + ", ".join(f"[{', '.join(f'{v:g}' for v in row)}]" for row in A) + "  (rows printed)",
        "",
        f"Row permutation (pivot swap in column 1: 3 > 1 > 2): perm = {perm}",
        "",
        "Matrix L (unit lower triangular):",
    ]
    lines += ["  [" + ", ".join(f"{v: .6f}" for v in row) + "]" for row in L]
    lines.append("Matrix U (upper triangular):")
    lines += ["  [" + ", ".join(f"{v: .6f}" for v in row) + "]" for row in U]
    lines += [
        "",
        "Verification B = L * U:",
    ]
    lines += ["  [" + ", ".join(f"{v: .6f}" for v in row) + "]" for row in reassembled]
    lines += [
        "",
        "P * A (the permuted original):",
    ]
    lines += ["  [" + ", ".join(f"{v: .6f}" for v in row) + "]" for row in PA]
    lines += [
        f"max |L U - P A| = {max_err:.3e}   ->  decomposition verified",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
