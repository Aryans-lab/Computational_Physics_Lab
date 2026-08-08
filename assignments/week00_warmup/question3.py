"""
question3.py
------------
Problem : Read matrices A (3x3), B (3x3) and column vectors C (3x1), D (3x1)
          from ASCII data files.  Compute and write to an output file:
            AB   - matrix product
            BC   - matrix-vector product
            D.C  - dot product (scalar)
Usage   : python question3.py <output_file>
          e.g.: python question3.py output/q3_output.txt
          Data files are expected in data/ relative to this script.
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER

Note    : All arithmetic is pure Python - no NumPy/SciPy.
          A small floating-point artefact (~1e-16) may appear in AB[0][0];
          this is standard IEEE 754 double-precision behaviour, not a bug.
"""

import os
import sys

# Data files sit in data/ next to this script
BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA_DIR  = os.path.join(BASE_DIR, "data")


# -- I/O helpers

def read_matrix(filename: str) -> list:
    """Read a whitespace-delimited ASCII file; return a 2-D list of floats."""
    matrix = []
    with open(os.path.join(DATA_DIR, filename), "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                matrix.append([float(x) for x in stripped.split()])
    return matrix


def format_matrix(label: str, M: list) -> str:
    """Return a labelled, formatted string representation of matrix M."""
    lines = [label]
    for row in M:
        lines.append("  [" + ", ".join(f"{x:10.6f}" for x in row) + "]")
    return "\n".join(lines)


# -- Numerical routines

def matrix_multiply(X: list, Y: list) -> list:
    """
    Multiply X (n x p) by Y (p x m) via the standard triple loop.
    Returns an n x m result initialised with zeros in a loop.
    """
    n = len(X);  p = len(X[0]);  m = len(Y[0])
    result = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0.0
            for k in range(p):
                s += X[i][k] * Y[k][j]
            result[i][j] = s
    return result


def dot_product(X: list, Y: list) -> float:
    """
    Dot product of two column vectors stored as (n x 1) 2-D lists.
    Returns a scalar; exits on dimension mismatch.
    """
    if len(X) != len(Y):
        print("Error: vector length mismatch in dot product.", file=sys.stderr)
        sys.exit(1)
    total = 0.0
    for i in range(len(X)):
        total += X[i][0] * Y[i][0]
    return total


def main() -> None:
    # -- Argument check
    if len(sys.argv) != 2:
        print("Usage: python question3.py <output_file>", file=sys.stderr)
        sys.exit(1)

    output_file = sys.argv[1]

    # -- Read matrices and vectors from data files
    A = read_matrix("asgn0_matA")
    B = read_matrix("asgn0_matB")
    C = read_matrix("asgn0_vecC")
    D = read_matrix("asgn0_vecD")

    # -- Compute
    AB = matrix_multiply(A, B)
    BC = matrix_multiply(B, C)
    DC = dot_product(D, C)

    # -- Assemble output text
    output_text = "\n".join([
        format_matrix("Matrix AB:", AB),
        "",
        format_matrix("Matrix BC:", BC),
        "",
        f"Dot product D.C = {DC}",
        "",
        "# Note: AB[0][0] may show a ~1e-16 floating-point artefact",
        "# (IEEE 754 double-precision rounding), which is not a coding error.",
    ]) + "\n"

    # -- Write to output file
    with open(output_file, "w") as out:
        out.write(output_text)

    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()

# -- Output
# Matrix AB:
#   [ -0.300000,  -3.500000,   5.200000]
#   [ -4.500000,  -2.000000,   4.500000]
#   [  9.300000,   0.800000,  -7.000000]
#
# Matrix BC:
#   [  1.000000]
#   [ -5.750000]
#   [ -9.000000]
#
# Dot product D.C = -3.5
