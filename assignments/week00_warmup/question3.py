"""
question3.py
------------
Problem : Read matrices A (3x3), B (3x3) and column vectors C (3x1), D (3x1)
          from ASCII data files in the data/ subdirectory.
          Compute and print:
            AB   - matrix product
            BC   - matrix-vector product
            D.C  - dot product (scalar)
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER

Note    : All arithmetic is pure Python - no NumPy/SciPy.
          A small floating-point artefact (~1e-16) may appear in AB[0][0];
          this is standard IEEE 754 double-precision behaviour, not a bug.
"""

import os
import sys

# Path to the data/ folder sitting next to this script
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


# -- I/O helpers

def read_matrix(filename: str) -> list:
    """Read a whitespace-delimited ASCII file and return a 2-D list of floats."""
    filepath = os.path.join(BASE_DIR, filename)
    matrix = []
    with open(filepath, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped:                          # skip blank lines
                row = [float(x) for x in stripped.split()]
                matrix.append(row)
    return matrix


def print_matrix(label: str, M: list) -> None:
    """Pretty-print a matrix with a descriptive label."""
    print(label)
    for row in M:
        print("  [" + ", ".join(f"{x:10.6f}" for x in row) + "]")


# -- Numerical routines

def matrix_multiply(X: list, Y: list) -> list:
    """
    Multiply matrix X (n x p) by matrix Y (p x m) using the standard
    O(n*p*m) triple loop.  Returns an n x m result matrix.
    """
    n = len(X);   p = len(X[0]);   m = len(Y[0])
    # initialise result with zeros (loop-based, as required by course rules)
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
    Returns a scalar.  Exits with an error if dimensions mismatch.
    """
    if len(X) != len(Y):
        print("Error: vectors must have equal length for a dot product.",
              file=sys.stderr)
        sys.exit(1)
    total = 0.0
    for i in range(len(X)):
        total += X[i][0] * Y[i][0]
    return total


# -- Main computation

A = read_matrix("asgn0_matA")
B = read_matrix("asgn0_matB")
C = read_matrix("asgn0_vecC")
D = read_matrix("asgn0_vecD")

AB = matrix_multiply(A, B)
BC = matrix_multiply(B, C)
DC = dot_product(D, C)

print_matrix("Matrix AB:", AB)
print_matrix("\nMatrix BC:", BC)
print(f"\nDot product D.C = {DC}")

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
