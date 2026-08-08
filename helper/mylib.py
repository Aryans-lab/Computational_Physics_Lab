"""
mylib.py
--------
Shared numerical library for PHY341/745 – Physics Computer Lab, NISER.

This file contains ALL numerical routines for the course.  The front code
(main.py) only performs I/O and delegates computation to functions here.
System-built routines (NumPy, SciPy, etc.) are NOT used for any core
algorithm, in compliance with course rules.

Author  : Aryan Bandyopadhyay
Roll No.: 2411014
"""

# ═══════════════════════════════════════════════════════════════════════════════
#  File I/O helpers
# ═══════════════════════════════════════════════════════════════════════════════

def read_tokens(filename: str) -> list:
    """
    Read an ASCII file and return all whitespace-separated tokens as strings.
    """
    with open(filename, "r") as f:
        return f.read().split()


def parse_numbers(tokens: list) -> list:
    """
    Convert a list of string tokens to floats.

    Auto-detection rule: if the very first token is an integer N and the
    total token count is exactly N + 1, treat the first token as a count
    header and skip it.  Otherwise parse every token.
    """
    if not tokens:
        return []

    try:
        n = int(tokens[0])
    except ValueError:
        return [float(x) for x in tokens]

    if len(tokens) == n + 1:
        return [float(x) for x in tokens[1:]]

    return [float(x) for x in tokens]


def read_second_column(filename: str) -> list:
    """
    Parse a two-column ASCII file of the form:
        <index>  <value>
    and return a list of the second-column values.
    """
    values = []
    with open(filename, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                values.append(float(parts[1]))
    return values


# ═══════════════════════════════════════════════════════════════════════════════
#  Scalar arithmetic
# ═══════════════════════════════════════════════════════════════════════════════

def sum_numbers(values: list) -> float:
    """Sum a list of floats using an explicit loop (no sum())."""
    total = 0.0
    for x in values:
        total += x
    return total


# ═══════════════════════════════════════════════════════════════════════════════
#  Matrix I/O parsers
# ═══════════════════════════════════════════════════════════════════════════════

def read_matrix_add(tokens: list) -> tuple:
    """
    Parse tokens for a matrix-addition problem.

    Expected token layout:
        n m               <- row and column count (both matrices are n x m)
        <n*m floats>      <- matrix A, row-major
        <n*m floats>      <- matrix B, row-major

    Returns: (n, m, A, B)
    """
    n = int(tokens[0])
    m = int(tokens[1])
    pos = 2

    A = []
    for i in range(n):
        row = []
        for j in range(m):
            row.append(float(tokens[pos]))
            pos += 1
        A.append(row)

    B = []
    for i in range(n):
        row = []
        for j in range(m):
            row.append(float(tokens[pos]))
            pos += 1
        B.append(row)

    return n, m, A, B


def read_matrix_multiply(tokens: list) -> tuple:
    """
    Parse tokens for a matrix-multiplication problem.

    Expected token layout:
        n p m             <- A is n x p, B is p x m
        <n*p floats>      <- matrix A, row-major
        <p*m floats>      <- matrix B, row-major

    Returns: (n, p, m, A, B)
    """
    n = int(tokens[0])
    p = int(tokens[1])
    m = int(tokens[2])
    pos = 3

    A = []
    for i in range(n):
        row = []
        for j in range(p):
            row.append(float(tokens[pos]))
            pos += 1
        A.append(row)

    B = []
    for i in range(p):
        row = []
        for j in range(m):
            row.append(float(tokens[pos]))
            pos += 1
        B.append(row)

    return n, p, m, A, B


# ═══════════════════════════════════════════════════════════════════════════════
#  Matrix operations
# ═══════════════════════════════════════════════════════════════════════════════

def matrix_add(A: list, B: list, n: int, m: int) -> list:
    """Return element-wise sum C = A + B for n x m matrices."""
    C = []
    for i in range(n):
        row = []
        for j in range(m):
            row.append(A[i][j] + B[i][j])
        C.append(row)
    return C


def matrix_multiply(A: list, B: list, n: int, p: int, m: int) -> list:
    """
    Multiply A (n x p) by B (p x m) using the standard triple loop.
    Returns C (n x m).
    """
    # initialise result matrix to zero using a loop (course requirement)
    C = []
    for i in range(n):
        C.append([0.0 for _ in range(m)])

    for i in range(n):
        for j in range(m):
            s = 0.0
            for k in range(p):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C


# ═══════════════════════════════════════════════════════════════════════════════
#  Output formatting
# ═══════════════════════════════════════════════════════════════════════════════

def matrix_to_string(C: list) -> str:
    """
    Convert a 2-D list to a whitespace-separated string suitable for file
    output.  Each row is on its own line; values are formatted to 6 d.p.
    """
    lines = []
    for row in C:
        lines.append("  ".join(f"{x:.6f}" for x in row))
    return "\n".join(lines) + "\n"
