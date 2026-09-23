"""
mylib.py
========

A computational-physics library implemented **entirely from first principles**
in pure standard-library Python.  No NumPy, no SciPy, no third-party numeric
packages: every routine below — direct and iterative linear solvers, non-linear
root finders, polynomial root isolation, and quadrature — is written by hand so
that the floating-point mechanics (IEEE 754 round-off, pivoting, conditioning,
convergence orders) stay visible.

Design rules
------------
1. **Library / front-code separation.**  This module contains *only* numerical
   algorithms.  Each driver script (``assignments/*/question*.py``) performs
   input/output and calls into this file — the standard convention of the
   NISER PHY341/745 lab this work was developed in.
2. **Verification by residuals.**  Every solver is designed to be checked:
   ``matrix_residual(A, x, b)`` reports the backward error, and the driver
   scripts print it alongside each solution.
3. **Honest numerics.**  Partial pivoting wherever a pivot is chosen,
   scale-relative singularity tolerances (``~eps * n * max|A|``), no silent
   fallbacks.
4. **Pure functions.**  Algorithms return values (and iteration counts); they
   do not print.  Formatting belongs to the front code.

Module map
----------
1.  Pseudo-random number generation & distributions (LCG, uniform, exponential)
2.  Complex arithmetic (MyComplex)
3.  Matrix/vector utilities & file I/O
4.  Direct linear solvers (Gauss-Jordan, inversion, LU, determinants)
5.  Symmetric positive-definite systems (Cholesky)
6.  Iterative linear solvers (Jacobi, Gauss-Seidel, SOR)
7.  Non-linear root finding (bisection, regula falsi, fixed point, Newton)
8.  Polynomial roots & deflation (Horner, Laguerre, synthetic division)
9.  Numerical integration & quadrature (midpoint, trapezoidal, Simpson 1/3,
    Monte Carlo, Gauss-Legendre, Gauss-Laguerre)
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

__all__ = [
    # 1. random
    "LCG_A", "LCG_C", "LCG_M", "myrand", "myrand_reset",
    "random_uniform", "random_exponential",
    # 2. complex
    "MyComplex",
    # 3. linear algebra utilities / I/O
    "read_matrix_from_file", "read_vector_from_file", "write_matrix_to_file",
    "print_matrix", "matrix_multiply", "matrix_transpose", "matrix_add",
    "matrix_sub", "matrix_vector_multiply", "dot_product_vector",
    "vector_norm", "matrix_residual", "is_symmetric", "is_diagonally_dominant",
    # 4. direct solvers
    "gauss_jordan_elimination_augmented", "gauss_jordan_inverse",
    "lu_decomposition", "lu_forback", "matrix_determinant_lu",
    # 5. Cholesky
    "cholesky_decomposition", "cholesky_forback",
    # 6. iterative solvers
    "jacobi_it", "gauss_seidel", "sor_gauss_seidel",
    # 7. root finding
    "bracket_root", "bisection", "regula_falsi", "fixed_point",
    "finite_difference_derivative", "partial_derivative", "jacobian",
    "newton_raphson", "newton_raphson_system",
    # 8. polynomials
    "polynomial_value", "polynomial_first_derivative",
    "polynomial_second_derivative", "laguerre", "synthetic_division",
    "laguerre_roots",
    # 9. integration
    "midpoint", "trapezoidal", "simpson", "monte_carlo",
    "integration_error_bound_N", "GAUSS_LEGENDRE", "GAUSS_LAGUERRE",
    "gaussian_quadrature",
]

# Machine epsilon of IEEE 754 double precision.
EPS = 2.220446049250313e-16

Matrix = Sequence[Sequence[float]]
Vector = Sequence[float]
ScalarFunction = Callable[[float], float]


# =============================================================================
# 1. RANDOM NUMBER GENERATION & DISTRIBUTIONS
# =============================================================================

# GCC / C runtime LCG parameters.  m = 2^15, a = 1103515245, c = 12345.
# These satisfy the Hull-Dobell full-period conditions for a power-of-two
# modulus (a = 1 mod 4, c odd, m a power of 2), so the generator cycles
# through all m states before repeating.
LCG_A = 1103515245
LCG_C = 12345
LCG_M = 32768

_LCG_STATE = {"x0": 0}


def myrand(seed: int | None = None) -> float:
    """
    Linear Congruential Generator producing floats in [0, 1).

    Recurrence:  x_{i+1} = (a * x_i + c) mod m,   return x_{i+1} / m.

    Passing ``seed`` (re)initialises the internal state; omitting it draws the
    next value of the current sequence.  The state persists across calls, so a
    sequence can be built up over many calls — call ``myrand_reset(seed)`` to
    start a fresh, reproducible sequence.
    """
    if seed is not None:
        _LCG_STATE["x0"] = int(seed) % LCG_M
    x1 = (LCG_A * _LCG_STATE["x0"] + LCG_C) % LCG_M
    _LCG_STATE["x0"] = x1
    return x1 / LCG_M


def myrand_reset(seed: int = 0) -> None:
    """Explicitly reseed the LCG state without drawing a number."""
    _LCG_STATE["x0"] = int(seed) % LCG_M


def random_uniform(a: float = 0.0, b: float = 1.0) -> float:
    """Uniform deviate in [a, b) via linear scaling X = a + (b-a)*xi."""
    return a + (b - a) * myrand()


def random_exponential(lam: float = 1.0) -> float:
    """
    Exponential deviate with PDF q(y) = lam * exp(-lam y) via the
    inverse-transform method:  y = -(1/lam) * ln(u),  u in (0, 1].
    """
    u = myrand()
    while u <= 0.0:  # guard the (measure-zero) case u = 0
        u = myrand()
    return -(1.0 / lam) * math.log(u)


# =============================================================================
# 2. COMPLEX ARITHMETIC (WITHOUT THE BUILT-IN complex TYPE)
# =============================================================================

class MyComplex:
    """
    Minimal complex-number class.  Supports +, -, *, and |z| via operator
    overloading; operands may be MyComplex or plain numbers.
    """

    __slots__ = ("real", "imag")

    def __init__(self, real: float, imag: float = 0.0) -> None:
        self.real = float(real)
        self.imag = float(imag)

    def __add__(self, other: MyComplex | complex | float | int) -> MyComplex:
        if isinstance(other, MyComplex):
            return MyComplex(self.real + other.real, self.imag + other.imag)
        return MyComplex(self.real + other, self.imag)

    def __sub__(self, other: MyComplex | complex | float | int) -> MyComplex:
        if isinstance(other, MyComplex):
            return MyComplex(self.real - other.real, self.imag - other.imag)
        return MyComplex(self.real - other, self.imag)

    def __mul__(self, other: MyComplex | complex | float | int) -> MyComplex:
        if isinstance(other, MyComplex):
            return MyComplex(
                self.real * other.real - self.imag * other.imag,
                self.real * other.imag + self.imag * other.real,
            )
        return MyComplex(self.real * other, self.imag * other)

    def __abs__(self) -> float:
        return math.hypot(self.real, self.imag)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, MyComplex):
            return self.real == other.real and self.imag == other.imag
        if isinstance(other, complex):
            return self.real == other.real and self.imag == other.imag
        if isinstance(other, (int, float)):
            return self.imag == 0.0 and self.real == other
        return NotImplemented

    def __repr__(self) -> str:
        sign = "+" if self.imag >= 0 else "-"
        return f"({self.real} {sign} {abs(self.imag)}j)"

    def __str__(self) -> str:
        sign = "+" if self.imag >= 0 else "-"
        return f"{self.real} {sign} {abs(self.imag)}j"


# =============================================================================
# 3. MATRIX / VECTOR UTILITIES & FILE I/O
# =============================================================================

def read_matrix_from_file(filename: str) -> list[list[float]]:
    """
    Read a whitespace-delimited ASCII file of numbers into a 2-D list of
    floats.  Blank lines and lines starting with '#' are ignored.
    ``filename`` is used verbatim (CWD-relative or absolute).
    """
    path = str(filename)
    try:
        handle = open(path)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Data file not found: {path!r}. Pass a path relative to the "
            f"current working directory, or an absolute path."
        ) from exc

    matrix: list[list[float]] = []
    with handle as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                matrix.append([float(x) for x in line.split()])
    if not matrix:
        raise ValueError(f"Data file {path!r} contains no numbers.")
    return matrix


def read_vector_from_file(filename: str) -> list[float]:
    """
    Read an ASCII file into a 1-D list of floats.  Each line contributes its
    first number, so both 'one value per line' and 'row of one' files work.
    """
    rows = read_matrix_from_file(filename)
    return [row[0] for row in rows]


def write_matrix_to_file(matrix: Matrix, filename: str, precision: int = 6) -> None:
    """Write a 1-D or 2-D list of numbers to ``filename`` in aligned columns."""
    if not matrix:
        return
    with open(filename, "w") as f:
        if not isinstance(matrix[0], (list, tuple)):
            for val in matrix:
                f.write(f"{float(val):.{precision}f}\n")
        else:
            for row in matrix:
                f.write("  ".join(f"{float(v):12.{precision}f}" for v in row) + "\n")


def print_matrix(matrix: Matrix, label: str | None = None, precision: int = 6) -> None:
    """Pretty-print a 1-D or 2-D list of numbers to stdout."""
    if label:
        print(f"--- {label} ---")
    if not matrix:
        print("[]")
        return
    if not isinstance(matrix[0], (list, tuple)):
        print("[" + ", ".join(f"{float(x):.{precision}f}" for x in matrix) + "]")
        return
    for row in matrix:
        print("  [" + "  ".join(f"{float(x):12.{precision}f}" for x in row) + "  ]")


def matrix_multiply(X: Matrix, Y: Matrix) -> list[list[float]]:
    """Multiply X (n x m) by Y (m x p); returns the n x p product."""
    rows_x, cols_x = len(X), len(X[0])
    rows_y, cols_y = len(Y), len(Y[0])
    if cols_x != rows_y:
        raise ValueError(
            f"Dimension mismatch: ({rows_x}x{cols_x}) @ ({rows_y}x{cols_y})"
        )
    result = [[0.0] * cols_y for _ in range(rows_x)]
    for i in range(rows_x):
        for k in range(cols_x):
            xik = X[i][k]
            if xik == 0.0:
                continue
            for j in range(cols_y):
                result[i][j] += xik * Y[k][j]
    return result


def matrix_transpose(A: Matrix) -> list[list[float]]:
    """Return the transpose of A."""
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def matrix_add(A: Matrix, B: Matrix) -> list[list[float]]:
    """Element-wise sum of two same-shape matrices."""
    if len(A) != len(B) or any(len(r) != len(B[i]) for i, r in enumerate(A)):
        raise ValueError("Shape mismatch in matrix_add.")
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def matrix_sub(A: Matrix, B: Matrix) -> list[list[float]]:
    """Element-wise difference of two same-shape matrices."""
    if len(A) != len(B) or any(len(r) != len(B[i]) for i, r in enumerate(A)):
        raise ValueError("Shape mismatch in matrix_sub.")
    return [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def matrix_vector_multiply(A: Matrix, x: Vector) -> list[float]:
    """Matrix-vector product A @ x (A: n x m, x: length m)."""
    if len(A[0]) != len(x):
        raise ValueError("Dimension mismatch in matrix_vector_multiply.")
    return [sum(a * v for a, v in zip(row, x)) for row in A]


def dot_product_vector(X: Vector, Y: Vector) -> float:
    """Dot product of two equal-length 1-D vectors."""
    if len(X) != len(Y):
        raise ValueError("Vector length mismatch in dot_product_vector.")
    return sum(x * y for x, y in zip(X, Y))


def vector_norm(v: Vector, p: int = 2) -> float:
    """The p-norm of a 1-D vector (default: Euclidean)."""
    return sum(abs(x) ** p for x in v) ** (1.0 / p)


def matrix_residual(A: Matrix, x: Vector, b: Vector) -> float:
    """
    Backward error of a linear-system solution:  ||A x - b||_2.
    A correct solution of an n x n system lands at ~sqrt(n) * eps * ||A|| * ||x||.
    """
    diff = [sum(a * v for a, v in zip(row, x)) - b_i for row, b_i in zip(A, b)]
    return vector_norm(diff, 2)


def is_symmetric(A: Matrix, tol: float = 1e-9) -> bool:
    """True if A is square and A[i][j] == A[j][i] within ``tol`` for all i, j."""
    n = len(A)
    if any(len(row) != n for row in A):
        return False
    return all(abs(A[i][j] - A[j][i]) <= tol for i in range(n) for j in range(n))


def is_diagonally_dominant(A: Matrix, strict: bool = False) -> bool:
    """
    True if every row satisfies  |A_ii|  (>= or >)  sum_{j != i} |A_ij|.
    Strict diagonal dominance guarantees convergence of Jacobi and Gauss-Seidel.
    """
    n = len(A)
    for i in range(n):
        off = sum(abs(A[i][j]) for j in range(n) if j != i)
        if strict:
            if abs(A[i][i]) <= off:
                return False
        else:
            if abs(A[i][i]) < off:
                return False
    return True


def _pivot_tolerance(A: Matrix) -> float:
    """
    Scale-relative threshold for treating an entry as zero.
    Uses ~eps * n * max|A_ij| so the test is independent of the matrix scale.
    """
    if not A or not A[0]:
        return EPS
    n = len(A)
    max_entry = max(abs(v) for row in A for v in row)
    return max(16.0 * EPS * n * max_entry, 1e-300)


# =============================================================================
# 4. DIRECT LINEAR SOLVERS
# =============================================================================

def gauss_jordan_elimination_augmented(
    augmented_matrix: Matrix, tol: float | None = None
) -> list[float]:
    """
    Solve A x = b with Gauss-Jordan elimination on the augmented [A | b]
    (n x (n+1)), using partial pivoting.  Reduces to RREF; the solution is the
    last column.

    Raises
    ------
    ValueError
        * "Inconsistent"   — a zero row with a non-zero right-hand side.
        * "Dependent"      — rank-deficient but consistent (infinitely many).
        * "Singular"       — any other rank-deficient case.
    """
    if not augmented_matrix or not augmented_matrix[0]:
        return []
    M = [[float(v) for v in row] for row in augmented_matrix]
    n = len(M)
    if len(M[0]) != n + 1:
        raise ValueError("Augmented matrix must have shape n x (n+1).")
    if tol is None:
        tol = _pivot_tolerance(M)

    for i in range(n):
        # Partial pivoting: largest magnitude in column i, rows i..n-1.
        pivot_row = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[pivot_row][i]) <= tol:
            continue  # numerically zero column; handled by the rank check below
        if pivot_row != i:
            M[i], M[pivot_row] = M[pivot_row], M[i]
        pivot = M[i][i]
        for j in range(i, n + 1):
            M[i][j] /= pivot
        for k in range(n):
            if k != i:
                factor = M[k][i]
                if factor == 0.0:
                    continue
                for j in range(i, n + 1):
                    M[k][j] -= factor * M[i][j]

    # Rank / consistency check after elimination.
    for i in range(n):
        coeff_zero = all(abs(M[i][j]) <= tol for j in range(n))
        if coeff_zero and abs(M[i][n]) > tol:
            raise ValueError("Inconsistent system: no solution exists.")
        if coeff_zero:
            # Zero row with zero RHS: the system is underdetermined.
            if any(abs(M[i][j]) > tol or abs(M[k][j]) > tol
                   for k in range(n) for j in range(n)):
                pass
            raise ValueError(
                "Dependent (rank-deficient) system: infinitely many solutions."
            )
    return [M[i][n] for i in range(n)]


def gauss_jordan_inverse(A: Matrix, tol: float | None = None) -> list[list[float]]:
    """
    Invert a square matrix by Gauss-Jordan elimination on [A | I_n].
    Uses partial pivoting; raises ValueError if A is singular.
    """
    if not A or not A[0]:
        return []
    M = [[float(v) for v in row] for row in A]
    n = len(M)
    if any(len(row) != n for row in M):
        raise ValueError("Matrix must be square.")
    if tol is None:
        tol = _pivot_tolerance(M)

    augmented = [M[i][:] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for i in range(n):
        pivot_row = max(range(i, n), key=lambda r: abs(augmented[r][i]))
        if abs(augmented[pivot_row][i]) <= tol:
            raise ValueError("Matrix is singular; inverse does not exist.")
        if pivot_row != i:
            augmented[i], augmented[pivot_row] = augmented[pivot_row], augmented[i]
        pivot = augmented[i][i]
        for j in range(i, 2 * n):
            augmented[i][j] /= pivot
        for k in range(n):
            if k != i:
                factor = augmented[k][i]
                if factor == 0.0:
                    continue
                for j in range(i, 2 * n):
                    augmented[k][j] -= factor * augmented[i][j]
    return [row[n:] for row in augmented]


def lu_decomposition(
    A: Matrix, tol: float | None = None, return_perm: bool = False
):
    """
    Doolittle LU factorisation **with partial pivoting**:  P A = L U, where L
    is unit lower triangular and U upper triangular.

    Parameters
    ----------
    A           : n x n matrix.
    tol         : pivot threshold (defaults to a scale-relative value).
    return_perm : if True, also return the row-permutation list ``perm`` such
                  that perm[i] is the source row of row i of P A.

    Returns (L, U) or (L, U, perm).  Raises ValueError on a zero pivot.
    """
    if not A or not A[0]:
        return []
    M = [[float(v) for v in row] for row in A]
    n = len(M)
    if any(len(row) != n for row in M):
        raise ValueError("Matrix A is not square.")
    if tol is None:
        tol = _pivot_tolerance(M)

    perm = list(range(n))
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        # Partial pivoting on column i (rows i..n-1).
        pivot_row = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[pivot_row][i]) <= tol:
            raise ValueError(f"LU decomposition failed: zero pivot at column {i}.")
        if pivot_row != i:
            M[i], M[pivot_row] = M[pivot_row], M[i]
            # The multipliers in columns < i belong to the rows, not the
            # positions — swap them along with the rows (diagonals stay put).
            for k in range(i):
                L[i][k], L[pivot_row][k] = L[pivot_row][k], L[i][k]
            perm[i], perm[pivot_row] = perm[pivot_row], perm[i]
        # Doolittle: row i of U, then column i of L.
        for j in range(i, n):
            U[i][j] = M[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        for j in range(i + 1, n):
            L[j][i] = (M[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    if return_perm:
        return L, U, perm
    return L, U


def lu_forback(A, b, method: str = "lu") -> list[float]:
    """
    Solve A x = b via (partially pivoted) LU decomposition and
    forward/backward substitution:  P A = L U,  L y = P b,  U x = y.

    ``A`` and ``b`` may be nested lists or paths to ASCII data files.
    ``method='Cholesky'`` routes SPD systems through ``cholesky_decomposition``.
    """
    C = read_matrix_from_file(A) if isinstance(A, str) else A
    if isinstance(b, str):
        d_mat = read_matrix_from_file(b)
        bv = [row[0] for row in d_mat]
    elif isinstance(b[0], (list, tuple)):
        bv = [float(row[0]) for row in b]
    else:
        bv = [float(x) for x in b]
    n = len(C)

    if method == "Cholesky":
        L = cholesky_decomposition(C)
        return cholesky_forback(L, bv)

    L, U, perm = lu_decomposition(C, return_perm=True)
    pb = [bv[i] for i in perm]

    # Forward substitution: L y = P b  (unit diagonal of L).
    y = [0.0] * n
    for i in range(n):
        y[i] = pb[i] - sum(L[i][j] * y[j] for j in range(i))

    # Backward substitution: U x = y.
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return x


def matrix_determinant_lu(A: Matrix) -> float:
    """
    det(A) from partially pivoted LU:  det(A) = sign(P) * prod_i U_ii.
    The permutation sign is tracked, so the result is exact in sign.
    """
    L, U, perm = lu_decomposition(A, return_perm=True)
    # sign = (-1)^{n - #cycles}: a cycle of length L is L-1 transpositions,
    # so it flips the sign when L is even.  ``count`` below ends at L + 1.
    sign = 1
    visited = [False] * len(perm)
    for i in range(len(perm)):
        if visited[i]:
            continue
        count = 1
        j = i
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            count += 1
        if count % 2 == 1:
            sign = -sign
    det = sign
    for i in range(len(U)):
        det *= U[i][i]
    return det


# =============================================================================
# 5. SYMMETRIC POSITIVE-DEFINITE SYSTEMS: CHOLESKY
# =============================================================================

def cholesky_decomposition(A: Matrix, tol: float = 1e-9) -> list[list[float]]:
    """
    Factor a real symmetric positive-definite matrix into A = L L^T with L
    lower triangular and positive diagonal.  Verifies symmetry and
    positive-definiteness explicitly (raises ValueError otherwise).
    """
    if not A:
        return []
    A_mat = [[float(v) for v in row] for row in A]
    n = len(A_mat)
    if any(len(row) != n for row in A_mat):
        raise ValueError("Matrix A is not square.")
    for i in range(n):
        for j in range(i):
            if abs(A_mat[i][j] - A_mat[j][i]) > tol:
                raise ValueError(
                    "Matrix A is not symmetric; Cholesky requires symmetry."
                )

    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            val = A_mat[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if val <= 0.0:
                    raise ValueError(
                        f"Matrix A is not positive-definite "
                        f"(non-positive pivot at index {i})."
                    )
                L[i][j] = math.sqrt(val)
            else:
                L[i][j] = val / L[j][j]
    return L


def cholesky_forback(L: Matrix, b) -> list[float]:
    """
    Solve A x = b given its Cholesky factor (A = L L^T):
    forward-substitute L y = b, then back-substitute L^T x = y.
    ``b`` may be a column vector (list of [x]) or a 1-D list.
    """
    n = len(L)
    bv = [float(b[i][0]) if isinstance(b[i], (list, tuple)) else float(b[i])
          for i in range(n)]

    y = [0.0] * n
    for i in range(n):
        y[i] = (bv[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(L[j][i] * x[j] for j in range(i + 1, n))) / L[i][i]
    return x


# =============================================================================
# 6. ITERATIVE LINEAR SOLVERS: JACOBI, GAUSS-SEIDEL, SOR
# =============================================================================

def _flatten_b(b) -> list[float]:
    return [float(row[0]) if isinstance(row, (list, tuple)) else float(row) for row in b]


def jacobi_it(A: Matrix, b, tol: float = 1e-8, max_iterations: int = 1000):
    """
    Jacobi iteration:  x_i^{(k+1)} = (b_i - sum_{j != i} A_ij x_j^{(k)}) / A_ii.

    All components are updated from the *previous* iterate (synchronous
    updates — trivially parallelisable).  Converges when A is strictly
    diagonally dominant (check with ``is_diagonally_dominant(A, strict=True)``).

    Returns (x, iterations).  Raises ValueError on non-convergence.
    """
    if not A or not A[0] or not b:
        return [], 0
    A_mat = [[float(v) for v in row] for row in A]
    bv = _flatten_b(b)
    n = len(bv)
    if any(A_mat[i][i] == 0.0 for i in range(n)):
        raise ValueError("Jacobi requires a non-zero diagonal.")

    x = [0.0] * n
    for k in range(1, max_iterations + 1):
        x_prev = x[:]
        for i in range(n):
            s = sum(A_mat[i][j] * x_prev[j] for j in range(n) if j != i)
            x[i] = (bv[i] - s) / A_mat[i][i]
        if max(abs(x[i] - x_prev[i]) for i in range(n)) < tol:
            return x, k
    raise ValueError(f"Jacobi did not converge within {max_iterations} iterations.")


def gauss_seidel(A: Matrix, b, tol: float = 1e-8, max_iterations: int = 500):
    """
    Gauss-Seidel iteration: like Jacobi, but each update uses the freshest
    available values (in-place sweep).  Roughly twice as fast as Jacobi for
    typical systems; same diagonal-dominance convergence guarantee.

    Returns (x, iterations).  Raises ValueError on non-convergence.
    """
    if not A or not A[0] or not b:
        return [], 0
    A_mat = [[float(v) for v in row] for row in A]
    bv = _flatten_b(b)
    n = len(bv)
    if any(A_mat[i][i] == 0.0 for i in range(n)):
        raise ValueError("Gauss-Seidel requires a non-zero diagonal.")

    x = [0.0] * n
    for k in range(1, max_iterations + 1):
        x_prev = x[:]
        for i in range(n):
            s = sum(A_mat[i][j] * x[j] for j in range(i)) \
                + sum(A_mat[i][j] * x[j] for j in range(i + 1, n))
            x[i] = (bv[i] - s) / A_mat[i][i]
        if max(abs(x[i] - x_prev[i]) for i in range(n)) < tol:
            return x, k
    raise ValueError(f"Gauss-Seidel did not converge within {max_iterations} iterations.")


def sor_gauss_seidel(
    A: Matrix, b, tol: float = 1e-8, max_iterations: int = 500, omega: float = 1.0
):
    """
    Successive Over-Relaxation:
        x_i^{(k+1)} = (1 - w) x_i^{(k)} + (w / A_ii) (b_i - sum_{j != i} A_ij x_j).
    w = 1 recovers Gauss-Seidel; 1 < w < 2 over-relaxes (faster for smooth
    problems, e.g. Poisson); 0 < w < 1 under-relaxes (stabilises oscillation).

    Returns (x, iterations).  Raises ValueError on non-convergence.
    """
    if not (0.0 < omega < 2.0):
        raise ValueError("SOR requires 0 < omega < 2.")
    if not A or not A[0] or not b:
        return [], 0
    A_mat = [[float(v) for v in row] for row in A]
    bv = _flatten_b(b)
    n = len(bv)
    if any(A_mat[i][i] == 0.0 for i in range(n)):
        raise ValueError("SOR requires a non-zero diagonal.")

    x = [0.0] * n
    for k in range(1, max_iterations + 1):
        x_prev = x[:]
        for i in range(n):
            s = sum(A_mat[i][j] * x[j] for j in range(n) if j != i)
            x_gs = (bv[i] - s) / A_mat[i][i]
            x[i] = (1.0 - omega) * x_prev[i] + omega * x_gs
        if max(abs(x[i] - x_prev[i]) for i in range(n)) < tol:
            return x, k
    raise ValueError(f"SOR did not converge within {max_iterations} iterations.")


# =============================================================================
# 7. ROOT FINDING FOR NON-LINEAR EQUATIONS
# =============================================================================

def bracket_root(f: ScalarFunction, a: float, b: float, max_iter: int = 100) -> tuple[float, float]:
    """
    Expand [a, b] until f(a) and f(b) have opposite signs (a valid bracket for
    the bisection family).  The side with the smaller |f| is pushed out, the
    step size growing by 10% each expansion.  Returns the bracketed (a, b).
    """
    m = 0.5
    for _ in range(max_iter):
        fa, fb = f(a), f(b)
        if fa * fb < 0.0:
            return a, b
        if abs(fa) <= abs(fb):
            a -= m * (b - a)
        else:
            b += m * (b - a)
        m += 0.1
    raise ValueError(f"Failed to bracket a root within {max_iter} expansions.")


def bisection(
    f: ScalarFunction, a: float, b: float,
    accuracy: float = 1e-8, max_iterations: int = 100,
) -> tuple[float, int]:
    """
    Bisection on [a, b]: guaranteed linear convergence (interval halves each
    step).  If the endpoints do not bracket a sign change, the bracket is
    expanded automatically (see ``bracket_root``).

    Returns (root, iterations).
    """
    if f(a) * f(b) >= 0.0:
        a, b = bracket_root(f, a, b)
    for i in range(1, max_iterations + 1):
        c = 0.5 * (a + b)
        fc = f(c)
        if abs(fc) < accuracy and abs(b - a) < accuracy:
            return c, i
        if fc * f(a) < 0.0:
            b = c
        else:
            a = c
    return 0.5 * (a + b), max_iterations


def regula_falsi(
    f: ScalarFunction, a: float, b: float,
    accuracy: float = 1e-8, max_iterations: int = 100,
) -> tuple[float, int]:
    """
    Regula falsi (false position): secant/linear-interpolation root estimate
    across the bracket — typically faster than bisection, but can stagnate if
    one endpoint sticks.  Auto-brackets like ``bisection``.

    Returns (root, iterations).
    """
    if f(a) * f(b) >= 0.0:
        a, b = bracket_root(f, a, b)
    for i in range(1, max_iterations + 1):
        fa, fb = f(a), f(b)
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        if abs(fc) < accuracy and abs(b - a) < accuracy:
            return c, i
        if fc * fa < 0.0:
            b = c
        else:
            a = c
    return c, max_iterations


def fixed_point(
    g: ScalarFunction, x0: float, accuracy: float = 1e-6, max_iterations: int = 30
) -> tuple[float, int]:
    """
    Picard fixed-point iteration x_{k+1} = g(x_k) for solving x = g(x).
    Converges when |g'(x)| < 1 in a neighbourhood of the fixed point.

    Returns (root, iterations).
    """
    x = float(x0)
    for i in range(1, max_iterations + 1):
        x_next = g(x)
        if abs(x_next - x) < accuracy:
            return x_next, i
        x = x_next
    return x, max_iterations


def finite_difference_derivative(
    f: ScalarFunction, x: float, h: float = 1e-5, order: int = 1
) -> float:
    """
    Central-difference derivatives (both O(h^2)):
      order=1:  f'(x)  ~ [f(x+h) - f(x-h)] / (2h)
      order=2:  f''(x) ~ [f(x+h) - 2f(x) + f(x-h)] / h^2
    """
    if order == 1:
        return (f(x + h) - f(x - h)) / (2.0 * h)
    if order == 2:
        return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)
    raise ValueError("order must be 1 or 2.")


def partial_derivative(
    f: Callable[..., float], var_index: int, point: Vector, h: float = 1e-5
) -> float:
    """
    Central-difference partial derivative of a multivariate scalar function
    with respect to ``point[var_index]``.
    """
    fwd, bwd = list(point), list(point)
    fwd[var_index] += h
    bwd[var_index] -= h
    return (f(*fwd) - f(*bwd)) / (2.0 * h)


def jacobian(f: Callable[..., list], point: Vector, h: float = 1e-6) -> list[list[float]]:
    """
    Numerical Jacobian of a vector-valued function f: R^n -> R^m at ``point``
    via central finite differences.  Returns an m x n matrix.
    """
    n = len(point)
    m = len(f(*point))
    J = [[0.0] * n for _ in range(m)]
    for i in range(m):
        def component(*args, _i=i) -> float:
            return f(*args)[_i]
        for j in range(n):
            J[i][j] = partial_derivative(component, j, point, h)
    return J


def newton_raphson(
    f: ScalarFunction, df: ScalarFunction | None = None, x0: float = 0.0,
    accuracy: float = 1e-6, max_iterations: int = 30, h: float = 1e-5,
) -> tuple[float, int]:
    """
    1-D Newton-Raphson:  x_{k+1} = x_k - f(x_k)/f'(x_k).  Quadratic
    convergence near a simple root (the error roughly squares each step).
    Uses the supplied analytical derivative ``df`` if given, otherwise an
    O(h^2) central difference.

    Returns (root, iterations).  Raises ValueError if the derivative vanishes.
    """
    x = float(x0)
    for i in range(max_iterations):
        fx = f(x)
        if abs(fx) < accuracy:
            return x, i
        dfx = df(x) if df is not None else (f(x + h) - f(x - h)) / (2.0 * h)
        if abs(dfx) < 1e-12:
            raise ValueError("Newton-Raphson failed: derivative is (near) zero.")
        x_next = x - fx / dfx
        if abs(x_next - x) < accuracy and abs(f(x_next)) < accuracy:
            return x_next, i + 1
        x = x_next
    return x, max_iterations


def newton_raphson_system(
    f: Callable[..., list],
    J: Callable[..., Matrix],
    initial_guess: Vector,
    accuracy: float = 1e-6,
    max_iterations: int = 30,
) -> tuple[list[float], int]:
    """
    Multivariate Newton-Raphson for a non-linear system F(x) = 0:
        x^{(k+1)} = x^{(k)} - J(x^{(k)})^{-1} F(x^{(k)}).
    The Newton step is obtained by *solving* J dx = F via Gauss-Jordan rather
    than forming J^{-1} explicitly (cheaper and better conditioned).
    Quadratic convergence near a simple root.

    Returns (x_solution, iterations).
    """
    x = [float(v) for v in initial_guess]
    n = len(x)
    for i in range(max_iterations):
        fx = [float(v) for v in f(*x)]
        if sum(abs(v) for v in fx) < accuracy:
            return x, i
        Jx = J(*x)
        # Solve Jx * dx = fx  by Gauss-Jordan on [Jx | fx].
        augmented = [Jx[k][:] + [fx[k]] for k in range(n)]
        dx = gauss_jordan_elimination_augmented(augmented)
        x = [x[k] - dx[k] for k in range(n)]
        if sum(abs(d) for d in dx) < accuracy:
            return x, i + 1
    return x, max_iterations


# =============================================================================
# 8. POLYNOMIAL ROOTS & DEFLATION: HORNER, LAGUERRE, SYNTHETIC DIVISION
# =============================================================================

def polynomial_value(coefficients: Sequence[float], x: float) -> float:
    """
    Evaluate P(x) = a_n x^n + ... + a_1 x + a_0 with Horner's scheme
    (coefficients in descending powers): O(n) operations, no powers of x.
    """
    result = float(coefficients[0])
    for c in coefficients[1:]:
        result = result * x + float(c)
    return result


def polynomial_first_derivative(coefficients: Sequence[float]) -> list[float]:
    """Coefficients of P'(x), descending powers:  d(a_k x^k) = k a_k x^{k-1}."""
    n = len(coefficients) - 1
    return [float(coefficients[i]) * (n - i) for i in range(n)]


def polynomial_second_derivative(coefficients: Sequence[float]) -> list[float]:
    """Coefficients of P''(x), descending powers."""
    return polynomial_first_derivative(polynomial_first_derivative(coefficients))


def laguerre(
    coefficients: Sequence[float], b0: float,
    accuracy1: float = 1e-8, accuracy2: float = 1e-6, max_iterations: int = 30,
) -> tuple[float, int]:
    """
    Laguerre's method for a single real root of P:  b_{k+1} = b_k - a, where
        a = n / (G ± sqrt((n-1)(nH - G^2))),   G = P'/P,  H = G^2 - P''/P.
    The sign in the denominator is chosen to maximise |denominator| (the
    numerically safe branch).  Near simple roots the convergence is cubic.

    Real roots only (negative radicands are clamped to zero).
    Returns (root, iterations).
    """
    n = len(coefficients) - 1
    first = polynomial_first_derivative(coefficients)
    second = polynomial_second_derivative(coefficients)
    b = float(b0)

    for i in range(1, max_iterations + 1):
        P = polynomial_value(coefficients, b)
        if abs(P) < accuracy1:
            return b, i
        if abs(P) < 1e-300:  # deeper than any meaningful root — call it done
            return b, i

        P1 = polynomial_value(first, b)
        P2 = polynomial_value(second, b)

        G = P1 / P
        H = G * G - P2 / P

        radicand = (n - 1) * (n * H - G * G)
        if radicand < 0.0:  # complex-root territory; stay on the real branch
            radicand = 0.0
        sqrt_term = math.sqrt(radicand)
        denom1, denom2 = G + sqrt_term, G - sqrt_term
        denom = denom1 if abs(denom1) >= abs(denom2) else denom2
        if abs(denom) < 1e-14:
            raise ValueError("Laguerre method: denominator vanished.")

        a = n / denom
        b_new = b - a
        if abs(b_new - b) < accuracy1 and abs(polynomial_value(coefficients, b_new)) < accuracy2:
            return b_new, i
        b = b_new

    raise ValueError(
        f"Laguerre method did not converge within {max_iterations} iterations."
    )


def synthetic_division(coefficients: Sequence[float], root: float,
                       accuracy: float = 1e-6) -> list[float]:
    """
    Deflate P by dividing by (x - root) via Horner's synthetic division.
    Returns the degree-(n-1) quotient; raises ValueError if the remainder is
    non-zero within ``accuracy`` (i.e. ``root`` is not a root).
    """
    quotient = [float(coefficients[0])]
    for c in coefficients[1:-1]:
        quotient.append(float(c) + root * quotient[-1])
    remainder = float(coefficients[-1]) + root * quotient[-1]
    if abs(remainder) > accuracy:
        raise ValueError(
            f"Synthetic division failed: remainder {remainder:.3e} "
            f"exceeds tolerance {accuracy:.1e}."
        )
    return quotient


def laguerre_roots(coefficients: Sequence[float], b0: float = 0.0):
    """
    Find all real roots of P: repeatedly apply Laguerre's method, deflate with
    synthetic division, and warm-start the next search from the root just
    found.  The final deflated linear factor a1 x + a0 gives the last root
    exactly.

    Returns (roots, per-root iteration counts).
    """
    roots: list[float] = []
    iterations: list[int] = []
    curr = [float(c) for c in coefficients]
    guess = float(b0)

    while len(curr) > 2:
        root, n_iter = laguerre(curr, guess)
        roots.append(root)
        iterations.append(n_iter)
        curr = synthetic_division(curr, root)
        guess = root  # warm start: roots tend to cluster in deflated polynomials

    final_root = -curr[1] / curr[0]
    roots.append(final_root)
    iterations.append(0)
    return roots, iterations


# =============================================================================
# 9. NUMERICAL INTEGRATION & QUADRATURE
# =============================================================================

def midpoint(f: ScalarFunction, a: float, b: float, N: int) -> float:
    """
    Composite midpoint (rectangle) rule:
        I ~ h * sum_i f(a + (i + 1/2) h),   h = (b-a)/N.
    Error:  |E| <= (b-a)^3 / (24 N^2) * max|f''|.
    """
    h = (b - a) / N
    return h * sum(f(a + (i + 0.5) * h) for i in range(N))


def trapezoidal(f: ScalarFunction, a: float, b: float, N: int) -> float:
    """
    Composite trapezoidal rule:
        I ~ (h/2) [ f(a) + 2 sum_i f(a + i h) + f(b) ],  h = (b-a)/N.
    Error:  |E| <= (b-a)^3 / (12 N^2) * max|f''|.
    """
    h = (b - a) / N
    total = f(a) + f(b) + 2.0 * sum(f(a + i * h) for i in range(1, N))
    return (h / 2.0) * total


def simpson(f: ScalarFunction, a: float, b: float, N: int) -> tuple[float, int]:
    """
    Composite Simpson's 1/3 rule (N must be even):
        I ~ (h/3) [ f(a) + 4 sum_{odd} f(x_i) + 2 sum_{even} f(x_i) + f(b) ].
    Error:  |E| <= (b-a)^5 / (180 N^4) * max|f''''|.
    Returns (integral, number of function evaluations).
    """
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    if N % 2 != 0:
        raise ValueError("Simpson's 1/3 rule requires an even number of intervals N.")
    h = (b - a) / N
    total = f(a) + f(b)
    evaluations = 2
    for i in range(1, N):
        x = a + i * h
        total += (2.0 if i % 2 == 0 else 4.0) * f(x)
        evaluations += 1
    return (h / 3.0) * total, evaluations


def monte_carlo(
    f: ScalarFunction, a: float, b: float, N: int, seed: int = 1
) -> tuple[float, float, float]:
    """
    Uniform Monte Carlo integration:  I ~ (b-a) * mean(f(x_i)),  x_i uniform
    on [a, b], drawn from the in-house LCG (seeded for reproducibility).
    Statistical error scales as 1/sqrt(N), independent of dimension.
    Returns (integral, sigma_f, sigma_integral).
    """
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    myrand(seed)
    sum_f = 0.0
    sum_f2 = 0.0
    for _ in range(N):
        val = f(a + (b - a) * myrand())
        sum_f += val
        sum_f2 += val * val
    mean_f = sum_f / N
    variance_f = max(0.0, sum_f2 / N - mean_f * mean_f)  # guard round-off negatives
    sigma_f = math.sqrt(variance_f)
    sigma_integral = (b - a) * sigma_f / math.sqrt(N)
    return (b - a) * mean_f, sigma_f, sigma_integral


def integration_error_bound_N(
    method: str, a: float, b: float, max_deriv: float, target_error: float = 1e-6
) -> int:
    """
    Minimum number of subintervals N from the classical error bounds:
      Midpoint:    N = ceil( sqrt( (b-a)^3 M2 / (24 eps) ) )
      Trapezoidal: N = ceil( sqrt( (b-a)^3 M2 / (12 eps) ) )
      Simpson:     N = ceil( ( (b-a)^5 M4 / (180 eps) )^(1/4) ), rounded up to even
    """
    L = abs(b - a)
    m = method.lower()
    if m in ("midpoint", "mid"):
        return math.ceil(math.sqrt(L**3 * max_deriv / (24.0 * target_error)))
    if m in ("trapezoidal", "trap", "trapezoid"):
        return math.ceil(math.sqrt(L**3 * max_deriv / (12.0 * target_error)))
    if m in ("simpson", "simp"):
        n = math.ceil((L**5 * max_deriv / (180.0 * target_error)) ** 0.25)
        return n if n % 2 == 0 else n + 1
    raise ValueError("method must be 'midpoint', 'trapezoidal', or 'simpson'.")


# Quadrature nodes/weights.  Legendre: interval [-1, 1].  Laguerre: interval
# [0, inf) with weight function w(x) = e^{-x}.  Values are the exact
# orthogonal-polynomial zeros/weights, listed to full double precision.
GAUSS_LEGENDRE: dict[int, tuple[list[float], list[float]]] = {
    1: ([0.0], [2.0]),
    2: ([-0.5773502691896258, 0.5773502691896258], [1.0, 1.0]),
    3: ([-0.7745966692414834, 0.0, 0.7745966692414834],
        [0.5555555555555556, 0.8888888888888888, 0.5555555555555556]),
    4: ([-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526],
        [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]),
    5: ([-0.9061798459386640, -0.5384693101056831, 0.0,
         0.5384693101056831, 0.9061798459386640],
        [0.2369268850561891, 0.4786286704993665, 0.5688888888888888,
         0.4786286704993665, 0.2369268850561891]),
    6: ([-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969,
         0.6612093864662645, 0.9324695142031521],
        [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910,
         0.3607615730481386, 0.1713244923791704]),
}

GAUSS_LAGUERRE: dict[int, tuple[list[float], list[float]]] = {
    1: ([1.0], [1.0]),
    2: ([0.5857864376269050, 3.4142135623730950],
        [0.8535533905932737, 0.1464466094067262]),
    3: ([0.4157745567834778, 2.2942803602790400, 6.2899450829374794],
        [0.7110930099291730, 0.2785177335692410, 0.0103892565015861]),
    4: ([0.3225476896193924, 1.7457611011583471, 4.5366202969211290, 9.3950709123011311],
        [0.6031541043416337, 0.3574186924377996, 0.0388879085150054, 0.0005392947055613]),
    5: ([0.2635603197181409, 1.4134030591065173, 3.5964257710407220, 7.0858100058588367, 12.6408008442757825],
        [0.5217556105828085, 0.3986668110831760, 0.0759424496817077, 0.0036117586799221, 0.0000233699723858]),
    6: ([0.2228466041792608, 1.1889321016726230, 2.9927363260593143, 5.7751435691045102, 9.8374674183825886, 15.9828739806017821],
        [0.4589646739499648, 0.4170008307721200, 0.1133733820740449, 0.0103991974531491, 0.0002610172028149, 0.0000008985479064]),
}


def gaussian_quadrature(
    f: ScalarFunction, a: float, b: float, N: int, method: str = "legendre"
) -> float:
    """
    Gaussian quadrature with a single routine for both families:
      * 'legendre' — integrates f over [a, b]; the standard nodes on [-1, 1]
        are mapped to [a, b] by  x = ((b-a)/2) t + (b+a)/2.
      * 'laguerre' — integrates e^{-t} f(t) over [0, inf); ``a``, ``b`` are
        ignored by convention (pass 0 and a large number).
    N must be one of {1..6} (the tabulated orders).
    Exact for polynomials of degree <= 2N - 1 (up to round-off).
    """
    m = method.lower()
    if m == "legendre":
        if N not in GAUSS_LEGENDRE:
            raise ValueError(f"Gauss-Legendre: N must be in {sorted(GAUSS_LEGENDRE)}")
        nodes, weights = GAUSS_LEGENDRE[N]
        scale = (b - a) / 2.0
        shift = (b + a) / 2.0
        return scale * sum(w * f(scale * t + shift) for t, w in zip(nodes, weights))
    if m == "laguerre":
        if N not in GAUSS_LAGUERRE:
            raise ValueError(f"Gauss-Laguerre: N must be in {sorted(GAUSS_LAGUERRE)}")
        nodes, weights = GAUSS_LAGUERRE[N]
        return sum(w * f(t) for t, w in zip(nodes, weights))
    raise ValueError("method must be 'legendre' or 'laguerre'.")
