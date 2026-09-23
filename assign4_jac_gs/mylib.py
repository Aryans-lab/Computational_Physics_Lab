# Computational Physics Laboratory (PHY341 / PHY745)
# Core Library - mylib.py
# Author: Aryan Bandyopadhyay | Roll Number: 2411014
# School of Physical Sciences, NISER Bhubaneswar

import os
import math

# Default base directory for lab data files
base_dir = "D:\\NISER\\Computational Physics Lab"


# =============================================================================
# 1. RANDOM NUMBER GENERATION & DISTRIBUTIONS
# =============================================================================

def myrand(seed=0):
    """
    Linear Congruential Generator (LCG) for pseudo-random numbers in [0, 1).
    Recurrence: x_{i+1} = (a * x_i + c) mod m
    Parameters follow GCC runtime standard:
        a = 1103515245 (multiplier)
        c = 12345      (increment)
        m = 32768      (modulus, 2^15)
    Seed is remembered across function calls via function attribute.
    """
    if not hasattr(myrand, 'x0'):
        myrand.x0 = 0
    a = 1103515245
    c = 12345
    m = 32768
    if seed:
        myrand.x0 = seed
    x1 = (a * myrand.x0 + c) % m
    rn = x1 / m
    myrand.x0 = x1
    return rn


def random_uniform(a=0.0, b=1.0):
    """
    Generate a pseudo-random floating point number uniformly distributed in [a, b).
    Linear scaling: X = a + (b - a) * xi, where xi in [0, 1).
    """
    return a + (b - a) * myrand()


def random_exponential(lam=1.0):
    """
    Generate an exponentially distributed random deviate with PDF q(y) = lambda * exp(-lambda * y).
    Uses the inverse-transform method: y = - (1 / lambda) * ln(u), where u in (0, 1].
    """
    u = myrand()
    while u <= 0.0:
        u = myrand()
    return - (1.0 / lam) * math.log(u)


# =============================================================================
# 2. COMPLEX NUMBER ALGEBRA
# =============================================================================

class MyComplex:
    """
    Self-contained complex number class implementing arithmetic and modulus.
    Designed without relying on Python's built-in complex type.
    """
    def __init__(self, real, image=0.0):
        self.real = float(real)
        self.image = float(image)

    def display_complex(self):
        sign = "+" if self.image >= 0 else "-"
        print(f"{self.real} {sign} {abs(self.image)}j")

    def add_complex(self, c1, c2):
        return MyComplex(c1.real + c2.real, c1.image + c2.image)

    def subtract_complex(self, c1, c2):
        return MyComplex(c1.real - c2.real, c1.image - c2.image)

    def multiply_complex(self, c1, c2):
        real_part = c1.real * c2.real - c1.image * c2.image
        image_part = c1.real * c2.image + c1.image * c2.real
        return MyComplex(real_part, image_part)

    def modulus_complex(self, c):
        return (c.real * c.real + c.image * c.image) ** 0.5

    # Operator overloading helpers for clean expressions in code
    def __add__(self, other):
        if isinstance(other, MyComplex):
            return MyComplex(self.real + other.real, self.image + other.image)
        return MyComplex(self.real + other, self.image)

    def __sub__(self, other):
        if isinstance(other, MyComplex):
            return MyComplex(self.real - other.real, self.image - other.image)
        return MyComplex(self.real - other, self.image)

    def __mul__(self, other):
        if isinstance(other, MyComplex):
            r = self.real * other.real - self.image * other.image
            im = self.real * other.image + self.image * other.real
            return MyComplex(r, im)
        return MyComplex(self.real * other, self.image * other)

    def __abs__(self):
        return (self.real**2 + self.image**2)**0.5

    def __repr__(self):
        sign = "+" if self.image >= 0 else "-"
        return f"MyComplex({self.real}, {self.image})"

    def __str__(self):
        sign = "+" if self.image >= 0 else "-"
        return f"{self.real} {sign} {abs(self.image)}j"


# =============================================================================
# 3. MATRIX & VECTOR UTILITIES / FILE I/O
# =============================================================================

def read_matrix_from_file(filename):
    """
    Reads an ASCII text file containing rows of space/tab-separated numbers into a 2D list.
    Checks current working directory, full path, or base_dir fallback.
    """
    if os.path.exists(filename):
        filepath = filename
    elif os.path.isabs(filename):
        filepath = filename
    else:
        filepath = os.path.join(base_dir, filename)

    matrix = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                row = [float(x) for x in line.split()]
                matrix.append(row)
    return matrix


def write_matrix_to_file(matrix, filename, precision=6):
    """
    Writes a 2D list or 1D list to an external text file in aligned column format.
    Ensures non-interactive output complies with lab submission requirements.
    """
    with open(filename, 'w') as f:
        # Check if 1D list or 2D list
        if not matrix:
            return
        if not isinstance(matrix[0], (list, tuple)):
            for val in matrix:
                f.write(f"{val:.{precision}f}\n")
        else:
            for row in matrix:
                row_str = "  ".join(f"{val:12.{precision}f}" for val in row)
                f.write(row_str + "\n")


def print_matrix(matrix, label=None, precision=4):
    """
    Neatly prints a matrix or vector to terminal with formatted column alignment.
    """
    if label:
        print(f"--- {label} ---")
    if not matrix:
        print("[]")
        return
    if not isinstance(matrix[0], (list, tuple)):
        # 1D vector
        print("[" + ", ".join(f"{x:.{precision}f}" for x in matrix) + "]")
    else:
        for row in matrix:
            print("  [" + "  ".join(f"{x:10.{precision}f}" for x in row) + " ]")


def matrix_multiply(X, Y):
    """
    Multiplies matrix X (shape n x m) by matrix Y (shape m x p).
    Returns product matrix of shape n x p.
    """
    rows_X = len(X)
    col_X = len(X[0])
    rows_Y = len(Y)
    col_Y = len(Y[0])

    if col_X != rows_Y:
        raise ValueError(f"Matrix dimension mismatch: ({rows_X}x{col_X}) and ({rows_Y}x{col_Y}) cannot be multiplied.")

    result = [[0.0 for _ in range(col_Y)] for _ in range(rows_X)]
    for i in range(rows_X):
        for j in range(col_Y):
            for k in range(col_X):
                result[i][j] += X[i][k] * Y[k][j]
    return result


def dot_product_vector(X, Y):
    """
    Computes dot product of two vectors X and Y.
    Supports 1D flat lists [x1, x2, ...] or column vectors [[x1], [x2], ...].
    """
    if len(X) != len(Y):
        print("Vectors are not of the same length")
        return None

    result = 0.0
    for i in range(len(X)):
        val_x = X[i][0] if isinstance(X[i], (list, tuple)) else X[i]
        val_y = Y[i][0] if isinstance(Y[i], (list, tuple)) else Y[i]
        result += val_x * val_y
    return result


def matrix_transpose(A):
    """
    Returns the transpose A^T of a matrix A.
    """
    rows = len(A)
    cols = len(A[0])
    return [[A[i][j] for i in range(rows)] for j in range(cols)]


def matrix_add(A, B):
    """
    Element-wise addition of two matrices of identical dimension.
    """
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def matrix_sub(A, B):
    """
    Element-wise subtraction (A - B) of two matrices of identical dimension.
    """
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def matrix_vector_multiply(A, x):
    """
    Matrix-vector multiplication A * x.
    Accepts x as a 1D list [x0, x1, ...] or 2D column vector [[x0], [x1], ...].
    Returns a 1D list.
    """
    n = len(A)
    m = len(A[0])
    xv = [x[i][0] if isinstance(x[i], (list, tuple)) else x[i] for i in range(len(x))]
    if len(xv) != m:
        raise ValueError(f"Dimension mismatch in matrix-vector multiply: ({n}x{m}) with vector length {len(xv)}")
    result = [0.0] * n
    for i in range(n):
        result[i] = sum(A[i][j] * xv[j] for j in range(m))
    return result


def vector_norm(v, p=2):
    """
    Computes the p-norm of vector v:
        p=2:        Euclidean norm sqrt(sum(v_i^2))
        p=1:        Manhattan norm sum(|v_i|)
        p='inf':    Infinity / Chebyshev norm max(|v_i|)
    """
    flat = [v[i][0] if isinstance(v[i], (list, tuple)) else v[i] for i in range(len(v))]
    if p == 2:
        return math.sqrt(sum(x * x for x in flat))
    elif p == 1:
        return sum(abs(x) for x in flat)
    elif p == 'inf' or p == float('inf'):
        return max(abs(x) for x in flat)
    else:
        return (sum(abs(x)**p for x in flat)) ** (1.0 / p)


def matrix_residual(A, x, b):
    """
    Computes the residual norm ||A*x - b||_2 to verify accuracy of linear system solutions.
    An exact solution yields residual on the order of machine precision (~ 1e-15).
    """
    Ax = matrix_vector_multiply(A, x)
    bv = [b[i][0] if isinstance(b[i], (list, tuple)) else b[i] for i in range(len(b))]
    r = [Ax[i] - bv[i] for i in range(len(bv))]
    return vector_norm(r, p=2)


def is_symmetric(A, tol=1e-9):
    """
    Checks whether square matrix A is symmetric (A == A^T) within numerical tolerance.
    """
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j] - A[j][i]) > tol:
                return False
    return True


def is_diagonally_dominant(A, strict=False):
    """
    Checks if square matrix A is diagonally dominant:
        |A[i][i]| >= sum_{j != i} |A[i][j]|  (weak)
        |A[i][i]| >  sum_{j != i} |A[i][j]|  (strict, if strict=True)
    Crucial check for guaranteed convergence of Jacobi and Gauss-Seidel iterations.
    """
    n = len(A)
    for i in range(n):
        off_diag_sum = sum(abs(A[i][j]) for j in range(n) if j != i)
        if strict:
            if abs(A[i][i]) <= off_diag_sum:
                return False
        else:
            if abs(A[i][i]) < off_diag_sum:
                return False
    return True


# =============================================================================
# 4. DIRECT LINEAR SYSTEM SOLVERS: GAUSS-JORDAN & LU DECOMPOSITION
# =============================================================================

def gauss_jordan_elimination_augmented(augmented_matrix):
    """
    Solves linear system A*x = b using Gauss-Jordan elimination with partial pivoting.
    Input: augmented_matrix of shape n x (n + 1) -> [A | b].
    Transforms system to Reduced Row Echelon Form (RREF).
    Returns solution vector x of length n.
    Raises ValueError for singular, inconsistent, or dependent systems.
    """
    if not augmented_matrix or not augmented_matrix[0]:
        return []

    augmented = [[float(x) for x in row] for row in augmented_matrix]
    n = len(augmented)
    m = len(augmented[0])

    if m != n + 1:
        raise ValueError("Augmented matrix must have shape n x (n+1).")

    for i in range(n):
        # Partial pivoting: select pivot row with largest magnitude in column i
        pivot_row = max(range(i, n), key=lambda r: abs(augmented[r][i]))

        if abs(augmented[pivot_row][i]) < 1e-14:
            continue

        # Swap rows to bring largest pivot to current diagonal
        if pivot_row != i:
            augmented[i], augmented[pivot_row] = augmented[pivot_row], augmented[i]

        # Normalize pivot row
        pivot = augmented[i][i]
        for j in range(i, n + 1):
            augmented[i][j] /= pivot

        # Eliminate column i in all other rows
        for k in range(n):
            if k != i:
                factor = augmented[k][i]
                for j in range(i, n + 1):
                    augmented[k][j] -= factor * augmented[i][j]

    # Check for consistency (no-solution or infinitely many solutions)
    for i in range(n):
        all_zero = all(abs(augmented[i][j]) < 1e-12 for j in range(n))
        if all_zero and abs(augmented[i][n]) > 1e-12:
            raise ValueError("Incompatible system: no solution exists.")
        if all_zero and abs(augmented[i][n]) <= 1e-12:
            raise ValueError("Dependent system: infinitely many solutions.")

    return [augmented[i][n] for i in range(n)]


def gauss_jordan_inverse(A):
    """
    Computes the inverse of an n x n square matrix A using Gauss-Jordan elimination
    on the augmented matrix [A | I_n].
    Returns inverse matrix A^{-1} as a 2D list.
    """
    if not A or not A[0]:
        return []

    A_mat = [[float(x) for x in row] for row in A]
    n = len(A_mat)
    m = len(A_mat[0])

    if m != n:
        raise ValueError("Gauss-Jordan inverse requires a square matrix.")

    # Augmented matrix [A | I]
    augmented = [A_mat[i][:] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for i in range(n):
        # Partial pivoting
        pivot_row = max(range(i, n), key=lambda r: abs(augmented[r][i]))
        if abs(augmented[pivot_row][i]) < 1e-14:
            raise ValueError("Matrix is singular; inverse does not exist.")

        if pivot_row != i:
            augmented[i], augmented[pivot_row] = augmented[pivot_row], augmented[i]

        pivot = augmented[i][i]
        for j in range(i, 2 * n):
            augmented[i][j] /= pivot

        for k in range(n):
            if k != i:
                factor = augmented[k][i]
                for j in range(i, 2 * n):
                    augmented[k][j] -= factor * augmented[i][j]

    inverse = [row[n:] for row in augmented]
    return inverse


def lu_decomposition(A):
    """
    Doolittle's LU Decomposition of square matrix A:
        A = L * U
    where L is unit lower triangular (L[i][i] = 1) and U is upper triangular.
    Returns: [L, U] as lists of 2D lists.
    """
    if not A or not A[0]:
        return []

    A_mat = [[float(x) for x in row] for row in A]
    n = len(A_mat)
    if any(len(row) != n for row in A_mat):
        raise ValueError("Matrix A is not a square matrix.")

    # Row swap if A[0][0] == 0 to avoid zero pivot in first step
    if abs(A_mat[0][0]) < 1e-14:
        for i in range(1, n):
            if abs(A_mat[i][0]) > 1e-14:
                A_mat[0], A_mat[i] = A_mat[i], A_mat[0]
                break

    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[0.0 for _ in range(n)] for _ in range(n)]

    # Doolittle algorithm: compute row i of U, then column i of L
    for i in range(n):
        for j in range(i, n):
            U[i][j] = A_mat[i][j] - sum(L[i][k] * U[k][j] for k in range(i))

        if abs(U[i][i]) < 1e-14:
            raise ValueError(f"LU decomposition failed: zero pivot encountered at diagonal index {i}.")

        for j in range(i + 1, n):
            L[j][i] = (A_mat[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return [L, U]


def lu_forback(A, b, method='lu'):
    """
    Solves A * x = b using LU decomposition (Doolittle) or Cholesky:
        1. L * y = b  (Forward substitution)
        2. U * x = y  (Backward substitution)
    Accepts filenames (strings) or nested lists for A and b.
    Returns solution vector x as a list of floats.
    """
    # Load A if filename provided
    if isinstance(A, str):
        C = read_matrix_from_file(A)
    else:
        C = [[float(x) for x in row] for row in A]

    # Load b if filename provided
    if isinstance(b, str):
        d_mat = read_matrix_from_file(b)
        bv = [d_mat[i][0] for i in range(len(d_mat))]
    elif isinstance(b[0], (list, tuple)):
        bv = [float(row[0]) for row in b]
    else:
        bv = [float(x) for x in b]

    n = len(C)

    if method == 'Cholesky':
        L = cholesky_decomposition(C)
        return cholesky_forback(L, bv)

    # Doolittle LU
    L, U = lu_decomposition(C)

    # Forward substitution: L * y = b  (L has 1 on diagonal)
    y = []
    for i in range(n):
        sum_Ly = sum(L[i][j] * y[j] for j in range(i))
        y.append(bv[i] - sum_Ly)

    # Backward substitution: U * x = y
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum_Ux = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - sum_Ux) / U[i][i]

    return x


def matrix_determinant_lu(A):
    """
    Calculates determinant of square matrix A from Doolittle LU decomposition:
        det(A) = det(L) * det(U) = 1.0 * prod(U[i][i])
    Returns float determinant value.
    """
    L, U = lu_decomposition(A)
    det = 1.0
    for i in range(len(U)):
        det *= U[i][i]
    return det


# =============================================================================
# 5. SYMMETRIC POSITIVE-DEFINITE SYSTEMS: CHOLESKY DECOMPOSITION
# =============================================================================

def cholesky_decomposition(A):
    """
    Factorizes a real symmetric positive-definite matrix A into A = L * L^T,
    where L is a lower-triangular matrix with positive diagonal entries.
    Includes explicit verification of symmetry and positive-definiteness.
    """
    if not A:
        return []

    A_mat = [[float(x) for x in row] for row in A]
    n = len(A_mat)

    if any(len(row) != n for row in A_mat):
        raise ValueError("Matrix A is not a square matrix.")

    # Symmetry check: A[i][j] == A[j][i]
    for i in range(n):
        for j in range(i):
            if abs(A_mat[i][j] - A_mat[j][i]) > 1e-8:
                raise ValueError("Matrix A is not symmetric; Cholesky requires a symmetric matrix.")

    L = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            val = A_mat[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if val <= 0.0:
                    raise ValueError(f"Matrix A is not positive-definite (leading principal minor <= 0 at index {i}).")
                L[i][j] = math.sqrt(val)
            else:
                L[i][j] = val / L[j][j]

    return L


def cholesky_forback(L, b):
    """
    Solves A * x = b given Cholesky factor L (such that A = L * L^T):
        1. L * y = b     (Forward substitution)
        2. L^T * x = y   (Backward substitution)
    Supports column vector [[b0], [b1], ...] or 1D list [b0, b1, ...].
    """
    n = len(L)
    bv = [b[i][0] if isinstance(b[i], (list, tuple)) else float(b[i]) for i in range(n)]

    # Forward substitution: L * y = b
    y = []
    for i in range(n):
        sum1 = sum(L[i][j] * y[j] for j in range(i))
        y.append((bv[i] - sum1) / L[i][i])

    # Backward substitution: L^T * x = y
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum2 = sum(L[j][i] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - sum2) / L[i][i]

    return x


# =============================================================================
# 6. ITERATIVE LINEAR SOLVERS: JACOBI, GAUSS-SEIDEL, AND SOR
# =============================================================================

def jacobi_it(A, b, tol=1e-8, max_iterations=1000):
    """
    Solves linear system A * x = b using the Jacobi iterative method:
        x_i^{(k+1)} = (1 / A_{ii}) * [ b_i - sum_{j != i} A_{ij} x_j^{(k)} ]
    Performs diagonal dominance row pivoting if needed.
    """
    if not A or not A[0] or not b:
        return []

    A_mat = [[float(x) for x in row] for row in A]
    bv = [float(row[0]) if isinstance(row, (list, tuple)) else float(row) for row in b]
    n = len(bv)
    x = [0.0] * n

    # Row swapping to achieve diagonal dominance if possible
    for i in range(n):
        off_sum = sum(abs(A_mat[i][j]) for j in range(n) if j != i)
        if abs(A_mat[i][i]) < off_sum:
            for j in range(i + 1, n):
                if abs(A_mat[j][i]) > abs(A_mat[i][i]):
                    A_mat[i], A_mat[j] = A_mat[j], A_mat[i]
                    bv[i], bv[j] = bv[j], bv[i]
                    break

    # Jacobi iterations (all updates computed from previous vector x_prev)
    for k in range(1, max_iterations + 1):
        x_prev = x.copy()
        for i in range(n):
            sum_other = sum(A_mat[i][j] * x_prev[j] for j in range(n) if j != i)
            x[i] = (bv[i] - sum_other) / A_mat[i][i]

        if all(abs(x[i] - x_prev[i]) < tol for i in range(n)):
            print(f"Converged after {k} iterations")
            return x

    raise ValueError(f"Jacobi method couldn't converge within {max_iterations} iterations.")


def gauss_seidel(A, b, tol=1e-8, max_iterations=500):
    """
    Solves linear system A * x = b using Gauss-Seidel iterations:
        x_i^{(k+1)} = (1 / A_{ii}) * [ b_i - sum_{j < i} A_{ij} x_j^{(k+1)} - sum_{j > i} A_{ij} x_j^{(k)} ]
    Uses freshly updated values in-place immediately, accelerating convergence.
    """
    if not A or not A[0] or not b:
        return []

    A_mat = [[float(x) for x in row] for row in A]
    bv = [float(row[0]) if isinstance(row, (list, tuple)) else float(row) for row in b]
    n = len(bv)
    x = [0.0] * n

    # Row swapping to enforce diagonal dominance
    for i in range(n):
        off_sum = sum(abs(A_mat[i][j]) for j in range(n) if j != i)
        if abs(A_mat[i][i]) < off_sum:
            for j in range(i + 1, n):
                if abs(A_mat[j][i]) > abs(A_mat[i][i]):
                    A_mat[i], A_mat[j] = A_mat[j], A_mat[i]
                    bv[i], bv[j] = bv[j], bv[i]
                    break

    # Gauss-Seidel iterations
    for k in range(1, max_iterations + 1):
        x_prev = x.copy()
        for i in range(n):
            sum_lower = sum(A_mat[i][j] * x[j] for j in range(i))
            sum_upper = sum(A_mat[i][j] * x[j] for j in range(i + 1, n))
            x[i] = (bv[i] - sum_lower - sum_upper) / A_mat[i][i]

        if all(abs(x[i] - x_prev[i]) < tol for i in range(n)):
            print(f"Converged after {k} iterations")
            return x

    raise ValueError(f"Gauss-Seidel method couldn't converge within {max_iterations} iterations.")


def sor_gauss_seidel(A, b, tol=1e-8, max_iterations=500, omega=1.0):
    """
    Successive Over-Relaxation (SOR) method:
        x_i^{(k+1)} = (1 - omega) * x_i^{(k)} + (omega / A_{ii}) * [ b_i - sum_{j < i} A_{ij} x_j^{(k+1)} - sum_{j > i} A_{ij} x_j^{(k)} ]
    omega = 1: Standard Gauss-Seidel
    1 < omega < 2: Over-relaxation (accelerates convergence for slow systems)
    0 < omega < 1: Under-relaxation (stabilizes oscillating or diverging systems)
    """
    if not A or not A[0] or not b:
        return []

    A_mat = [[float(x) for row_x in row for x in [row_x]] if isinstance(row, list) else [float(row)] for row in A]
    bv = [float(row[0]) if isinstance(row, (list, tuple)) else float(row) for row in b]
    n = len(bv)
    x = [0.0] * n

    # Row swapping to improve diagonal dominance
    for i in range(n):
        off_sum = sum(abs(A_mat[i][j]) for j in range(n) if j != i)
        if abs(A_mat[i][i]) < off_sum:
            for j in range(i + 1, n):
                if abs(A_mat[j][i]) > abs(A_mat[i][i]):
                    A_mat[i], A_mat[j] = A_mat[j], A_mat[i]
                    bv[i], bv[j] = bv[j], bv[i]
                    break

    # SOR iteration loop
    for k in range(1, max_iterations + 1):
        x_prev = x.copy()
        for i in range(n):
            sum_other = sum(A_mat[i][j] * x[j] for j in range(n) if j != i)
            x_gs = (bv[i] - sum_other) / A_mat[i][i]
            x[i] = omega * x_gs + (1.0 - omega) * x_prev[i]

        if all(abs(x[i] - x_prev[i]) < tol for i in range(n)):
            print(f"Converged after {k} iterations")
            return x

    raise ValueError(f"SOR method couldn't converge within {max_iterations} iterations.")


# =============================================================================
# 7. ROOT FINDING FOR NON-LINEAR EQUATIONS
# =============================================================================

def bracket_root(f, a, b, beta=0.1, max_iter=100):
    """
    Expands interval [a, b] until f(a) * f(b) < 0 (intermediate value theorem bracket).
    Follows lecture slide algorithm:
        If |f(a)| < |f(b)|: a = a - beta * (b - a)
        If |f(b)| < |f(a)|: b = b + beta * (b - a)
    """
    for _ in range(max_iter):
        fa = f(a)
        fb = f(b)
        if fa * fb < 0.0:
            return a, b
        if abs(fa) < abs(fb):
            a -= beta * (b - a)
        else:
            b += beta * (b - a)
            beta += 0.05
    raise ValueError(f"Failed to bracket a root in {max_iter} expansions.")


def bisection(f, a, b, accuracy=1e-8, max_iterations=100):
    """
    Bisection method for finding a root of f(x) = 0 in [a, b].
    Shrinks interval by half each step: c = (a + b) / 2.
    Guaranteed linear convergence if f(a) * f(b) < 0.
    """
    # Auto-bracket if boundary signs match
    if f(a) * f(b) >= 0:
        m = 0.5
        while f(a) * f(b) >= 0:
            if abs(f(a)) < abs(f(b)):
                a = a - m * (b - a)
            else:
                b = b + m * (b - a)
                m += 0.1

    for i in range(max_iterations):
        c = (a + b) / 2.0
        fc = f(c)
        if abs(fc) < accuracy and abs(b - a) < accuracy:
            print(f"The root is approximately {c}")
            print(f"The number of iterations taken is {i + 1}")
            return c
        if fc * f(a) < 0:
            b = c
        else:
            a = c

    print(f"The root is approximately {c}")
    print(f"The number of iterations taken is {max_iterations}")
    return c


def regula_falsi(f, a, b, accuracy=1e-8, max_iterations=100):
    """
    Regula Falsi (False Position) method:
    Uses linear interpolation across bracket [a, b]:
        c = b - f(b) * (b - a) / (f(b) - f(a))
    Faster than bisection for smooth functions.
    """
    if f(a) * f(b) >= 0:
        m = 0.5
        while f(a) * f(b) >= 0:
            if abs(f(a)) < abs(f(b)):
                a = a - m * (b - a)
            else:
                b = b + m * (b - a)
                m += 0.2

    for i in range(max_iterations):
        fa = f(a)
        fb = f(b)
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)

        if abs(fc) < accuracy and abs(b - a) < accuracy:
            print(f"The root is approximately {c}")
            print(f"The number of iterations taken is {i + 1}")
            return c

        if fc * fa < 0:
            b = c
        else:
            a = c

    print(f"The root is approximately {c}")
    print(f"The number of iterations taken is {max_iterations}")
    return c


def fixed_point(g, x0, accuracy=1e-6, max_iterations=30):
    """
    Picard Fixed-Point Iteration: solves x = g(x).
    Convergence requires |g'(x)| < 1 near the root.
    """
    x = x0
    for i in range(max_iterations):
        x_next = g(x)
        if abs(x_next - x) < accuracy:
            print(f"The root is: {x_next}")
            print(f"The number of iterations taken is {i + 1}")
            return x_next
        x = x_next

    print("Fixed-point method did not converge within the maximum iterations.")
    print(f"The last computed value is: {x}")
    return x, max_iterations


def finite_difference_derivative(f, x, h=1e-5, order=1):
    """
    Numerical derivative via central differences:
        order=1: f'(x)  ~ [f(x + h) - f(x - h)] / (2h)          O(h^2)
        order=2: f''(x) ~ [f(x + h) - 2f(x) + f(x - h)] / (h^2) O(h^2)
    """
    if order == 1:
        return (f(x + h) - f(x - h)) / (2.0 * h)
    elif order == 2:
        return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)
    else:
        raise ValueError("order must be 1 or 2.")


def partial_derivative(f, var_index, point, h=1e-5):
    """
    Computes partial derivative of multivariate scalar function f with respect to x_{var_index}
    using the central difference scheme:
        df/dx_i ~ [f(..., x_i + h, ...) - f(..., x_i - h, ...)] / (2h)
    """
    point_forward = list(point)
    point_backward = list(point)
    point_forward[var_index] += h
    point_backward[var_index] -= h
    return (f(*point_forward) - f(*point_backward)) / (2.0 * h)


def jacobian(f, point, h=1e-6):
    """
    Computes numerical Jacobian matrix J of a vector-valued function f: R^n -> R^m
    at a given point using central finite differences:
        J_{ij} = d(f_i) / d(x_j)
    Returns J as a 2D list of shape m x n.
    """
    n = len(point)
    f0 = f(*point)
    m = len(f0)
    J = [[0.0 for _ in range(n)] for _ in range(m)]
    for i in range(m):
        for j in range(n):
            def fi(*args):
                return f(*args)[i]
            J[i][j] = partial_derivative(fi, j, point, h)
    return J


def newton_raphson(f, df=None, x0=0.0, accuracy=1e-6, max_iterations=30, h=1e-5):
    """
    Newton-Raphson method for 1D root finding:
        x_{n+1} = x_n - f(x_n) / f'(x_n)
    If analytical derivative df is not supplied, uses central difference derivative.
    Exhibits quadratic convergence near simple roots.
    """
    x = float(x0)
    for i in range(max_iterations):
        fx = f(x)
        if abs(fx) < accuracy:
            print(f"The root is: {x}")
            print(f"The number of iterations taken is {i}")
            return x

        if df is not None:
            dfx = df(x)
        else:
            dfx = (f(x + h) - f(x - h)) / (2.0 * h)

        if abs(dfx) < 1e-12:
            print("Error: Derivative too close to zero; Newton-Raphson fails.")
            return None

        x_next = x - (fx / dfx)

        if abs(x_next - x) < accuracy and abs(f(x_next)) < accuracy:
            print(f"The root is: {x_next}")
            print(f"The number of iterations taken is {i + 1}")
            return x_next

        x = x_next

    print("Newton-Raphson did not converge within maximum iterations.")
    return x, max_iterations


def newton_raphson_system(f, J, initial_guess, accuracy=1e-6, max_iterations=30):
    """
    Multivariate Newton-Raphson method for non-linear systems f(x) = 0:
        x^{(k+1)} = x^{(k)} - [J(x^{(k)})]^{-1} * f(x^{(k)})
    Inversion of Jacobian matrix J is performed via Gauss-Jordan elimination.
    Returns: (x_solution, iterations_taken).
    """
    x = [float(v) for v in initial_guess]
    for i in range(max_iterations):
        fx = f(*x)
        norm_fx = sum(abs(val) for val in fx)
        if norm_fx < accuracy:
            return x, i

        Jx = J(*x)
        J_inv = gauss_jordan_inverse(Jx)
        delta_x = [sum(J_inv[j][k] * fx[k] for k in range(len(fx))) for j in range(len(x))]
        x_next = [x[j] - delta_x[j] for j in range(len(x))]

        # Convergence on step size
        if sum(abs(delta_x[j]) for j in range(len(x))) < accuracy:
            return x_next, i + 1

        x = x_next

    return x, max_iterations


# =============================================================================
# 8. POLYNOMIAL ROOTS & DEFLATION: HORNER & LAGUERRE METHODS
# =============================================================================

def polynomial_value(coefficients, x):
    """
    Evaluates polynomial P(x) = a_n * x^n + ... + a_1 * x + a_0 using Horner's scheme:
    Coefficients passed in descending order: [a_n, a_{n-1}, ..., a_1, a_0].
    Requires only n multiplications and n additions (O(n)).
    """
    result = coefficients[0]
    for i in range(1, len(coefficients)):
        result = result * x + coefficients[i]
    return result


def polynomial_first_derivative(coefficients):
    """
    Returns coefficients of P'(x) given coefficients of P(x) in descending powers:
        d/dx [a_n * x^n + ... + a_1 * x + a_0] = n*a_n * x^{n-1} + ... + a_1
    """
    n = len(coefficients) - 1
    return [coefficients[i] * (n - i) for i in range(n)]


def polynomial_second_derivative(coefficients):
    """
    Returns coefficients of P''(x) given coefficients of P(x) in descending powers.
    """
    first_deriv = polynomial_first_derivative(coefficients)
    return polynomial_first_derivative(first_deriv)


def laguerre(coefficients, b0, accuracy1=1e-8, accuracy2=1e-6, max_iterations=30):
    """
    Laguerre's method to find one real root of polynomial P(x):
        G = P'(x) / P(x)
        H = G^2 - P''(x) / P(x)
        a = n / [ G +/- sqrt((n - 1) * (n*H - G^2)) ]
    The denominator sign is chosen to maximize |denominator| for stability.
    Trial update: b_{k+1} = b_k - a.
    Returns: (root, iterations_taken).
    """
    n = len(coefficients) - 1
    first_deriv = polynomial_first_derivative(coefficients)
    second_deriv = polynomial_second_derivative(coefficients)
    b = float(b0)

    for i in range(max_iterations):
        P = polynomial_value(coefficients, b)
        if abs(P) < accuracy1:
            return b, i + 1

        P1 = polynomial_value(first_deriv, b)
        P2 = polynomial_value(second_deriv, b)

        G = P1 / P
        H = G * G - P2 / P

        radicand = (n - 1) * (n * H - G * G)
        if radicand < 0:
            radicand = 0.0  # Restricted to real roots

        sqrt_term = math.sqrt(radicand)
        denom1 = G + sqrt_term
        denom2 = G - sqrt_term
        denom = denom1 if abs(denom1) >= abs(denom2) else denom2

        if abs(denom) < 1e-14:
            raise ZeroDivisionError("Laguerre denominator vanished.")

        a = n / denom
        b_new = b - a

        if abs(b_new - b) < accuracy1 and abs(polynomial_value(coefficients, b_new)) < accuracy2:
            return b_new, i + 1

        b = b_new

    raise ValueError("Laguerre method did not converge within the maximum iterations.")


def synthetic_division(coefficients, root, accuracy=1e-6):
    """
    Deflates polynomial P(x) by dividing by (x - root) using Horner's synthetic division.
    Given P(x) = [a_n, a_{n-1}, ..., a_0], computes quotient Q(x) of degree n - 1.
    Verifies that the remainder is zero within specified tolerance.
    """
    quotient = [coefficients[0]]
    for i in range(1, len(coefficients) - 1):
        quotient.append(coefficients[i] + root * quotient[-1])

    remainder = coefficients[-1] + root * quotient[-1]
    if abs(remainder) > accuracy:
        raise ValueError(f"Remainder {remainder} is non-zero within tolerance {accuracy}; division failed.")
    return quotient


def laguerre_roots(coefficients, b0=0.0):
    """
    Finds ALL real roots of polynomial P(x) by iteratively applying Laguerre's
    method and synthetic division deflation until degree reduces to 1.
    Returns: (list_of_roots, list_of_iteration_counts).
    """
    roots = []
    total_iterations = []
    curr_coeffs = coefficients.copy()
    guess = float(b0)

    while len(curr_coeffs) > 2:
        root, n_iter = laguerre(curr_coeffs, guess)
        roots.append(root)
        total_iterations.append(n_iter)
        curr_coeffs = synthetic_division(curr_coeffs, root)
        guess = root  # Use current root as warm start for next root

    # Linear remainder: a1 * x + a0 = 0 -> x = -a0 / a1
    final_root = -curr_coeffs[1] / curr_coeffs[0]
    roots.append(final_root)
    total_iterations.append(0)

    return roots, total_iterations


# =============================================================================
# 9. NUMERICAL INTEGRATION & GAUSSIAN QUADRATURE
# =============================================================================

def midpoint(f, a, b, N):
    """
    Midpoint numerical integration (composite rectangle rule):
        I = h * sum_{i=0}^{N-1} f(a + (i + 0.5) * h), where h = (b - a) / N
    Error bound: E_M <= (b - a)^3 / (24 * N^2) * max|f''(x)|.
    """
    h = (b - a) / N
    total = 0.0
    for i in range(N):
        x_mid = a + (i + 0.5) * h
        total += f(x_mid)
    return h * total


def trapezoidal(f, a, b, N):
    """
    Trapezoidal numerical integration:
        I = (h / 2) * [ f(a) + 2 * sum_{i=1}^{N-1} f(a + i*h) + f(b) ]
    Error bound: E_T <= (b - a)^3 / (12 * N^2) * max|f''(x)|.
    """
    h = (b - a) / N
    total = f(a) + f(b)
    for i in range(1, N):
        total += 2.0 * f(a + i * h)
    return (h / 2.0) * total


def simpson(f, a, b, N):
    """
    Simpson's 1/3-rule for numerical integration:
        I = (h / 3) * [ f(a) + 4 * sum_{odd} f(x_i) + 2 * sum_{even} f(x_i) + f(b) ]
    Requirement: N must be an even positive integer (odd number of grid points).
    Error bound: E_S <= (b - a)^5 / (180 * N^4) * max|f''''(x)|.
    Returns: (integral_value, total_function_evaluations).
    """
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    if N % 2 != 0:
        raise ValueError("For Simpson's 1/3 method, N must be even.")

    h = (b - a) / N
    total = f(a) + f(b)
    evaluations = 2

    for i in range(1, N):
        x = a + i * h
        if i % 2 == 0:
            total += 2.0 * f(x)
        else:
            total += 4.0 * f(x)
        evaluations += 1

    integral = (h / 3.0) * total
    return integral, evaluations


def monte_carlo(f, a, b, N, seed=1):
    """
    Uniform Monte Carlo integration in interval [a, b]:
        I = (b - a) * <f>
    Uses self-defined LCG generator myrand().
    Returns: (integral, standard_deviation_sigma_f, standard_error_of_integral).
    Statistical uncertainty decreases as 1 / sqrt(N).
    """
    if N <= 0:
        raise ValueError("N must be a positive integer.")

    myrand(seed)
    sum_f = 0.0
    sum_f2 = 0.0

    for _ in range(N):
        x = a + (b - a) * myrand()
        val = f(x)
        sum_f += val
        sum_f2 += val * val

    mean_f = sum_f / N
    mean_f2 = sum_f2 / N

    variance_f = mean_f2 - mean_f * mean_f
    if variance_f < 0.0:
        variance_f = 0.0

    sigma_f = math.sqrt(variance_f)
    integral = (b - a) * mean_f
    sigma_integral = (b - a) * sigma_f / math.sqrt(N)

    return integral, sigma_f, sigma_integral


def integration_error_bound_N(method, a, b, max_deriv, target_error=1e-6):
    """
    Calculates the minimum number of intervals N required to guarantee that the
    integration error is bounded by target_error (from lecture error bounds):
        Midpoint:    E <= (b - a)^3 / (24 * N^2) * M2  -> N = ceil(sqrt((b - a)^3 * M2 / (24 * eps)))
        Trapezoidal: E <= (b - a)^3 / (12 * N^2) * M2  -> N = ceil(sqrt((b - a)^3 * M2 / (12 * eps)))
        Simpson:     E <= (b - a)^5 / (180 * N^4) * M4 -> N = ceil(((b - a)^5 * M4 / (180 * eps))^(1/4)) (even)
    """
    L = abs(b - a)
    m = method.lower()
    if m in ('midpoint', 'mid'):
        n_exact = math.sqrt((L**3 * max_deriv) / (24.0 * target_error))
        return math.ceil(n_exact)
    elif m in ('trapezoidal', 'trap', 'trapezoid'):
        n_exact = math.sqrt((L**3 * max_deriv) / (12.0 * target_error))
        return math.ceil(n_exact)
    elif m in ('simpson', 'simp'):
        n_exact = ((L**5 * max_deriv) / (180.0 * target_error)) ** 0.25
        n_ceil = math.ceil(n_exact)
        return n_ceil if n_ceil % 2 == 0 else n_ceil + 1
    else:
        raise ValueError("method must be 'midpoint', 'trapezoidal', or 'simpson'.")


# Legendre orthogonal polynomial roots (nodes) and weights on [-1, 1] for n = 1 to 6
GAUSS_LEGENDRE = {
    1: ([0.0], [2.0]),
    2: ([-0.5773502691896257, 0.5773502691896257],
        [1.0, 1.0]),
    3: ([-0.7745966692414834, 0.0, 0.7745966692414834],
        [0.5555555555555556, 0.8888888888888889, 0.5555555555555556]),
    4: ([-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526],
        [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]),
    5: ([-0.9061798459386640, -0.5384693101056831, 0.0,
          0.5384693101056831, 0.9061798459386640],
        [0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
         0.4786286704993665, 0.2369268850561891]),
    6: ([-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969,
          0.6612093864662645, 0.9324695142031521],
        [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910,
         0.3607615730481386, 0.1713244923791704])
}

# Laguerre orthogonal polynomial roots (nodes) and weights on [0, infinity) for n = 1 to 6
# Weight function w(x) = exp(-x)
GAUSS_LAGUERRE = {
    1: ([1.0], [1.0]),
    2: ([0.5857864376269049, 3.414213562373095],
        [0.8535533905932737, 0.1464466094067262]),
    3: ([0.4157745567834791, 2.294280360279041, 6.289945082937479],
        [0.7110930099291730, 0.2785177335692408, 0.0103892565015861]),
    4: ([0.3225476896193923, 1.7457611011583466, 4.536620296921128, 9.395070912301133],
        [0.6031541043416336, 0.3574186924377997, 0.0388879085150054, 0.0005392947055613]),
    5: ([0.2635603197181409, 1.413403059106517, 3.596425771040722, 7.085810005858837, 12.640800844275782],
        [0.5217556105828087, 0.3986668110831759, 0.0759424496817076, 0.0036117586799220, 0.0000233699723858]),
    6: ([0.2228466041792607, 1.188932101672623, 2.992736326059314, 5.775143569104510, 9.837467418382589, 15.982873980601701],
        [0.4589646739499636, 0.4170008307721209, 0.1133733820740449, 0.0103991974531491, 0.0002610172028149, 0.0])
}


def gaussian_quadrature(f, a, b, N, method='legendre'):
    """
    Gaussian Quadrature numerical integration:
        - method='legendre': Integrates f(x) over [a, b] using Gauss-Legendre quadrature.
          Transforms interval [-1, 1] to [a, b] via:
              x = ((b - a)/2) * t + (b + a)/2,  dx = ((b - a)/2) * dt
        - method='laguerre': Integrates f(t) over [0, infinity) with weight exp(-t):
              integral_0^inf e^{-t} * f(t) dt ~ sum_{i=1}^N w_i * f(t_i)
    Supports orders N in {1, 2, 3, 4, 5, 6}.
    """
    m = method.lower()
    if m == 'legendre':
        if N not in GAUSS_LEGENDRE:
            raise ValueError(f"For Gauss-Legendre, N must be one of: {list(GAUSS_LEGENDRE.keys())}")
        nodes, weights = GAUSS_LEGENDRE[N]
        scale = (b - a) / 2.0
        shift = (b + a) / 2.0
        total = 0.0
        for i in range(N):
            x = scale * nodes[i] + shift
            total += weights[i] * f(x)
        return scale * total

    elif m == 'laguerre':
        if N not in GAUSS_LAGUERRE:
            raise ValueError(f"For Gauss-Laguerre, N must be one of: {list(GAUSS_LAGUERRE.keys())}")
        nodes, weights = GAUSS_LAGUERRE[N]
        total = 0.0
        for i in range(N):
            total += weights[i] * f(nodes[i])
        return total

    else:
        raise ValueError("method must be 'legendre' or 'laguerre'.")