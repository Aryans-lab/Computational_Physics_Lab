# Assignment 9B: Question 2
# Write a library function for least square fitting of data using special polynomials,
# Chebyshev and Legendre to n-th order. For now use the modified Chebyshev given below:
# phi_0(x) = 1, phi_1(x) = 2x - 1, phi_2(x) = 8x^2 - 8x + 1, phi_3(x) = 32x^3 - 48x^2 + 18x - 1
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

from mylib import lu_forback

def modified_chebyshev_basis(m, x):
    # hardcoded modified Chebyshev basis functions on [0, 1] (for now)
    if m == 0:
        return 1.0
    elif m == 1:
        return 2.0 * x - 1.0
    elif m == 2:
        return 8.0 * (x ** 2) - 8.0 * x + 1.0
    elif m == 3:
        return 32.0 * (x ** 3) - 48.0 * (x ** 2) + 18.0 * x - 1.0
    else:
        raise ValueError(f"Modified Chebyshev basis currently defined for m in [0, 3], got {m}.")


def chebyshev_fit(X, k=3, sigma=None):

    # Fits data X using modified Chebyshev (for now) basis functions up to order k (till k=3 for now) 

    if not X:
        raise ValueError("Dataset cannot be empty.")

    if isinstance(X[0], (list, tuple)) and not isinstance(X[0][0], (list, tuple)) and len(X) == 2:
        x_vals = list(X[0])
        y_vals = list(X[1])
    else:
        x_vals = [row[0] for row in X]
        y_vals = [row[1] for row in X]

    N = len(x_vals)
    if N < k + 1:
        raise ValueError(f"At least {k + 1} data points are required for order {k}.")

    if sigma is None:
        w = [1.0] * N
    elif isinstance(sigma, (int, float)):
        w = [1.0 / (float(sigma) ** 2)] * N
    else:
        w = [1.0 / (s ** 2) for s in sigma]

    dim = k + 1
    # Hardcoded basis values
    phi_vals = [[modified_chebyshev_basis(i, j) for i in range(dim)] for i in x_vals]

    # Build regression matrix A and RHS vector b
    A = [[0.0] * dim for i in range(dim)]
    b = [0.0] * dim

    for i in range(dim):
        b[i] = sum(w[a] * y_vals[a] * phi_vals[a][i] for a in range(N))
        for j in range(dim):
            A[i][j] = sum(w[a] * phi_vals[a][i] * phi_vals[a][j] for a in range(N))

    # Solve A * a = b using LU decomposition
    coeffs = lu_forback(A, b)
    return coeffs, A

# End of Code