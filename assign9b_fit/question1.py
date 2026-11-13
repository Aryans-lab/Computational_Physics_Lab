# Assignment 9B: Question 1
# Code a library function for standard polynomial of order k for least square fitting
# of data. The order k should be a parameter to be passed to the function.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

from mylib import lu_forback, gauss_jordan_inverse, frobenius_norm

def poly_fit(X, k, sigma=None):

    if not X:
        raise ValueError("Dataset cannot be empty.")

    # Unpack X
    if isinstance(X[0], (list, tuple)) and not isinstance(X[0][0], (list, tuple)) and len(X) == 2:
        x_vals = list(X[0])
        y_vals = list(X[1])
    else:
        x_vals = [row[0] for row in X]
        y_vals = [row[1] for row in X]

    N = len(x_vals)  # Total number of data points
    if N < k + 1:
        raise ValueError(f"At least {k + 1} data points are required for order {k}.")

    # Error weights: w_i = 1 / sigma_i^2 (default is 1)
    if sigma is None:
        w = [1.0] * N
    elif isinstance(sigma, (int, float)):
        w = [1.0 / (float(sigma) ** 2)] * N
    else:
        w = [1.0 / (s ** 2) for s in sigma]

    dim = k + 1
    # Build Hankel matrix (dim x dim) and RHS vector (dim)
    A = [[0.0] * dim for i in range(dim)]
    b = [0.0] * dim

    for i in range(dim):
        b[i] = sum(w[a] * y_vals[a] * (x_vals[a] ** i) for a in range(N))
        for j in range(dim):
            A[i][j] = sum(w[b] * (x_vals[b] ** (i + j)) for b in range(N))

    # Solve A * a = b using LU decomposition
    coeffs = lu_forback(A, b)

    # Compute Condition Number
    B= gauss_jordan_inverse(A)
    cond_num = frobenius_norm(A) * frobenius_norm(B)
    return coeffs, A, cond_num

# End of Code