
# Prepare a library function for Power Iteration method to find eigenvalues and eigenvectors of a matrix when it is passed to it.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from mylib import matrix_vector_multiply, vector_norm, dot_product_vector, myrand

def power_iteration(A, tolerance=1e-2, max_iterations=1000, x0=None, y=None):
    if not A or any(len(row) != len(A) for row in A):
        raise ValueError("Matrix A must be a non-empty square matrix.")

    n = len(A)
    if x0 is None:
        # Fill x0 with random values in (0.5, 1.5)
        x0 = [0.5 + myrand() for i in range(n)]
    if y is None:
        # Defaulting y to x0 guarantees a non-zero initial dot product
        y = list(x0)

    # Normalize initial vector
    norm_x0 = vector_norm(x0)
    if norm_x0 == 0:
        raise ValueError("Initial vector x0 cannot be the zero vector.")
    xk = [val / norm_x0 for val in x0]

    old_lambda = None
    for i in range(1, max_iterations + 1):
        x_next = matrix_vector_multiply(A, xk)

        denominator = dot_product_vector(xk, y)
        if abs(denominator) < 1e-15:
            raise ValueError("The denominator (xk . y) became zero during iteration.")

        new_lambda = dot_product_vector(x_next, y) / denominator

        # Normalize x_next at EVERY step to prevent
        norm = vector_norm(x_next)
        if norm == 0:
            return 0.0, xk, i
        xk = [value / norm for value in x_next]

        if old_lambda is not None and abs(new_lambda - old_lambda) < tolerance:
            return new_lambda, xk, i

        old_lambda = new_lambda

    raise RuntimeError("Power iteration did not converge within max_iterations.")

# End of Code