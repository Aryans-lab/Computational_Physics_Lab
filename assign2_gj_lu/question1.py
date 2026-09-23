# Prepare library routines for Gauss-Jordan elimination method using row pivoting
# and LU decomposition. For LU mention your choice of Doolittle or Crout clearly.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import read_matrix_from_file

def gauss_jordan_elimination_augmented(augmented_matrix):
    if not augmented_matrix or not augmented_matrix[0]:
        return []
    
    # Making Augmented matrix of float type to avoid integer division
    augmented = [[float(x) for x in row] for row in augmented_matrix]
    
    n = len(augmented)
    m = len(augmented[0])

    if m != n + 1:
        raise ValueError("Augmented matrix must have shape n x (n+1).")

    # Pivoting
    for i in range(n):
        if i == 0:
            pivot_row = max(range(i, n), key=lambda row: abs(augmented[row][i]))
        else:
            pivot_row = i
            while pivot_row < n and augmented[pivot_row][i] == 0:
                pivot_row += 1

        if pivot_row == n or augmented[pivot_row][i] == 0:
            continue

        if pivot_row != i:
            augmented[i], augmented[pivot_row] = augmented[pivot_row], augmented[i]

        pivot = augmented[i][i]
        for j in range(i, n + 1):
            augmented[i][j] /= pivot

        for k in range(n):
            if k != i:
                factor = augmented[k][i]
                for j in range(i, n +1):
                    augmented[k][j] -= factor * augmented[i][j]

    for i in range(n):
        all_zero = True
        for j in range(n):
            if augmented[i][j] != 0:
                all_zero = False
                break

        if all_zero and augmented[i][n] != 0:
            raise ValueError("Incompatible system: no solution exists.")

        if all_zero and augmented[i][n] == 0:
            raise ValueError("Dependent system: infinitely many solutions.")

    return [augmented[i][n] for i in range(n)]


def lu_decomposition(A):
    if not A or not A[0]:
            return []
    
    A = [[float(x) for x in row] for row in A]
    if len(A) != len(A[0]) or any(len(row) != len(A) for row in A):
            raise ValueError("Matrix A is not a square matrix")

    if A[0][0] == 0:
        for i in range(1, len(A)):
            if A[i][0] != 0:
                A[0], A[i] = A[i], A[0]
                break
    
    # Doolittle's method for LU decomposition
    L = [[float(i == j) for j in range(len(A))] for i in range(len(A))] 

    # Empty U matrix initialized with zeros
    U = [[float(0) for j in range(len(A))] for i in range(len(A))]

    # For Doolittle's method, the first row of U is the same as the first row of A
    U[0]=A[0]
    for i in range(len(A)):
        for j in range(i, len(A)):
            U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        if U[i][i] == 0:
            raise ValueError("LU decomposition failed: zero pivot")
        for j in range(i + 1, len(A)):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return [L, U]


# End of Code