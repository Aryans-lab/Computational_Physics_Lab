# Prepare library routines for Cholesky decomposition, which must include a check for symmetric matrix.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import mylib
from mylib import read_matrix_from_file
from mylib import lu_decomposition


def lu_forback(A, b):

    C = read_matrix_from_file(A)
    d = read_matrix_from_file(b)
    L=lu_decomposition(A)[0]
    U=lu_decomposition(A)[1]

    #foreword substitution
    y=[]
    for i in range(len(C)):
        sum1=0
        for j in range(i):
            sum1+=L[i][j]*y[j]
        y.append(d[i]-sum1)

    #backward substitution
    x=[]
    for i in range(len(C)-1,-1,-1):
        sum2=0
        for j in range(i+1,len(C)):
            sum2+=U[i][j]*x[len(C)-1-j]
        x.append((y[i]-sum2)/U[i][i])

    return x[::-1]


def cholesky_decomposition(A):
    if not A:
        return []

    A = [[float(x) for x in row] for row in A]

    if any(len(row) != len(A) for row in A):
        raise ValueError("Matrix A is not a square matrix")

    # Check for symmetric Matrix
    for i in range(len(A)):
        for j in range(i):
            if A[i][j] != A[j][i]:
                raise ValueError("Matrix A is not a symmetric matrix")
            
    # Construct the lower-triangular factor row by row.
    L = [[0.0 for i in range(len(A))] for i in range(len(A))]
    for i in range(len(A)):
        for j in range(i + 1):
            value = A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if value <= 0.0:
                    raise ValueError("Matrix A is not positive definite")
                L[i][j] = value**0.5
            else:
                L[i][j] = value / L[j][j]

    return L

# End of Code