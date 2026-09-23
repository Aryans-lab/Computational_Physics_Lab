# Updated for Jacobi, Gauss-Seidel, and SOR methods in the library mylib.py

import mylib
from mylib import read_matrix_from_file

def jacobi_it(A, b, tol=1e-8, max_iterations=1000):
    if not A or not A[0] or not b:
       return []

    A = [[float(x) for x in row] for row in A]
    b = [float(row[0]) for row in b] 

    x=[0 for i in range(len(b))]

    # Swapping and diagonal dominance check (accounts for 0 diagonal)
    for i in range(len(A)):
        if abs(A[i][i]) < sum(abs(A[i][j]) for j in range(len(A)) if j != i):
            for j in range(i + 1, len(A)):
                if abs(A[j][i]) > abs(A[i][i]):
                    A[i], A[j] = A[j], A[i]
                    b[i], b[j] = b[j], b[i]
                    break

    # The main iterative method loop
    for k in range(1, max_iterations+1):
        x_prev = x.copy()
        for i in range(len(b)):
            x[i]=(b[i]-sum(A[i][j]*x[j] for j in range(len(b)) if j!=i))/A[i][i]
        if all(abs(x[i]-x_prev[i]) < tol for i in range(len(b))):
            print(f"Converged after {k} iterations")
            break
        if k == max_iterations:
            raise ValueError("Jacobi method couldn't converge within max iteration limit")

    return x
   
# Gauss Seidel iterative method for solving linear equations Ax = b

def gauss_seidel(A, b, tol=1e-8, max_iterations=500):
    if not A or not A[0] or not b:
       return []

    A = [[float(x) for x in row] for row in A]
    b = [float(row[0]) for row in b] 

    x=[0 for i in range(len(b))]

    # Swapping and diagonal dominance check (accounts for 0 diagonal)
    for i in range(len(A)):
        if abs(A[i][i]) < sum(abs(A[i][j]) for j in range(len(A)) if j != i):
            for j in range(i + 1, len(A)):
                if abs(A[j][i]) > abs(A[i][i]):
                    A[i], A[j] = A[j], A[i]
                    b[i], b[j] = b[j], b[i]
                    print (A)
                    break

    # The main iterative method loop
    for k in range(1, max_iterations+1):
        x_prev = x.copy()
        for i in range(len(b)):
            x[i]=(b[i]-sum(A[i][j]*x[j] for j in range(i))-sum(A[i][j]*x[j] for j in range(i+1, len(b))))/A[i][i]
        if all(abs(x[i]-x_prev[i]) < tol for i in range(len(b))):
            print(f"Converged after {k} iterations")
            break
        if k == max_iterations:
            raise ValueError("Gauss-Seidel method couldn't converge within max iteration limit")

    return x


# Successive Over-relaxation method

def sor_gauss_seidel(A, b, tol=1e-8, max_iterations=500, omega=1):
    if not A or not A[0] or not b:
       return []

    A = [[float(x) for x in row] for row in A]
    b = [float(row[0]) for row in b] 

    x=[0.0 for i in range(len(b))]

    # Swapping and diagonal dominance check
    for i in range(len(A)):
        if abs(A[i][i]) < sum(abs(A[i][j]) for j in range(len(A)) if j != i):
            for j in range(i + 1, len(A)):
                if abs(A[j][i]) > abs(A[i][i]):
                    A[i], A[j] = A[j], A[i]
                    b[i], b[j] = b[j], b[i]
                    break

    # The main iterative method loop
    for k in range(1, max_iterations+1):
        x_prev = x.copy()
        for i in range(len(b)):
            x[i]=omega*(b[i]-sum(A[i][j]*x[j] for j in range(len(b)) if j!=i))/A[i][i] + (1-omega)*x_prev[i]
        if all(abs(x[i]-x_prev[i]) < tol for i in range(len(b))):
            print(f"Converged after {k} iterations")
            break
        if k == max_iterations:
            raise ValueError("Gauss-Seidel method couldn't converge within max iteration limit")

    return x
