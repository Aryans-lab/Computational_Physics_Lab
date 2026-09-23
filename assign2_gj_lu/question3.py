# Find LU decomposition of the matrix and verify. 
# Matrix: A =[[1,2,4], [3,8,14], [2,6,13]]
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import read_matrix_from_file
from mylib import matrix_multiply


def lu_decomposition(A):
    if not A or not A[0]:
            return []
    
    A = [[float(x) for x in row] for row in A]
    if len(A) != len(A[0]):
            raise ValueError("Matrix A is not a square matrix")
    L = [[float(i == j) for j in range(len(A))] for i in range(len(A))] 
    U = [[float(0) for j in range(len(A))] for i in range(len(A))]
    U[0]=A[0]
    for i in range(len(A)):
        for j in range(i, len(A)):
            U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        if U[i][i] == 0:
            raise ValueError("LU decomposition failed: zero pivot")
        for j in range(i + 1, len(A)):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U

if __name__ == "__main__":
    A = read_matrix_from_file('asgn2_mat2')
    L, U = lu_decomposition(A)
    print("Matrix L:")
    for row in L:
        print(row)
    print("Matrix U:")
    for row in U:
        print(row)
    B = matrix_multiply(L, U)
    print("Final Matrix B is:")
    for row in B:
        print(row)


# End of Code

#############################################################

# Output:
# Matrix L:
# [1.0, 0.0, 0.0]
# [3.0, 1.0, 0.0]
# [2.0, 1.0, 1.0]
# Matrix U:
# [1.0, 2.0, 4.0]
# [0.0, 2.0, 2.0]
# [0.0, 0.0, 3.0]
# Final Matrix B is:
# [1.0, 2.0, 4.0]
# [3.0, 8.0, 14.0]
# [2.0, 6.0, 13.0]

###############################################################