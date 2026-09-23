# Use Cholesky factorization for the matrix A and solve for Ax = b for the following
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import mylib
from mylib import cholesky_decomposition
from mylib import cholesky_forback
from mylib import read_matrix_from_file

if __name__ == "__main__":
    A=read_matrix_from_file('asgn3_mat2')
    b=read_matrix_from_file('asgn3_vec2')
    L= cholesky_decomposition(A)
    print("The Cholesky factorization of the matrix A is:")
    for row in L:
        print(row)  
    print("The solution of the system of linear equations is:")
    solution= cholesky_forback(L, b)
    for i in range(len(solution)):
        print(f"a{i+1} = {solution[i]}")

# End of Code