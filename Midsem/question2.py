# Find the inverse of the matrix using LU decompositions. Round off your result to 3 decimal places
# Name: Aryan Bandyopadhyay, Roll No: 2411014

from mylib import read_matrix_from_file, matrix_determinant_lu, lu_forback

def lu_inverse(A):
   
    if isinstance(A, str):
        A_mat = read_matrix_from_file(A)
    else:
        A_mat = [[float(x) for x in row] for row in A]

    n = len(A_mat)
    if any(len(row) != n for row in A_mat):
        raise ValueError("Matrix A must be square to compute inverse.")

    det = matrix_determinant_lu(A_mat)
    if abs(det) < 1e-12: # Setting the tolerance  
        raise ValueError(f"Matrix is singular (det(A) = {det}); inverse does not exist.")
    inv_cols = []
    for j in range(n):
        ej = [1.0 if i == j else 0.0 for i in range(n)]
        xj = lu_forback(A_mat, ej, method='lu')
        inv_cols.append(xj)

    # Transpose column vectors into rows
    A_inv = [[inv_cols[j][i] for j in range(n)] for i in range(n)]
    return A_inv

if __name__ == "__main__":
    A = read_matrix_from_file('midsem_mat1')

    det = matrix_determinant_lu(A)
    print(f"Determinant of A: {det:.3f}")

    if abs(det) < 1e-9:
        print("Determinant is 0, so inverse does not exist.")
    else:
        print("Determinant is non-zero, so inverse exists.\n")
        A_inv = lu_inverse(A)

        print("Inverse of matrix A (rounded to 3 decimal places):")
        for row in A_inv:
            print([round(x, 3) for x in row])

# End of Code

##########################################################
# Output:
# Determinant of A: 68.7168
# Determinant is non-zero, so inverse exists.
#
# Inverse of matrix A (rounded to 3 decimal places):
# [-0.708, 2.531, 2.431, 0.967, -3.902]
# [-0.193, 0.31, 0.279, 0.058, -0.294]
# [0.022, 0.365, 0.286, 0.051, -0.29]
# [0.273, -0.13, 0.132, -0.141, 0.449]
# [0.782, -2.875, -2.679, -0.701, 4.234]
############################################################