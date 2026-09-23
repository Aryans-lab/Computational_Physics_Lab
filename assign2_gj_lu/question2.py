# Check Gauss-Jordon with the system of equations for three variables.
# 2y + 5z = 1, 3x − y + 2z = −2, x − y + 3z = 3. The matrix is defined in asgn2_mat1 file
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import read_matrix_from_file

def gauss_jordan_elimination_augmented(augmented_matrix):
    if not augmented_matrix or not augmented_matrix[0]:
        return []

    augmented = [[float(x) for x in row] for row in augmented_matrix]

    # Making Augmented matrix of float type to avoid integer division
    augmented = [[float(x) for x in row] for row in augmented_matrix]

    # Pivoting
    if augmented[0][0] == 0:
        for i in range(len(augmented)):
            if augmented[i][0] != 0:
                augmented[0], augmented[i] = augmented[i], augmented[0]
                break

    n = len(augmented)
    m = len(augmented[0])

    if m != n + 1:
        raise ValueError("Augmented matrix must have shape n x (n+1).")

    for i in range(n):
        pivot_row = i
        while pivot_row < n and augmented[pivot_row][i] == 0:
            pivot_row += 1

        if pivot_row == n:
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


if __name__ == "__main__":
    A = read_matrix_from_file('asgn2_mat1')
    solution = gauss_jordan_elimination_augmented(A)
    print ("The values of x, y, z are: ", solution[0], solution[1], solution[2])

# End of code

#############################################################
# Output:
# The values of x, y, z are:  -2.0 -2.0 1.0
#############################################################