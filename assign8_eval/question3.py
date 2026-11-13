# Consider the 5 × 5 tridiagonal matrix given in the file assign8.txt. Using Power iteration method, verify that the first two largest 
# eigenvalues and their corresponding eigenvectors of the matrix satisfy λk = b + 2√(ac) cos(kπ/(n + 1)), and the i-th component of the k-th 
# eigenvector is given by vk,i = sin(ikπ/(n + 1)), where a = c = −1, b = 2, n = 5, k = 1, 2, . . . 5 and i is the i-th component of
# the k-th eigenvector. In case of any discrepancies, discuss its possible source(s).
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from mylib import power_iteration, read_matrix_from_file

if __name__ == "__main__":
    # Read matrix from file
    A = read_matrix_from_file("asgn8_mat2")

    # Perform power iteration for the first largest eigenvalue
    eigenvalue, eigenvector, iterations = power_iteration(A, tolerance=1e-8)

    # Output results for the first largest eigenvalue
    print("First dominant eigenvalue =", round(eigenvalue, 6))
    print("Corresponding eigenvector =", [round(val, 6) for val in eigenvector])

    # Perform power iteration for the second largest eigenvalue
    # To find the second largest eigenvalue, we can deflate the matrix

    n = len(A)
    outer_product = [[eigenvector[i] * eigenvector[j] for j in range(n)] for i in range(n)]
    A_deflated = [[A[i][j] - eigenvalue * outer_product[i][j] for j in range(n)] for i in range(n)]

    # Now perform power iteration on the deflated matrix to find the second largest eigenvalue
    eigenvalue2, eigenvector2, iterations2 = power_iteration(A_deflated, tolerance=1e-8)

    # Output results for the second largest eigenvalue
    print("Second dominant eigenvalue =", round(eigenvalue2, 6))
    print("Corresponding eigenvector =", [round(val, 6) for val in eigenvector2])

# End of Code