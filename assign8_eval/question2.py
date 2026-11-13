# Find the largest eigenvalue and its corresponding eigenvector of the following matrix using Power iteration method.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from mylib import power_iteration, read_matrix_from_file

if __name__ == "__main__":
    # Read matrix from file
    A = read_matrix_from_file("asgn8_mat1")

    # Perform power iteration
    eigenvalue, eigenvector, iterations = power_iteration(A, tolerance=1e-8)

    # Output results
    print("Dominant eigenvalue =", round(eigenvalue, 6))
    print("Corresponding eigenvector =", [round(val, 6) for val in eigenvector])

# End of Code