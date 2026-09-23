# Jacobi and Gauss-Seidel methods for solving linear equations
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import mylib
from mylib import read_matrix_from_file
from mylib import jacobi_it
from mylib import gauss_seidel

if __name__ == "__main__":
    A = read_matrix_from_file('asgn4_mat1')
    b = read_matrix_from_file('asgn4_vec1')
    # Solve the system of linear equations with Jacobi method
    solution = jacobi_it(A, b, tol=1e-6)
    print(f"With jacobi iteration we get the solution: ",solution)
    solution2 = gauss_seidel(A, b, tol=1e-6)
    print(f"With gauss-seidel iteration we get the solution: ",solution2)

# End of Code

###################################################################
# Output:
# Converged after 16 iterations
# With jacobi iteration we get the solution:  [0.9999997530614102, 0.9999997892247294, 0.9999999100460266, 0.9999998509593768, 0.9999998727858708, 0.9999999457079743]
# Converged after 16 iterations
# With gauss-seidel iteration we get the solution:  [0.9999997530614102, 0.9999997892247294, 0.9999999100460266, 0.9999998509593769, 0.9999998727858708, 0.9999999457079743]
####################################################################