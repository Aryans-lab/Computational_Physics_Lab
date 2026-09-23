# Solve the following linear equation by Gauss-Seidel and SOR with ω = 1.57 for
# n = 10 and compare the convergence rate.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import mylib
from mylib import read_matrix_from_file
from mylib import sor_gauss_seidel

if __name__ == "__main__":
    A =[]
    order=10
    for i in range(order):
        A.append([])
        for j in range(order):
            if i==j:
                A[i].append(2)
            elif j==i+1 or j==i-1:
                A[i].append(-1)
            else:
                A[i].append(0)
    b= [[1] for i in range(order)]
    solution = sor_gauss_seidel(A, b, omega=1.57, tol=1e-16)
    print(f"With SOR method we get the solution: ",solution)


# End of Code

###################################################
# Output: 
# Converged after 71 iterations
# With SOR method we get the solution:  [5.0, 9.0, 12.0, 14.0, 15.0, 15.0, 14.0, 12.0, 9.0, 5.0]
###################################################