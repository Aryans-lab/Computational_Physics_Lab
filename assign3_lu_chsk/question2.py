# Given below is a system of linear equations. Use the LU decomposition to solve
# it by forward-backward substitution. Use either Doolittle or Crout.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import mylib
from mylib import lu_forback

if __name__ == "__main__":
    solution= lu_forback('asgn3_mat1','asgn3_vec1')
    print("The solution of the system of linear equations is:")
    for i in range(len(solution)):
        print(f"a{i+1} = {solution[i]}")

# End of Code