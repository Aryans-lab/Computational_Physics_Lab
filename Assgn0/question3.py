# Problem: Take matrix elements from the files asgn0_matA, asgn0_matB, asgn0_vecC and asgn0_vecD (all are ASCII text files) and find AB, BC and D.c
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import os

#manual path where the my matrix files are

base_dir= "D:\\NISER\\Computational Physics Lab\\Warmup" 

def read_matrix_from_file(filename):
    filepath = os.path.join(base_dir, filename)
    with open(filepath, 'r') as f:
        matrix = []
        for line in f:
            line = line.strip()
            if line:
                row = [float(x) for x in line.split()]
                matrix.append(row)
        return matrix

# Reading the matrices from ASCII files
A = read_matrix_from_file('asgn0_matA')
B = read_matrix_from_file('asgn0_matB')
C = read_matrix_from_file('asgn0_vecC')
D = read_matrix_from_file('asgn0_vecD')

def matrix_multiply(X, Y):
    rows_X = len(X)
    col_X = len(X[0])
    rows_Y = len(Y)
    col_Y = len(Y[0])
    #here in the nested loop, I am making all the elementss zero (all elements of the matrix)
    result = [[0.0 for i in range(col_Y)] for i in range(rows_X)]
    for i in range(rows_X):
        for j in range(col_Y):
            for k in range(col_X):
                result[i][j] += X[i][k] * Y[k][j]
    return result

def dot_product_vector(X, Y):
    if len(X) != len(Y):
        print ("Vectors are not of the same length")
    else:
        result=0
        for i in range(len(X)):
            result = result + X[i][0]*Y[i][0]
        return result

# The results:
AB = matrix_multiply(A, B)
BC = matrix_multiply(B, C)
DC = dot_product_vector(D, C)

print("Matrix AB:")
for row in AB:
    print(row)

print ("matrix BC:")
for row in BC:
    print (row)

print ("Dot product of D.C is:", DC)

# End of my code


# Additional note: In the 1,1 position of the matrix AB it seems like I am getting the value error due to float precision. 

# Output:
#####################################################################
# Matrix AB:
# [-0.3000000000000007, -3.5, 5.2]
# [-4.5, -2.0, 4.5]
# [9.3, 0.8, -7.0]
# Matrix BC:
# [1.0]
# [-5.75]
# [-9.0]
# Dot product of D.C is: -3.5
######################################################################
