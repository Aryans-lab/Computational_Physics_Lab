# Use Laguerre’s method and deflation to find the roots (all real) of the following
# three polynomials.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import laguerre_roots

if __name__ == "__main__":
    # P1(x) = x^4 - x^3 - 7x^2 + x + 6
    P1 = [1, -1, -7, 1, 6]
    # P2(x) = x^4 - 5x^2 + 4
    P2 = [1, 0, -5, 0, 4]
    # P3(x) = 2x^5 - 19.5x^3 + 0.5x^2 + 13.5x - 4.5
    P3 = [2, 0, -19.5, 0.5, 13.5, -4.5]

    #First Polynomial
    roots1, iterations1 = laguerre_roots(P1, b0=1)
    print("P1 roots:")
    for i in range(len(roots1)):
        print("Root", i + 1, "=", roots1[i],
              "Iterations =", iterations1[i])

    #Second Polynomial
    roots2, iterations2 = laguerre_roots(P2, b0=1)
    print("\nP2 roots:")
    for i in range(len(roots2)):
        print("Root", i + 1, "=", roots2[i],
              "Iterations =", iterations2[i])


    # Third Polynomial
    roots3, iterations3 = laguerre_roots(P3, b0=1)
    print("\nP3 roots:")
    for i in range(len(roots3)):
        print("Root", i + 1, "=", roots3[i],
              "Iterations =", iterations3[i])

# End of Code

#################################################
# Output:
# P1 roots:
# Root 1 = 1.0 Iterations = 1
# Root 2 = -0.9999999999999651 Iterations = 5
# Root 3 = -2.0000000000000275 Iterations = 2
# Root 4 = 2.999999999999993 Iterations = 0

# P2 roots:
# Root 1 = 1.0 Iterations = 1
# Root 2 = 1.9999999997991764 Iterations = 3
# Root 3 = -0.9999999991967063 Iterations = 2
# Root 4 = -2.00000000060247 Iterations = 0

# P3 roots:
# Root 1 = 0.5000174745251814 Iterations = 9
# Root 2 = 0.4999825253061549 Iterations = 2
# Root 3 = -0.999999999777354 Iterations = 4
# Root 4 = -3.0000000000272604 Iterations = 2
# Root 5 = 2.9999999999732783 Iterations = 0
######################################################