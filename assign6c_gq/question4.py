# Integrate the following with appropriate n-point Gauss-Laguerre.
# Integral 0 to infinity of e^(-x) x^4 dx, compare with the analytical value of n! = 24
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import gaussian_quadrature

if __name__ == "__main__":
    # Define the main function to integrate (excluding the e^(-x) factor)
    def f(x):
        return x**4

    # We have a polynomial of degree 4, so if we use 3 point polynomial rule
    # if we use N=3, we will get 2N-1=5, which is greater than 4 so we will get the exact value
    N = 3

    # Using a large upper limit to approximate infinity but it doesn't matter for Laguerre since it is defined for [0, infinity)
    result = gaussian_quadrature(f, 0, 1000, N, method='laguerre') 
    result = result
    result = round(result, 8)  # Round to 8 decimal places for comparison

    print ("Integral of e^(-x)x^4 from 0 to infinity:")
    print("Gaussian Laguerre Quadrature (N=3) =", result) 
    print("Analytical value (4!) =", 24)
    print("Difference =", abs(result - 24))
    print("Relative error =", abs(result - 24)/24)

# End of Code

#################################################################
# Output:
# Integral of e^(-x)x^4 from 0 to infinity:
# Gaussian Laguerre Quadrature (N=3) = 24.0
# Analytical value (4!) = 24
# Difference = 0.0
# Relative error = 0.0
###################################################################