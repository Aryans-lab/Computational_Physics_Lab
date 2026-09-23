# Use appropriate Gaussian quadrature to evaluate the following integral with n =5
# Integral 0 to infinity of e^(-x) / (1 + x) dx, analytical value ≈ 0.59635
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import gaussian_quadrature

if __name__ == "__main__":
    # Define the main function to integrate (excluding the e^(-x) factor)
    def f(x):
        return 1.0 / (1.0 + x)

    # Use 5-point Gaussian quadrature (laguerre method)
    N = 5

    # Using a large upper limit to approximate infinity but it doesn't matter for Laguerre since it is defined for [0, infinity)
    result = gaussian_quadrature(f, 0, 1000, N, method='laguerre') 

    print ("Integral of e^(-x)/(1+x) from 0 to infinity:")
    print("Gaussian Laguerre Quadrature (N=5) =", result) 

# End of Code

################################################################
# Output:
# Integral of e^(-x)/(1+x) from 0 to infinity:
# Gaussian Laguerre Quadrature (N=5) = 0.5950840879689522
#################################################################