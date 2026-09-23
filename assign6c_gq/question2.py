# Use appropriate n = 4-point Gaussian quadrature to evaluate
# Integral -1 to 1 of x^2/(1+x^4)dx ~ 0.487495494
# Compare this with Simpson’s 1/3-rule up to the same precision.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import simpson, gaussian_quadrature

if __name__ == "__main__":
    # Define the function to integrate
    def f(x):
        return (x**2) / (1 + x**4)

    # Use 4-point Gaussian quadrature (legendre method)
    N = 4
    result_gaussian = gaussian_quadrature(f, -1.0, 1.0, N, method='legendre')

    # Use Simpson's 1/3-rule
    N_simpson = 8  # Adjusted for the desired precision (If I make it to 10, it will be precise to 3rd decimal place which we don't want)
    result_simpson = simpson(f, -1.0, 1.0, N_simpson)

    print("Integral of x^2/(1+x^4) from -1 to 1:")
    print("Gaussian Quadrature (N=4) =", result_gaussian)
    print("Simpson's Rule (N=8) =", result_simpson[0], "with function evaluations =", result_simpson[1])

# End of Code


#####################################################
# Output:
# Integral of x^2/(1+x^4) from -1 to 1:
# Gaussian Quadrature (N=4) = 0.4816354816354816
# Simpson's Rule (N=8) = 0.4881357142840971 with function evaluations = 9
######################################################
