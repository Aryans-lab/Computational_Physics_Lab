# Find the root of the following function in the interval [1.5, 3.0] to an accuracy of
# 1e−6 using both Bisection and Regula Falsi method. f(x) = log(x/2) − sin(5x/2)
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import numpy as np
from mylib import bisection, regula_falsi

if __name__ == '__main__':

    # Define the function
    def f(x):
        return np.log(x/2)-np.sin(5*x/2)

    # Define the interval
    a = 0.5
    b = 1.0

    # Find the root using Bisection method
    root_bisection = bisection(f, a, b, accuracy=1e-6)

    # Find the root using Regula Falsi method
    root_regula_falsi = regula_falsi(f, a, b, accuracy=1e-6)

# End of Code

###############################################################
# Output:
# The root is approximately 1.4019301055908207
# The number of iterations taken is 22
# The root is approximately 1.4019299316146128
# The number of iterations taken is 9
#################################################################