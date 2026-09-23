# Find the root of the following function in the interval [−1.5, 1.5] to an accuracy of 10−6 using Bisection, 
# Regula Falsi and Newton-Raphson methods. Compare their convergence rate. f(x) = 3x + sin(x) − e^x
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import numpy as np
from mylib import bisection, regula_falsi, newton_raphson

if __name__ == '__main__':
    print("1. Bisection Method:")
    # Define the function
    def f(x):
        return 3*x+np.sin(x)-np.exp(x)

    # Define the interval
    a = -1.5
    b = 1.5

    # Find the root using Bisection method
    root_bisection = bisection(f, a, b, accuracy=1e-6)

    print("\n2. Regula Falsi Method:")
    # Find the root using Regula Falsi method
    root_regula_falsi = regula_falsi(f, a, b, accuracy=1e-6)

    print("\n3. Newton-Raphson Method")
    def df(x):
        return 3 + np.cos(x) - np.exp(x) 

    initial_guess_nr = 0.0
    # Solving with analytical derivative
    root_nr = newton_raphson(f, df=df, x0=initial_guess_nr, accuracy=1e-6)

    # Solving with finite-difference derivative approximation
    print("\nNewton-Raphson using numerical derivative:")
    root_nr_num = newton_raphson(f, df=None, x0=initial_guess_nr, accuracy=1e-6)

# End of Code

#######################################################
# Output:
# 1. Bisection Method:
# The root is approximately 0.3604220151901245
# The number of iterations taken is 23
#
# 2. Regula Falsi Method:
# The root is approximately 0.3604217029603244
# The number of iterations taken is 25
#
# 3. Newton-Raphson Method
# The root is: 0.36042168047601975
# The number of iterations taken is 3
#
# Newton-Raphson using numerical derivative:
# The root is: 0.3604216804760228
# The number of iterations taken is 3
#########################################################