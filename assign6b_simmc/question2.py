# Evaluate the following integrals using Midpoint and Simpson’s 1/3-rule accurate
# up to 6 places in decimal. Use the error bound formulae to determine N.
# f1: integral of 1/x from 1 to 2, f2: integral from 0 to pi/2 of xcosx
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import math
import mylib
from mylib import simpson, midpoint

if __name__ == "__main__":
    # Example 1:
    # integral of 1/x from 1 to 2
    def f1(x):
        return 1.0 / x

    # The value of N for Midpoint and Simpson's rule is determined based on the error bound formulae. 
    # For this example, we have chosen N=289 for Midpoint and N=20 for Simpson's rule to achieve the desired accuracy of 1e-6
    
    midpoint_result = midpoint(f1, 1.0, 2.0, 289)
    simpson_result = simpson(f1, 1.0, 2.0, 20)

    print("Integral of 1/x from 1 to 2")
    print("Midpoint =", midpoint_result)
    print("Simpson =", simpson_result[0], "with function evaluations =", simpson_result[1])
    print()

    # Example 2:
    # integral from 0 to pi/2 of x*cos(x) dx
    def f2(x):
        return x * math.cos(x)

    # For this example, we have chosen N=610 for Midpoint and N=22 for 
    # Simpson's rule to achieve the desired accuracy of 1e-6

    midpoint_result = midpoint(f2, 0.0, math.pi / 2.0, 610)
    simpson_result = simpson(f2, 0.0, math.pi / 2.0, 22)

    print("Integral of x*cos(x) from 0 to pi/2")
    print("Midpoint =", midpoint_result)
    print("Simpson =", simpson_result[0], "with function evaluations =", simpson_result[1])

# End of Code