# Consider the following non-linear equation f(x) = −x − cos(x)
# Find an appropriate interval bracket in which the root possible is, starting with the interval [2, 3].
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import numpy as np
from mylib import regula_falsi

if __name__ == '__main__':
    # Define the function
    def f(x):
        return -x - np.cos(x)

    # Define the interval
    a = 2.0
    b = 3.0 

    # Find the interval bracket where the root is possible

    if f(a) * f(b) >= 0:
        while f(a) * f(b) >= 0:
            m=0.5
            if abs(f(a)) < abs(f(b)):
                a = a - m*(b-a)
            if abs(f(b)) < abs(f(a)):
                b = b + m*(b-a)
                m +=0.1

    print (f"The root is in the interval [{a}, {b}]")    

    # Find the root using Regula Falsi method
    root_regula_falsi = regula_falsi(f, a, b, accuracy=1e-6)

# End of Code

############################################
# Output:
# The root is in the interval [-2.0625, 5.53125]
# The root is approximately -0.7390851332151607
# The number of iterations taken is 22
############################################