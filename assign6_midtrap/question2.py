# Using both Midpoint and Trapezoidal method evaluate the following integrals
# and compare with the analytical results, for N = 4, 8, 15, 20. Present the result
# in tabular form. Integrate: 1/x from 1 to 2, xcos x 0 to pi/2, x tan inv x 0 to 1
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import numpy as np
from mylib import midpoint, trapezoidal

if __name__ == "__main__":

    # Defining all the three functions
    def f1(x):
        return 1 / x
    def f2(x):
        return x * np.cos(x)
    def f3(x):
        return x * np.arctan(x)
    
    # Now compute the integrals using both methods for N = 4, 8, 15, 20 and present in tabular form. 
    print("Function 1: 1/x from 1 to 2")
    print("N\tMidpoint\tTrapezoidal")
    for N in [4, 8, 15, 20]:
        mid = midpoint(f1, 1, 2, N)
        trap = trapezoidal(f1, 1, 2, N)
        print(f"{N}\t{mid:.8f}\t{trap:.8f}")

    print("\nFunction 2: xcos(x) from 0 to pi/2")
    print("N\tMidpoint\tTrapezoidal")
    for N in [4, 8, 15, 20]:
        mid = midpoint(f2, 0, np.pi/2, N)
        trap = trapezoidal(f2, 0, np.pi/2, N)
        print(f"{N}\t{mid:.8f}\t{trap:.8f}")

    print("\nFunction 3: xtan^-1(x) from 0 to 1")
    print("N\tMidpoint\tTrapezoidal")
    for N in [4, 8, 15, 20]:
        mid = midpoint(f3, 0, 1, N)
        trap = trapezoidal(f3, 0, 1, N)
        print(f"{N}\t{mid:.8f}\t{trap:.8f}")


# End of Code


#################################################################
# Output:
# Function 1: 1/x from 1 to 2
# N       Midpoint        Trapezoidal
# 4       0.69121989      0.69702381
# 8       0.69266055      0.69412185
# 15      0.69300843      0.69342480
# 20      0.69306910      0.69330338
#
# Function 2: xcos(x) from 0 to pi/2
# N       Midpoint        Trapezoidal
# 4       0.58744792      0.53760713
# 8       0.57493427      0.56252752
# 15      0.57197166      0.56844624
# 20      0.57145729      0.56947459
# 
# Function 3: xtan^-1(x) from 0 to 1
# N       Midpoint        Trapezoidal
# 4       0.28204605      0.29209835
# 8       0.28456102      0.28707220
# 15      0.28516010      0.28587426
# 20      0.28526426      0.28566596
######################################################################