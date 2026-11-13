# Wein's displacement law
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import math
import mylib
from mylib import newton_raphson


if __name__ == "__main__":

    # Given Function is- f(x) = (x - 5) * e^x + 5 = 0
    def f(x):
       return (x - 5.0) * math.exp(x) + 5.0

    # First derivative: f'(x) = (x - 4) * e^x
    def df(x):
        return (x - 4.0) * math.exp(x)
    
    # Given constants
    h = 6.626e-34 # in m^2 kg / s
    k = 1.381e-23 # in m^2 kg / s^2 K
    c = 3.0e8    # m / s

    # Initial guess
    x0 = 6.0
    x = newton_raphson(f, df, x0=x0, accuracy=1e-4)

    # Calculate Wien's displacement constant b = hc / (k * x)
    b = (h * c) / (k * x)

    print(f"\nRoot x = {x:.6f}")
    print(f"Calculated Wien's constant b (within the tolerance of 1e-4)= {b:.6e} m K")

# End of code

#######################################################################################
# Output:
# The root is: 4.965114942713781
# The number of iterations taken is 5

# Root x = 4.965115
# Calculated Wien's constant b (within the tolerance of 1e-4)= 2.899010e-03 m K
########################################################################################