# In QM, radial wave function for ground state H atom is R10=(2/a0^(3/2))*exp(-r/a0) where a0 is Bohr radius. 
# Find the expectattion value <r> of electron's distance from its nucleus
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import math
import mylib
from mylib import gaussian_quadrature

def f(u):
    return 4.0 * (u**3) * math.exp(-2.0 * u)

if __name__ == "__main__":
    # Upper integration limit in units of a0 where the wave function has effectively decayed to 0
    b = 6.0
    N = 6       # 6-point Gauss-Legendre quadrature

    # Numerical integration using Gauss-Legendre quadrature
    result = gaussian_quadrature(f, 0.0, b, N=N, method='legendre')

    print("Numerical evaluation of <r> using Gauss-Legendre Quadrature:")
    print(f"Gauss-Legendre (N = {N}, upper limit = {b} a0):")
    print(f"<r> = {result:.4f} a0")


# End of code

##############################################
# Output:
# Numerical evaluation of <r> using Gauss-Legendre Quadrature:
# Gauss-Legendre (N = 6, upper limit = 6.0 a0):
# <r> = 1.5004 a0
#################################################