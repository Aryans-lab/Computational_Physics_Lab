# A 2m long beam has a linear mass density lambda (x)=x^2 where x is measured from one its ends. Find the centre of mass 
# of the beam numerically, correct upto 4 decimal places using a method of your choice.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import mylib
from mylib import midpoint

if __name__ == "__main__":
    L = 2.0     # Length of beam in meters
    N1 = 116   # Number of intervals, computed using error bound formula for midpoint, for density function
    N2= 200   # Number of intervals, computed using error bound formula for midpoint, for moment fn

    # Linear mass density: lambda(x) = x^2
    def density(x):
        return x**2
    
    # Moment integration: x * lambda(x) = x^3
    def moment(x):
        return x**3

    # Total mass M = integral of lambda(x) dx from 0 to L
    mass = midpoint(density, 0, L, N1)

    # Moment of mass = integral of x * lambda(x) dx from 0 to L
    Moment = midpoint(moment, 0, L, N2)

    # Centre of mass should be x_cm = (integral x * lambda dx) / (integral lambda dx)
    x_cm = Moment / mass

    print("Numerical calculation using Midpoint method (N1 = 116, N2 = 200):")
    print(f"Centre of mass (x_cm)= {x_cm:.4f} m")


# End of code

##########################################
# Output:
# Numerical calculation using Midpoint method (N1=116, N2=200):
# Centre of mass (x_cm)= 1.5000 m
############################################

#
