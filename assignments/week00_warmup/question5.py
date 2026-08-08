"""
question5.py
------------
Problem : Define a MyComplex class and use it to compute the sum, difference,
          product, and modulus of
            c1 = (1.3 - 2.2j)
            c2 = (-0.8 + 1.7j)
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER

Note    : Python's built-in complex type is NOT used; all arithmetic is
          implemented from first principles inside MyComplex.
"""

import math


class MyComplex:
    """
    A minimal complex number class supporting basic arithmetic.

    Attributes
    ----------
    real  : float  - real part
    imag  : float  - imaginary part
    """

    def __init__(self, real: float, imag: float = 0.0) -> None:
        self.real = real
        self.imag = imag

    # -- Display

    def __repr__(self) -> str:
        sign = "+" if self.imag >= 0 else "-"
        return f"({self.real} {sign} {abs(self.imag)}j)"

    # -- Arithmetic

    def __add__(self, other: "MyComplex") -> "MyComplex":
        """Return self + other."""
        return MyComplex(self.real + other.real, self.imag + other.imag)

    def __sub__(self, other: "MyComplex") -> "MyComplex":
        """Return self - other."""
        return MyComplex(self.real - other.real, self.imag - other.imag)

    def __mul__(self, other: "MyComplex") -> "MyComplex":
        """Return self * other using (a+bj)(c+dj) = (ac-bd) + (ad+bc)j."""
        real_part = self.real * other.real - self.imag * other.imag
        imag_part = self.real * other.imag + self.imag * other.real
        return MyComplex(real_part, imag_part)

    # -- Modulus

    def modulus(self) -> float:
        """Return |self| = sqrt(real^2 + imag^2)."""
        return math.sqrt(self.real ** 2 + self.imag ** 2)


# -- Driver

c1 = MyComplex(1.3, -2.2)
c2 = MyComplex(-0.8, 1.7)

print("c1 =", c1)
print("c2 =", c2)
print()
print("Sum        :", c1 + c2)
print("Difference :", c1 - c2)
print("Product    :", c1 * c2)
print(f"|c1|       : {c1.modulus():.10f}")
print(f"|c2|       : {c2.modulus():.10f}")

# -- Output
# c1 = (1.3 - 2.2j)
# c2 = (-0.8 + 1.7j)
#
# Sum        : (0.5 - 0.5j)
# Difference : (2.1 - 3.9j)
# Product    : (2.7 + 3.97j)
# |c1|       : 2.5553864678
# |c2|       : 1.8788294228
