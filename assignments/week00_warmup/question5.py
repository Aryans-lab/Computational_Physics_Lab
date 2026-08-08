"""
question5.py
------------
Problem : Define a MyComplex class and compute the sum, difference, product,
          and modulus of two complex numbers read from an external input file.
Usage   : python question5.py <input_file> <output_file>
          e.g.: python question5.py data/q5_input.txt output/q5_output.txt
Input file format:
    real1  imag1    <- c1 = real1 + imag1*j
    real2  imag2    <- c2 = real2 + imag2*j
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER

Note    : Python's built-in complex type is NOT used; all arithmetic is
          implemented from first principles inside MyComplex.
"""

import math
import sys


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

    def __repr__(self) -> str:
        sign = "+" if self.imag >= 0 else "-"
        return f"({self.real} {sign} {abs(self.imag)}j)"

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

    def modulus(self) -> float:
        """Return |self| = sqrt(real^2 + imag^2)."""
        return math.sqrt(self.real ** 2 + self.imag ** 2)


def main() -> None:
    # -- Argument check
    if len(sys.argv) != 3:
        print("Usage: python question5.py <input_file> <output_file>",
              file=sys.stderr)
        sys.exit(1)

    input_file  = sys.argv[1]
    output_file = sys.argv[2]

    # -- Read two complex numbers from file (skip comment lines)
    rows = []
    with open(input_file, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                parts = stripped.split()
                rows.append((float(parts[0]), float(parts[1])))

    c1 = MyComplex(rows[0][0], rows[0][1])
    c2 = MyComplex(rows[1][0], rows[1][1])

    # -- Compute all operations
    c_sum  = c1 + c2
    c_diff = c1 - c2
    c_prod = c1 * c2
    mod1   = c1.modulus()
    mod2   = c2.modulus()

    # -- Assemble output text
    output_text = (
        f"c1 = {c1}\n"
        f"c2 = {c2}\n"
        f"\n"
        f"Sum        : {c_sum}\n"
        f"Difference : {c_diff}\n"
        f"Product    : {c_prod}\n"
        f"|c1|       : {mod1:.10f}\n"
        f"|c2|       : {mod2:.10f}\n"
    )

    # -- Write to output file
    with open(output_file, "w") as out:
        out.write(output_text)

    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()

# -- Output (for c1=1.3-2.2j, c2=-0.8+1.7j)
# c1 = (1.3 - 2.2j)
# c2 = (-0.8 + 1.7j)
#
# Sum        : (0.5 - 0.5j)
# Difference : (2.1 - 3.9j)
# Product    : (2.7 + 3.97j)
# |c1|       : 2.5553864678
# |c2|       : 1.8788294228
