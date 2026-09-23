# Problem: Define MyComplex class / structure and calculate the sum, difference, product and modulus of (1.3 − 2.2j) and (−0.8 + 1.7j).
# Name: Aryan Bandyopadhyay, Roll number: 2411014

class MyComplex():
    def __init__(self, real, image=0.0):
        self.real = real
        self.image = image
    def display_complex(self):
        print(f"{self.real} + {self.image}j")
    def add_complex(self, c1, c2):
        return MyComplex(c1.real + c2.real, c1.image + c2.image)
    def subtract_complex(self, c1, c2):
        return MyComplex(c1.real - c2.real, c1.image - c2.image)
    def multiply_complex(self, c1, c2):
        real_part = c1.real * c2.real - c1.image * c2.image
        image_part = c1.real * c2.image + c1.image * c2.real
        return MyComplex(real_part, image_part)
    def modulus_complex(self, c):
        return (c.real*c.real + c.image*c.image)**0.5

# now creating the complex numbers:
c1 = MyComplex(1.3, -2.2)
c2 = MyComplex(-0.8, 1.7)

# Output:
print ("Sum of the complex numbers:")
c_sum = c1.add_complex(c1, c2)
c_sum.display_complex()

print("Difference of the complex numbers:")
c_diff = c1.subtract_complex(c1, c2)
c_diff.display_complex()

print("Product of the complex numbers:")
c_prod = c1.multiply_complex(c1, c2)
c_prod.display_complex()

print("Modulus of the first complex number:")
print(c1.modulus_complex(c1))

print("Modulus of the second complex number:")
print(c2.modulus_complex(c2))


# End of my code


# Output:
#####################################################################
# Sum of the complex numbers:
# 0.5 + -0.5000000000000002j
# Difference of the complex numbers:
# 2.1 + -3.9000000000000004j
# Product of the complex numbers:
# 2.7 + 3.97j
# Modulus of the first complex number:
# 2.5553864678361276
# Modulus of the second complex number:
# 1.8788294228055935
######################################################################