"""
question1.py
------------
Problem : Compute the sum of the first 20 even numbers and the factorial
          of 8 using an explicit loop (no math.factorial or built-ins for
          the core computation).
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER
"""


def sum_of_first_n_even_numbers(n: int) -> int:
    """Return the sum 2 + 4 + 6 + ... + 2n using a loop."""
    total = 0
    for i in range(1, n + 1):
        total += 2 * i
    return total


def factorial(n: int) -> int:
    """Return n! using a loop (no math.factorial)."""
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result


# -- Driver
print("Sum of first 20 even numbers:", sum_of_first_n_even_numbers(20))
print("Factorial of 8              :", factorial(8))

# -- Output
# Sum of first 20 even numbers: 420
# Factorial of 8              : 40320
