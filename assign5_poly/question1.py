# Write the library functions for Laguarre’s method and deflation for real roots
# only. Pass only the coefficients to the subroutines.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

# ---------------------------------------------------------
# coefficients are given in descending powers:
# [a_n, a_(n-1), ..., a_1, a_0]
# ---------------------------------------------------------

def polynomial_value(coefficients, x):
    result = coefficients[0]
    for i in range(1, len(coefficients)):
        result = result * x + coefficients[i]
    return result

# ---------------------------------------------------------
# First derivative of polynomial with only the coefficients
# ---------------------------------------------------------

def polynomial_first_derivative(coefficients):
    n = len(coefficients) - 1
    derivative = []
    for i in range(n):
        derivative.append(coefficients[i] * (n - i))
    return derivative

# ---------------------------------------------------------
# Second derivative of polynomial
# ---------------------------------------------------------

def polynomial_second_derivative(coefficients):
    first_derivative = polynomial_first_derivative(coefficients)
    second_derivative = polynomial_first_derivative(first_derivative)
    return second_derivative

# ---------------------------------------------------------
# Laguerre's Method. Just finding one root with a guess beta0.
# This formula for Laguerre's method is for real roots only.
# ---------------------------------------------------------

def laguerre(coefficients, b0, accuracy1=1e-8, accuracy2=1e-6, max_iterations=30):
    n = len(coefficients) - 1
    # Calculate the 1st and 2nd derivatives

    first_derivative = polynomial_first_derivative(coefficients)
    second_derivative = polynomial_second_derivative(coefficients)
    # Initialize the guess
    b = float(b0)

    for i in range(max_iterations):
        # Evaluate P(b)
        P = polynomial_value(coefficients, b)

        # If P(b) is already approximately zero,
        # b is the root.
        if abs(P) < accuracy1:
            return b, i+1

        # Evaluate P'(b) and P''(b)
        P1 = polynomial_value(first_derivative, b)
        P2 = polynomial_value(second_derivative, b)

        # Calculate G and H
        G = P1 / P
        H = G * G - P2 / P

        D = (n - 1) * (n * H - G * G)

        # if D < 0:
        #    raise ValueError("The root is complex.")

        # Two possible denominators
        denominator1 = G + D**0.5
        denominator2 = G - D**0.5

        # Choose denominator with larger absolute value
        if abs(denominator1) >= abs(denominator2):
            denominator = denominator1
        else:
            denominator = denominator2

        a = n / denominator
        # New Trial
        b_new = b - a

        # Convergence condition
        if abs(b_new - b) < accuracy1 and abs(polynomial_value(coefficients, b_new)) < accuracy2:
            return b_new, i + 1
        b = b_new


    raise ValueError("Laguerre method did not converge within the maximum number of iterations.")


# ---------------------------------------------------------
# Synthetic Division
# Divides polynomial P(x) by (x - root)
# ---------------------------------------------------------

def synthetic_division(coefficients, root, accuracy=1e-6):
    quotient = []
    # Bring down first coefficient
    quotient.append(coefficients[0])

    # Continue synthetic division
    for i in range(1, len(coefficients) - 1):
        value = coefficients[i] + root * quotient[-1] 
        quotient.append(value)

    # Final value is the remainder
    if abs(coefficients[-1] + root * quotient[-1]) < accuracy:
        remainder = 0.0
    else:
        raise ValueError("Remainder is not zero within the tolerance")
    return quotient

# ---------------------------------------------------------
# Find all real roots using Laguerre's method and deflation
# ---------------------------------------------------------

def laguerre_roots(coefficients, b0=0.0):
    roots = []
    total_iterations = []

    while len(coefficients) > 2:

        # Find one root
        root, n_iter = laguerre(coefficients, b0)

        # Store root and number of iterations
        roots.append(root)
        total_iterations.append(n_iter)

        # Deflate polynomial using the root
        coefficients = synthetic_division(coefficients, root)

        # Use the root just found as the next initial guess
        b0 = root

    # Remaining polynomial is linear: ax + b = 0
    root = -coefficients[1] / coefficients[0]
    roots.append(root)
    total_iterations.append(0)

    return roots, total_iterations

# End of Code