# Write library functions for Midpoint and Trapezoidal numerical methods for integration.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

# Midpoint integration

def midpoint(f, a, b, N):
    h = (b - a) / N
    total = 0.0

    # Keep adding the value to total

    for i in range(N):
        x_mid = a + (i + 0.5) * h
        total = total + f(x_mid)
    integral = h * total
    return integral

# Trapezoidal integration

def trapezoidal(f, a, b, N):

    # Dividing into N subintervals
    h = (b - a) / N
    total = 0.0

    # Adding the first and last points
    total = total + f(a) + f(b)

    # Interior points
    for i in range(1, N):
        x = a + i * h
        total = total + 2.0 * f(x)

    integral = (h / 2.0) * total
    return integral

# End of Code