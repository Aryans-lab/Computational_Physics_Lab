# Write library function for Simpson’s 1/3-rule and Monte Carlo integration by passing function.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import mylib
from mylib import myrand

# Simpson's 1/3 method
# N must be even

def simpson(f, a, b, N):
    if N <= 0:
        raise ValueError("N must be a positive integer.")

    if N % 2 != 0:
        raise ValueError("For Simpson's 1/3 method, N must be even.")

    #Computing the step size
    h = (b - a)/N

    # Defining the total sum of function evaluations  with the first and last terms
    total = f(a) + f(b)
    function_evaluations = 2

    # Interior points calculation
    for i in range(1, N):
        x = a + i * h
        if i % 2 == 0:
            total = total + 2.0 * f(x)
        else:
            total = total + 4.0 * f(x)
        function_evaluations += 1

    integral = (h / 3.0) * total

    return integral, function_evaluations


# Monte Carlo integration
# Uses the LCG random number generator myrand() from mylib.py

def monte_carlo(f, a, b, N, seed=1):
    if N <= 0:
        raise ValueError("N must be a positive integer.")

    # Initialise the pRNG
    myrand(seed)

    sum_f = 0.0
    sum_f2 = 0.0

    for i in range(N):
        # myrand() gives a uniform random number in [0, 1)
        xi = myrand()

        # Convert it to [a, b]
        x = a + (b - a) * xi

        value = f(x)

        sum_f = sum_f + value
        sum_f2 = sum_f2 + value * value

    mean_f = sum_f / N
    mean_f2 = sum_f2 / N

    sigma_f_squared = mean_f2 - mean_f * mean_f

    # Protect against tiny negative values caused by round-off
    if sigma_f_squared < 0.0:
        sigma_f_squared = 0.0

    sigma_f = sigma_f_squared ** 0.5

    integral = (b - a) * mean_f

    # Standard error of the integral estimate
    sigma_integral = (b - a) * sigma_f / (N ** 0.5)

    return integral, sigma_f, sigma_integral

# End of Code