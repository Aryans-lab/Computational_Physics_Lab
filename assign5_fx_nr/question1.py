# Write the library functions for Fixed-point and Newton-Raphson methods.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from mylib import gauss_jordan_inverse

# 1. Fixed Point Iteration Method

def fixed_point(g, x0, accuracy=1e-6, max_iterations=30):
    """
    Finds the root of f(x) = 0 by finding the fixed point of x = g(x).
    Convergence requires |g'(x)| < 1 near the root.
    """
    x = x0
    for i in range(max_iterations):
        x_next = g(x)
        
        # Stopping criteria: change between successive iterations
        if abs(x_next - x) < accuracy:
            print(f"The root is: {x_next}")
            print(f"The number of iterations taken is {i + 1}")
            return x_next
        x = x_next

    print("Fixed-point method did not converge within the maximum iterations.")
    print(f"The last computed value is: {x}")
    return x, max_iterations

# 2. Newton-Raphson Method

def newton_raphson(f, df=None, x0=0.0, accuracy=1e-6, max_iterations=30, h=1e-5):
    """
    Finds the root of f(x) = 0 using the Newton-Raphson formula:
        x_{n+1} = x_n - f(x_n) / f'(x_n)
    If analytical derivative df is not provided, central difference is used.
    """
    x = x0
    for i in range(max_iterations):
        fx = f(x)
        
        # Check if current guess already satisfies accuracy
        if abs(fx) < accuracy:
            print(f"The root is: {x}")
            print(f"The number of iterations taken is {i}")
            return x

        # Calculate derivative: analytical or central difference
        if df is not None:
            dfx = df(x)
        else:
            dfx = (f(x + h) - f(x - h)) / (2.0 * h)

        if abs(dfx) < 1e-12:
            print("Error: Derivative too close to zero; method fails.")
            return None

        # Newton-Raphson update step
        x_next = x - (fx / dfx)

        # Convergence criteria
        if abs(x_next - x) < accuracy and abs(f(x_next)) < accuracy:
            print(f"The root is: {x_next}")
            print(f"The number of iterations taken is {i + 1}")
            return x_next

        x = x_next

    print("Newton-Raphson did not converge within maximum iterations.")
    return x, max_iterations

# 3. Newton-Raphson for multivariable systems

def newton_raphson_system(f, J, initial_guess, accuracy=1e-6, max_iterations=30):
        x = initial_guess
        for i in range(max_iterations):
            fx = f(*x)
            norm_fx = sum([abs(val) for val in fx])
            if norm_fx < accuracy:
                return x, i
            Jx = J(*x)
            J_inv = gauss_jordan_inverse(Jx)
            x_next = [x[j] - sum(J_inv[j][k] * fx[k] for k in range(len(fx))) for j in range(len(x))]
            x = x_next
        return x, max_iterations

# End of Code