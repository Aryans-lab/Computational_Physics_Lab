# Solve the following using either Newton-Raphson or Fixed-point. 
# f1(x, y) = 3x^2 − 2xy − 3
# f2(x, y) = 3y^2 − 4xy
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from mylib import gauss_jordan_inverse

def partial_derivative(f, var_index, point, h=1e-5):
    point_forward = list(point)
    point_backward = list(point)
    point_forward[var_index] += h
    point_backward[var_index] -= h
    return (f(*point_forward) - f(*point_backward)) / (2 * h)

# Compute the Jacobian matrix of a vector-valued function f at a given point using central difference method. Store in a nested list the jacobian
def jacobian(f, point, h=1e-6):
    n = len(point)
    m = len(f(*point))
    J = [[0.0 for _ in range(n)] for _ in range(m)]
    for i in range(m):
        for j in range(n):
            def fi(*args):
                return f(*args)[i]
            J[i][j] = partial_derivative(fi, j, point, h)
    return J

# Newton Raphson multivariable
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

if __name__ == '__main__':

    # Define the system of equations
    def f(x, y):
        return [3*x**2 - 2*x*y - 3, 3*y**2 - 4*x*y]

    # Define the Jacobian matrix
    def J(x, y):
        return jacobian(f, [x, y])

    # Initial guess
    initial_guess = [2.5, 3.1]
    x, max_iterations = newton_raphson_system(f, J, initial_guess)
    print(f"Final root approximation: ({x[0]:.4f}, {x[1]:.4f}) after {max_iterations} iterations.")

# End of Code

#################################################################
# Output:
# Final root approximation: (3.0000, 4.0000) after 4 iterations.
##################################################################
