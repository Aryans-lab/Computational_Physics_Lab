# Solve the following steady state heat conduction problem using finite element method
# y′′(x) = −2 where x ∈ [0, 1] and y(x = 0) = 1, y(x = 1) = 3
# Compare your solution with the analytical result y(x) = −x^2 + 3x + 1.
# Name: Aryan Bandyopadhyay, Roll number: 2411014. 

import numpy as np
import matplotlib.pyplot as plt
from mylib import cholesky_decomposition, cholesky_forback


def fem_1d(f_rhs, a, b, ya, yb, n_elements):
    h = (b - a) / n_elements
    nodes = [a + i * h for i in range(n_elements + 1)]
    n_int = n_elements - 1

    K = [[0.0] * n_int for i in range(n_int)]
    F = [0.0] * n_int

    for i in range(n_int):
        K[i][i] = 2.0 / h
        if i > 0:
            K[i][i - 1] = -1.0 / h
        if i < n_int - 1:
            K[i][i + 1] = -1.0 / h
        
        # Load vector
        xi = nodes[i + 1]
        F[i] = -f_rhs(xi) * h

    # Use Dirichlet boundary conditions
    F[0] -= (-1.0 / h) * ya
    F[-1] -= (-1.0 / h) * yb

    # Solve K * y_int = F using Cholesky decomposition from mylib
    L = cholesky_decomposition(K)
    y_int = cholesky_forback(L, F)

    # Combine boundary values with interior solutions
    y_sol = [ya] + list(y_int) + [yb]
    return nodes, y_sol


def analytical_function(x):
    return -x**2 + 3.0 * x + 1.0

if __name__ == '__main__':
    a, b = 0.0, 1.0
    ya, yb = 1.0, 3.0
    n_elements = 20  # 20 elements (21 nodes)

    # RHS function f(x) where y''(x) = f(x) = -2
    def f_rhs(x):
        return -2.0

    nodes, y_fem = fem_1d(f_rhs, a, b, ya, yb, n_elements)
    y_exact = [analytical_function(x) for x in nodes]

    print("Finite Element method vs Analytical solution:")
    print("Equation: y''(x) = -2, with y(0) = 1, y(1) = 3")
    print(f"{'Node x':>10} | {'FEM y(x)':>14} | {'Exact y(x)':>14} | {'Absolute Error':>16}")

    for x_val, y_f, y_e in zip(nodes, y_fem, y_exact):
        err = abs(y_f - y_e)
        print(f"{x_val:10.4f} | {y_f:14.6f} | {y_e:14.6f} | {err:15.2e}")

    # Plotting comparison
    plt.figure(figsize=(8, 5))
    x_dense = np.linspace(a, b, 200)
    plt.plot(x_dense, analytical_function(x_dense), 'k-', label='Analytical Solution: $y = -x^2 + 3x + 1$')
    plt.plot(nodes, y_fem, 'ro--', markersize=5, linewidth=1.3, label=f'FEM Solution ({n_elements} elements)')

    plt.title("Steady State Heat Conduction via Finite Element Method (FEM)")
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("finite_element_heat_conduction.png")
    plt.show()