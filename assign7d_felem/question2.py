# Question 1. Crank-Nicolson scheme for the PDE: u_t = a u_xx + b u_x + c u
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import matplotlib.pyplot as plt
from mylib import gauss_jordan_inverse, matrix_multiply, matrix_vector_multiply


def crank_nicolson(Nx, Nt, L_min, L_max, T_min, T_max, a, b, c, u0):
    hx = (L_max - L_min) / (Nx +1)
    ht = (T_max - T_min) / Nt

    alpha = (a * ht) / (2.0 * (hx ** 2))
    beta = (b * ht) / (4.0 * hx)
    gamma = (c * ht) / 2.0

    interior_points = Nx - 1

    # Construct Matrix A and B (blank as of now)
    A = [[0.0] * interior_points for i in range(interior_points)]
    B = [[0.0] * interior_points for i in range(interior_points)]
    # Filling both A and B with the matrix elements I got after calculations
    for i in range(interior_points):
        # Diagonal elements
        A[i][i] = 1.0 + 2.0 * alpha - gamma
        B[i][i] = 1.0 - 2.0 * alpha + gamma
        # Lower triangular elements (i-1, j)
        if i > 0:
            A[i][i - 1] = -alpha + beta
            B[i][i - 1] = alpha - beta
        # Upper triangular (i+1, j)
        if i < interior_points - 1:
            A[i][i + 1] = -alpha - beta
            B[i][i + 1] = alpha + beta

    # Find A inverse
    A_inv = gauss_jordan_inverse(A)

    # Calculating A_inv * B
    C = matrix_multiply(A_inv, B)

    # Boundary vector
    b_bc = [0.0] * interior_points
    b_bc[0] = 4.0 * (alpha - beta)
    b_bc[-1] = 6.0 * (alpha + beta)

    # Precompute constant boundary contribution: D = A_inv * b_bc
    D = matrix_vector_multiply(A_inv, b_bc)

    # Initialize temperature array
    u = [[0.0] * (Nx + 1) for i in range(Nt + 1)]
    u[0] = list(u0)
    u[0][0] = 2.0
    u[0][-1] = 3.0

    # Matrix evolution: V_j = M_evol * V_{j-1} + D
    for n in range(Nt):
        previous_interior = u[n][1:-1]
        next_interior = []
        for i in range(interior_points):
            value = 0.0
            for j in range(interior_points):
                value += C[i][j] * previous_interior[j]
            value += D[i]
            next_interior.append(value)

        u[n + 1][1:-1] = next_interior
        u[n + 1][0] = 2.0
        u[n + 1][-1] = 3.0
    return u

if __name__ == "__main__":
    L_min = 0.0
    L_max = 1.0
    T_min = 0.0
    T_max = 1.0
    Nx = 20
    Nt = 5000
    steps = (0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)

    # Params
    a = 2.0
    b = -1.0
    c = 1.0

    hx = (L_max - L_min) / Nx
    x = [L_min + i * hx for i in range(Nx + 1)]

    # Initial condition: u(x, 0) = 2 + x^2
    init_u = [2.0 + position ** 2 for position in x]

    # Solve using Crank-Nicolson scheme
    u_solution= crank_nicolson(Nx, Nt, L_min, L_max, T_min, T_max, a, b, c, init_u)

    for step in steps:
        plt.plot(x, u_solution[step], label=f"t = {step * (T_max - T_min) / Nt:g}")

    plt.xlabel("Position, x")
    plt.ylabel("u(x, t)")
    plt.title("Crank-Nicolson solution of $u_t = a u_{xx} + b u_x + c u$")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("crank_nicolson_lines.png", dpi=600)
    plt.show()

# End of Code
