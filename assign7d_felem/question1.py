# Solve 1+1 dimension heat equation ∂^2u(x, t)/∂x2=∂u(x, t)/∂t on a metal rod of length 2 units, with the initial conditions
# u(0, t) = 0 deg C = u(2, t) for 0 ≤ t ≤ 4
# u(x, 0) = 20 |sin (πx)| deg C for 0 ≤ x ≤ 2
# Use explicit scheme taking number of position grid Nx = 20 and time grid
# Nt = 5000. Show the temperature profile across the length of the rod at time
# steps 0, 10, 20, 50 and 100 in a plot.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import math
import matplotlib.pyplot as plt


def heat_equation(Nx, Nt, L_min, L_max, T_min, T_max, u0):
    # Nx+1 is the total number of spatial points, so there are Nx intervals.
    hx = (L_max - L_min) / (Nx + 1)
    ht = (T_max - T_min) / Nt
    alpha = ht / (hx ** 2)

    if alpha > 0.5:
        raise ValueError(
            f"Explicit scheme is unstable for alpha={alpha:.3g}; "
            "it must be at most 0.5."
        )

    interior_points = Nx - 1
    A = [[0.0] * interior_points for _ in range(interior_points)]
    for i in range(interior_points):
        A[i][i] = 1 - 2 * alpha
        if i > 0:
            A[i][i - 1] = alpha
        if i < interior_points - 1:
            A[i][i + 1] = alpha

    # Initialize the temperature array
    u = [[0.0] * (Nx+1) for i in range(Nt + 1)]
    u[0] = list(u0)
    u[0][0] = 0.0
    u[0][-1] = 0.0

    # Matrix evolution: V_(j+1) = A V_j
    for n in range(Nt):
        previous_interior = u[n][1:-1]
        next_interior = []
        for i in range(interior_points):
            value = 0.0
            for j in range(interior_points):
                value += A[i][j] * previous_interior[j]
            next_interior.append(value)
        u[n + 1][1:-1] = next_interior

    return u

if __name__ == "__main__":
    L_min = 0.0
    L_max = 2.0
    T_min = 0.0
    T_max = 4.0
    Nx = 20
    Nt = 2000
    steps = (0, 10, 20, 50, 100, 200, 500, 1000)

    hx = (L_max - L_min) / Nx
    x = [L_min + i * hx for i in range(Nx + 1)]
    init_temp = [20.0 * abs(math.sin(math.pi * position)) for position in x]
    temperature = heat_equation(Nx, Nt, L_min, L_max, T_min, T_max, init_temp)

    for step in steps:
        plt.plot(x, temperature[step], label=f"t = {step * (T_max - T_min) / Nt:g}")

    plt.xlabel("Position, x")
    plt.ylabel("Temperature, u(x, t) (°C)")
    plt.title("Explicit finite-difference solution of the 1-D heat equation")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# End of Code
# 