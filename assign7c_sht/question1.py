# Consider a cantilever beam of length L = 2.0 m, square cross-section α = 0.05 × 0.05 m2, uniformly distributed load 
# q = 5000N/m and Young’s modulus Y =2.0 × 10^11 N/m2. The depression of the beam y(x), when fixed at the two ends, is given by
# d^2y/dx^2= −(q/2YI)x(x − L) where I is the moment of inertia. Determine the maximum and half-maximum depressions. I=bh^3/12
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import numpy as np
import matplotlib.pyplot as plt
import mylib
from mylib import rk4


def beam_deflection(t, y):
    L = 2.0  # Length of the beam
    q = 5000.0  # Uniformly distributed load in N/m
    Y = 2.0e11  # Young's modulus in N/m^2
    b = 0.05  # Width of the square cross-sec
    h = 0.05  # Height of the square cross-sec
    I = (b * h**3) / 12.0  # Moment of inertia

    dydt1 = y[1]
    dydt2 = -(q / (2 * Y * I)) * t * (t - L)
    
    return [dydt1, dydt2]


def lagrange_linear_interpolation(x1, y1, x2, y2, y_target):
    return x1 * (y_target - y2) / (y1 - y2) + x2 * (y_target - y1) / (y2 - y1)


if __name__ == '__main__':
    L = 2.0
    h_step = 0.005

    # Shooting method, two initial slope guesses
    z1 = 0.0
    z2 = 0.08

    # Boundary conditions: y(0) = 0 and y(L) = 0

    # Shot 1
    x, sol1 = rk4(beam_deflection, 0.0, L, [0.0, z1], h=h_step)
    y1 = [s[0] for s in sol1]
    y1_L = y1[-1]

    # Shot 2
    x, sol2 = rk4(beam_deflection, 0.0, L, [0.0, z2], h=h_step)
    y2 = [s[0] for s in sol2]
    y2_L = y2[-1]

    # Lagrange linear interpolation for better slope. 
    z_correct = lagrange_linear_interpolation(z1, y1_L, z2, y2_L, 0.0)

    # Final shot with interpolated slope
    x, sol_final = rk4(beam_deflection, 0.0, L, [0.0, z_correct], h=h_step)
    y_final = [s[0] for s in sol_final]

    # Maximum depression (beam deflects downwards, depression = -y)
    min_y = min(y_final)
    max_depression = abs(min_y)
    x_max = x[y_final.index(min_y)]

    # Half max depression
    half_depression = max_depression / 2
    y_half = -half_depression

    # Find the two x positions where depression is half-maximum
    x_half = []
    for i in range(len(y_final) - 1):
        if (y_final[i] >= y_half >= y_final[i + 1]) or (y_final[i] <= y_half <= y_final[i + 1]):
            xh = lagrange_linear_interpolation(x[i], y_final[i], x[i + 1], y_final[i + 1], y_half)
            x_half.append(xh)

    print("Shooting method output and iterations")
    print(f"Shot 1: Initial slope z1 = {z1:+.6f} and Boundary y(L) = {y1_L:+.6f} meters")
    print(f"Shot 2: Initial slope z2 = {z2:+.4f} and Boundary y(L) = {y2_L:+.6f} meters")
    print(f"Interpolated correct slope z = {z_correct:+.6f}")
    print(f"Final Shot: Boundary y(L) = {y_final[-1]:.2e} m (Target = 0.0 m)\n")
    print(f"Maximum depression = {max_depression:.6f} m at x = {x_max:.4f} m")
    print(f"Half-maximum depression = {half_depression * 1000:.4f} mm ({half_depression:.6f} m)")
    print(f"Half-maximum locations = x = {x_half[0]:.4f} m and x = {x_half[1]:.4f} m")

    # All the plots
    plt.figure(figsize=(10, 6))
    plt.plot(x, y1, 'r--', label=f'Shot 1 (slope = {z1:+.2f}, y(L) = {y1_L:+.3f} m)')
    plt.plot(x, y2, 'm--', label=f'Shot 2 (slope = {z2:+.2f}, y(L) = {y2_L:+.3f} m)')
    plt.plot(x, y_final, 'b-', linewidth=2, label=f'Final Solution (slope = {z_correct:+.4f}, y(L) ≈ 0)')

    # Mark max and half-max depressions on final solution
    plt.plot(x_max, min_y, 'ro', markersize=7, label=f'Max Depression ({max_depression*1000:.1f} mm)')
    plt.plot(x_half, [y_half, y_half], 'gs', markersize=6, label=f'Half-Max Depression ({half_depression*1000:.1f} mm)')

    # Boundary and reference lines
    plt.axhline(0, color='k', linestyle='-', linewidth=0.8, alpha=0.7)
    plt.axhline(y_half, color='gray', linestyle=':', alpha=0.7)
    plt.scatter([0, L], [0, 0], color='black', zorder=5, label='Fixed Ends (Boundary Conditions)')

    plt.title("Shooting Method for Beam Deflection (All Shots Superimposed)")
    plt.xlabel("x (m)")
    plt.ylabel("Deflection y (m)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("beam_deflection.png")
    plt.show()
