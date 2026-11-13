# Equation for heat conduction in a thin, un-insulated rod of length L = 10 m is d^2T/dx^2 + α(Ta − T) = 0 
# where the heat transfer coefficient α = 0.01 m^−2 parameterizes heat dissipated to the surrounding air and Ta = 20 deg C is 
# the ambient temperature. If T(x = 0) = 40 deg C and T(x = L) = 200 deg C, solve the boundary value problem to determine at what x the temperature is T = 100 deg C.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import matplotlib.pyplot as plt
from mylib import rk4, lagrange_linear_interpolation


def heat_rod(x, state):
    T, dT = state
    alpha = 0.01
    Ta = 20.0
    return [dT, -alpha * (Ta - T)]


if __name__ == '__main__':
    L = 10.0
    T0 = 40.0
    TL = 200.0
    h_step = 0.01

    # Two initial slope guesses
    z1 = 0.0
    z2 = 2.0

    # Shot 1
    x, sol1 = rk4(heat_rod, 0.0, L, [T0, z1], h=h_step)
    T1 = [s[0] for s in sol1]
    T1_L = T1[-1]

    # Shot 2
    x, sol2 = rk4(heat_rod, 0.0, L, [T0, z2], h=h_step)
    T2 = [s[0] for s in sol2]
    T2_L = T2[-1]

    # Lagrange linear interpolation for correct initial slope
    z_correct = lagrange_linear_interpolation(z1, T1_L, z2, T2_L, TL)

    # Final integration with the correct initial slope
    x, sol_final = rk4(heat_rod, 0.0, L, [T0, z_correct], h=h_step)
    T_final = [s[0] for s in sol_final]

    # Find x where T(x) = 100 deg C using interpolation
    T_target = 100.0
    x_target = None
    for i in range(len(T_final) - 1):
        if (T_final[i] <= T_target <= T_final[i + 1]) or (T_final[i] >= T_target >= T_final[i + 1]):
            x_target = lagrange_linear_interpolation(x[i], T_final[i], x[i + 1], T_final[i + 1], T_target)
            break

    print("Shooting method heat conduction output and iterations")
    print(f"Shot 1: Initial slope z1 = {z1:+.3f} : T(L) = {T1_L:.3f} deg C")
    print(f"Shot 2: Initial slope z2 = {z2:+.2f} : T(L) = {T2_L:.3f} deg C")
    print(f"Interpolated correct slope dT/dx(0) = {z_correct:.4f} deg C/m")
    print(f"Final Shot: T(L) = {T_final[-1]:.3f} deg C (Target = {TL:.2f} deg C)\n")
    print(f"Temperature T = 100 deg C is reached at x = {x_target:.4f} m")

    # Superimposed plot of all
    plt.figure(figsize=(9, 5))
    plt.plot(x, T1, 'r--', label=f'Shot 1 (slope = {z1:+.1f}, T(L) = {T1_L:.1f} °C)')
    plt.plot(x, T2, 'm--', label=f'Shot 2 (slope = {z2:+.1f}, T(L) = {T2_L:.1f} °C)')
    plt.plot(x, T_final, 'b-', linewidth=2, label=f'Final Solution (slope = {z_correct:.2f})')

    # Markers for target temperature and boundaries
    plt.plot(x_target, T_target, 'go', markersize=7, label=f'T = 100 °C at x = {x_target:.4f} m')
    plt.scatter([0, L], [T0, TL], color='black', zorder=5, label='Boundary Conditions')
    plt.axhline(T_target, color='gray', linestyle=':', alpha=0.7)

    plt.title("Temperature Distribution in Rod (Shooting Method + RK4)")
    plt.xlabel("Position x (m)")
    plt.ylabel("Temperature T (°C)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("heat_rod_temperature.png")
    plt.show()