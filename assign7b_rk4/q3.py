# 3. Use RK4 to solve the damped simple harmonic oscillator using the following
# initial values and parameters over a time range [0, 40]. Plot variation of x, v and total energy E with time t.
# x'' + mu * x' + omega^2 * x = 0 with x(0) = 1.0, v(0) = 0.0 and k = 1.0, m = 1.0, mu = 0.15
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import numpy as np
import matplotlib.pyplot as plt
from mylib import rk4

# Params
k = 1.0
m = 1.0
mu = 0.15
omega_sq = k / m

# System of first-order ODEs: dx/dt = v, dv/dt = -mu*v - omega^2*x
def damped_sho(t, state):
    x, v = state
    return [v, -mu * v - omega_sq * x]


if __name__ == '__main__':
    t0 = 0.0
    tf = 40.0
    init_state = [1.0, 0.0]
    h = 0.05

    t, sol = rk4(damped_sho, t0, tf, init_state, h=h)

    x = [s[0] for s in sol]
    v = [s[1] for s in sol]

    # Total energy E = 1/2 * m * v^2 + 1/2 * k * x^2
    E = [0.5 * m * (v[i]**2) + 0.5 * k * (x[i]**2) for i in range(len(t))]

    # Plot variation of x, v, and E with time t
    plt.figure(figsize=(10, 7))

    plt.subplot(3, 1, 1)
    plt.plot(t, x, 'b-', label='Position x(t)')
    plt.ylabel('x')
    plt.title('Damped Harmonic Oscillator (RK4)')
    plt.grid()
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(t, v, 'r-', label='Velocity v(t)')
    plt.ylabel('v')
    plt.grid()
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(t, E, 'g-', label='Total Energy E(t)')
    plt.xlabel('Time t')
    plt.ylabel('Energy E')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.savefig("q3_damped_sho.png")
    plt.show()
