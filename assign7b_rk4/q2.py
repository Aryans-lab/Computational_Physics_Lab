# Use Runge-Kutta-4 method to solve the following differential equation and
# compare with analytical form over the range x in [0, pi/5]. Use step sizes 0.05, 0.15.
# dy/dx = (x+y)^2 with y(0) = 1 and y(x) = tan(x + pi/4) - x
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import numpy as np
import matplotlib.pyplot as plt
from mylib import rk4

def f(x, y):
    return (x + y)**2

def y_analytical(x):
    return np.tan(x + np.pi / 4.0) - x

if __name__ == '__main__':
    a = 0.0
    b = np.pi / 5.0
    y0 = 1.0
    step_sizes = [0.05, 0.15]

    plt.figure(figsize=(8, 5))
    x_dense = np.linspace(a, b, 100)
    plt.plot(x_dense, y_analytical(x_dense), 'k:', label='Analytical')

    for h in step_sizes:
        x_num, y_num = rk4(f, a, b, y0, h=h)
        
        print(f"\nRK4 with step size h = {h}: ")
        plt.plot(x_num, y_num, 'o--', label=f'RK4 h={h}')

    plt.title("dy/dx = (x + y)^2")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("q2_plot.png")
    plt.show()
