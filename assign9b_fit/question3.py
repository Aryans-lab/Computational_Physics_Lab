# Assignment 9B: Question 3
# The data from an experiment is given in the file assign9fit.txt. Try a cubic least square fit using
# y = f(x) = a0 + a1*x + a2*x^2 + a3*x^3
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import os
import matplotlib.pyplot as plt
from mylib import poly_fit

def read_data(filename):
    if not os.path.exists(filename):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        alt_path = os.path.join(script_dir, filename)
        if os.path.exists(alt_path):
            filename = alt_path
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                data.append([float(parts[0]), float(parts[1])])
    return data


if __name__ == "__main__":
    data = read_data("assign9fit.txt")
    x_data = [pt[0] for pt in data]
    y_data = [pt[1] for pt in data]

    # Perform cubic polynomial fit (n = 3)
    coeffs, A_matrix, cond_num = poly_fit(data, 3)

    print("Cubic Polynomial Least Square Fit:")
    for i, a in enumerate(coeffs):
        print(f"  a{i} = {a:.6f}")

    # Condition number of the Hankel matrix
    print(f"Condition number of Hankel matrix: {cond_num:.4f}")

    # Polynomial evaluation function
    def f_cubic(x):
        return sum(coeffs[i] * (x ** i) for i in range(len(coeffs)))

    # Smooth curve for plotting
    x_min, x_max = min(x_data), max(x_data)
    n_pts = 300
    x_fit = [x_min + i * (x_max - x_min) / (n_pts - 1) for i in range(n_pts)]
    y_fit = [f_cubic(xi) for xi in x_fit]

    # Plotting
    plt.figure(figsize=(8, 5))
    plt.scatter(x_data, y_data, color='red', s=40, zorder=5, label="Experimental data")
    plt.plot(x_fit, y_fit, 'b-', linewidth=2,
             label=f"Cubic fit: $y = {coeffs[0]:.3f} + {coeffs[1]:.3f}x {coeffs[2]:.3f}x^2 + {coeffs[3]:.3f}x^3$")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Cubic Least Square Fit ($y = a_0 + a_1x + a_2x^2 + a_3x^3$)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("q3_plot.png", dpi=150)
    plt.show()
