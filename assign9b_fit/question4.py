# Assignment 9B: Question 4
# Try fitting the data set with a modified Chebyshev basis functions
# y = f(x) = a0*phi_0(x) + a1*phi_1(x) + a2*phi_2(x) + a3*phi_3(x)
# where, phi_0(x) = 1, phi_1(x) = 2x - 1, phi_2(x) = 8x^2 - 8x + 1, phi_3(x) = 32x^3 - 48x^2 + 18x - 1
# Can you explain the difference by comparing the relative merit of using different functions
# calculating the condition number? (You may use online condition number or Python's in-built condition number calculator.)
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import os
import matplotlib.pyplot as plt
from mylib import poly_fit, chebyshev_fit, modified_chebyshev_basis, gauss_jordan_inverse, frobenius_norm


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

    # Standard cubic polynomial fit (Q3)
    coeffs_poly, A_poly, cond_poly = poly_fit(data, 3)

    # Modified Chebyshev fit (Q4)
    coeffs_cheb, A_cheb = chebyshev_fit(data, 3)
    A_cheb_inv = gauss_jordan_inverse(A_cheb)
    cond_cheb = frobenius_norm(A_cheb) * frobenius_norm(A_cheb_inv)

    print("1. Modified Chebyshev Basis Coefficients:")
    for i, a in enumerate(coeffs_cheb):
        print(f"  a{i} = {a:.6f}")

    print("\n2. Standard Monomial Coefficients:")
    for i, a in enumerate(coeffs_poly):
        print(f"  a{i} = {a:.6f}")

    print("\n3. Condition Numbers Comparison:")
    print(f"  Condition number (Standard Polynomial Hankel Matrix): {cond_poly:.4f}")
    print(f"  Condition number (Modified Chebyshev Matrix): {cond_cheb:.4f}")

    # Evaluation function for Chebyshev fit
    def f_cheb(x):
        return sum(coeffs_cheb[m] * modified_chebyshev_basis(m, x) for m in range(len(coeffs_cheb)))

    def f_poly(x):
        return sum(coeffs_poly[m] * (x ** m) for m in range(len(coeffs_poly)))

    # curve
    x_min, x_max = min(x_data), max(x_data)
    n_pts = 300
    x_fit = [x_min + i * (x_max - x_min) / (n_pts - 1) for i in range(n_pts)]
    y_fit_cheb = [f_cheb(xi) for xi in x_fit]
    y_fit_poly = [f_poly(xi) for xi in x_fit]

    # Plot comparison
    plt.figure(figsize=(9, 5.5))
    plt.scatter(x_data, y_data, color='red', s=45, zorder=5, label="Experimental data")
    plt.plot(x_fit, y_fit_cheb, 'b-', linewidth=2.5,
             label=f"Modified Chebyshev fit (cond = {cond_cheb:.2f})")
    plt.plot(x_fit, y_fit_poly, 'k--', linewidth=1.5,
             label=f"Standard Cubic fit (cond = {cond_poly:.2f})")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Comparison: Standard Polynomial vs Modified Chebyshev Fit")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("q4_plot.png", dpi=600)
    plt.show()

# End of Code