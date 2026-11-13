# Fit the data given in the table below with power law y = ax^b and exponential y = ae^(-bx) models.
# Based on Pearson's r^2, determine which model gives better fit.
# Data: (2.5,13.0),(3.5,11.0),(5.0,8.5),(6.0,8.2),(7.5,7.0),(10.0,6.2),(12.5,5.2),(15.0,4.8),(17.5,4.6),(20.5,4.3)
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

from mylib import linear_regression
import matplotlib.pyplot as plt
import math

if __name__ == "__main__":
    # Dataset
    data = [(2.5,13.0),(3.5,11.0),(5.0,8.5),(6.0,8.2),(7.5,7.0),
            (10.0,6.2),(12.5,5.2),(15.0,4.8),(17.5,4.6),(20.5,4.3)]

    x_data = [pt[0] for pt in data]
    y_data = [pt[1] for pt in data]

    X_pow = [[math.log(xi), math.log(yi)] for xi, yi in zip(x_data, y_data)]
    b_pow, lna_pow, r2_pow = linear_regression(X_pow)
    a_pow = math.exp(lna_pow)

    print("MODEL 1: Power Law  y = a * x^b")
    print(f"  a = {a_pow:.6f}")
    print(f"  b = {b_pow:.6f}")
    print(f"  R² = {r2_pow:.6f}")

    X_exp = [[xi, math.log(yi)] for xi, yi in zip(x_data, y_data)]
    slope_exp, lna_exp, r2_exp = linear_regression(X_exp)
    a_exp = math.exp(lna_exp)
    b_exp = -slope_exp        

    print()
    print("MODEL 2: Exponential  y = a * e^(-b*x)")
    print(f"  a = {a_exp:.6f}")
    print(f"  b = {b_exp:.6f}")
    print(f"  R² = {r2_exp:.6f}")

    print()
    if r2_pow >= r2_exp:
        print(f"=> Power law fits better  (R²_pow = {r2_pow:.4f}  >  R²_exp = {r2_exp:.4f})")
    else:
        print(f"=> Exponential fits better  (R²_exp = {r2_exp:.4f}  >  R²_pow = {r2_pow:.4f})")

    x_min = min(x_data) * 0.9
    x_max = max(x_data) * 1.05
    n_pts = 400
    h = (x_max - x_min) / (n_pts - 1)
    x_curve = [x_min + i * h for i in range(n_pts)]

    y_pow = [a_pow * xi ** b_pow               for xi in x_curve]
    y_exp = [a_exp * math.exp(-b_exp * xi)     for xi in x_curve]

    plt.figure(figsize=(9, 6))
    plt.scatter(x_data, y_data, color='red', zorder=10, s=60, label="Data points")
    plt.plot(x_curve, y_pow, 'b-', linewidth=2,
             label=f"Power law: $y = {a_pow:.3f}\\,x^{{{b_pow:.3f}}}$  (R² = {r2_pow:.4f})")
    plt.plot(x_curve, y_exp, 'g--', linewidth=2,
             label=f"Exponential: $y = {a_exp:.3f}\\,e^{{-{b_exp:.4f}x}}$  (R² = {r2_exp:.4f})")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Power Law vs Exponential Fit")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("question4_plot.png", dpi=150)
    plt.show()

# End of Code
