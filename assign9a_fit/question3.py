# Use Lagrange's interpolation formula to find the best estimate for y(x = 6.7) from the data table below.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

from mylib import poly_fit, poly_eval, lagrange_interpolation
import matplotlib.pyplot as plt

if __name__ == "__main__":
    X = [[2,10],[3,15],[5,25],[8,40],[12,60]]
    x_data = [row[0] for row in X]
    y_data = [row[1] for row in X]

    y_est, coeffs = lagrange_interpolation(X, 6.7)
    print(f"y(6.7) = {y_est:.6f}")
    print("Polynomial coefficients (c0 + c1*x + ...):")
    for i, c in enumerate(coeffs):
        print(f"  c{i} = {c:.6e}")

    x_min, x_max = min(x_data) - 0.5, max(x_data) + 0.5
    x_curve = [x_min + i * (x_max - x_min) / 299 for i in range(300)]
    y_curve = [lagrange_interpolation(X, xi)[0] for xi in x_curve]

    plt.figure(figsize=(8, 5))
    plt.scatter(x_data, y_data, color='red', zorder=5, label="Data points")
    plt.plot(x_curve, y_curve, 'b-', label="Lagrange polynomial")
    plt.scatter([6.7], [y_est], color='green', marker='*', s=200, zorder=6,
                label=f"y(6.7) = {y_est:.4f}")
    coeff_str = "\n".join(f"$c_{i}$ = {c:.3e}" for i, c in enumerate(coeffs))
    plt.text(0.05, 0.95, coeff_str, transform=plt.gca().transAxes,
             verticalalignment='top', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Lagrange Interpolation with Polynomial Coefficients")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("question3_plot.png", dpi=150)
    plt.show()

# End of Code
