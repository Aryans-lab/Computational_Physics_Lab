# Use Monte Carlo integration scheme (using your own pRNG) to estimate the integral accurate
# up to 4 places in decimal. Integral of -1 to 1 sin^2x. PLot the Integral value vs N and deviation vs N
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import math
import matplotlib.pyplot as plt
from mylib import monte_carlo

if __name__ == "__main__":
    # Function to integrate
    def f(x):
        return math.sin(x) ** 2

    # Exact value:
    # Integral from -1 to 1 of sin^2(x) dx
    exact_value = 1 - math.sin(2) / 2

    N_values = [100 * i for i in range(11, 501)]

    integral_values = []
    deviations = []
    sigma_values = []

    for i, N in enumerate(N_values):
        # Different seed for every Monte Carlo run
        seed = i + 1
        integral, sigma_f, sigma_integral = monte_carlo(
            f, -1.0, 1.0, N, seed
        )

        integral_values.append(integral)

        # Absolute deviation from exact value
        deviation = abs(integral - exact_value)
        deviations.append(deviation)
        sigma_values.append(sigma_integral)

        print(
            f"N = {N:5d}, "
            f"Seed = {seed:3d}, "
            f"Integral = {integral:.8f}, "
            f"Deviation = {deviation:.8f}, "
            f"Sigma = {sigma_integral:.8f}"
        )

    # -------------------------------------------------------
    # Find the first N for which accuracy is within 10^-4
    # -------------------------------------------------------
    accuracy = 1e-4

    first_N = None
    first_integral = None
    first_deviation = None

    for i in range(len(N_values)):
        if deviations[i] < accuracy:
            first_N = N_values[i]
            first_integral = integral_values[i]
            first_deviation = deviations[i]
            break

    print()
    print("==================================================")
    print("RESULT")
    print("==================================================")
    print(f"Exact value = {exact_value:.6f}")

    if first_N is not None:
        print(f"Required accuracy = {accuracy}")
        print(f"First N satisfying the accuracy = {first_N}")
        print(f"Integral = {first_integral:.6f}")
        print(f"Deviation = {first_deviation:.6f}")
    else:
        print("Accuracy of 10^-4 was not reached")
        print("within the tested range of N = 100 to 10000.")

    # -------------------------------------------------------
    # Integral value vs N
    # -------------------------------------------------------
    plt.figure(figsize=(8, 6))
    plt.plot(
        N_values, integral_values, marker="o", markersize=4, label="Monte Carlo"
    )
    plt.axhline(
        y=exact_value, linestyle="--", label=f"Exact = {exact_value:.8f}"
    )
    plt.xlabel("N")
    plt.ylabel("Integral Value")
    plt.title("Monte Carlo Integration: Integral Value vs N")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("integral_vs_N.png", dpi=600)
    plt.show()

    # -------------------------------------------------------
    # Deviation vs N
    # -------------------------------------------------------
    plt.figure(figsize=(8, 6))
    plt.plot(
        N_values, deviations, marker="o", markersize=4, label="Absolute Deviation"
    )
    plt.axhline(
        y=1e-4, linestyle="--", label=r"Required accuracy = $10^{-4}$"
    )
    plt.xlabel("N")
    plt.ylabel("Absolute Deviation")
    plt.title("Monte Carlo Integration: Deviation vs N")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("deviation_vs_N.png", dpi=600)
    plt.show()


# End of code