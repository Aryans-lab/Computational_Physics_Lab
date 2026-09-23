# Determine the value of π using throwing method by choosing a quarter circle of unit radius in the first quadrant. 
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import matplotlib.pyplot as plt
from mylib import myrand
import math

def estimate_pi(n):
    x_random = [myrand() for _ in range(n)]
    y_random = [myrand() for _ in range(n)]

    inside_circle = [1 for x, y in zip(x_random, y_random) if x**2 + y**2 <= 1]

    pi_estimate = (sum(inside_circle) / n) * 4
    return pi_estimate

# Generate π estimates for N from 20 to 5000
N_values = list(range(20, 5001))
pi_estimates = [estimate_pi(n) for n in N_values]
residuals = [p - math.pi for p in pi_estimates]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

ax1.scatter(N_values, pi_estimates, s=1)
ax1.axhline(math.pi, label='True π')
ax1.set_xlabel('Number of throws (N)')
ax1.set_ylabel('Estimated π')
ax1.set_title('Estimation of π')
ax1.legend()

ax2.scatter(N_values, residuals, s=1)
ax2.set_xlabel('Number of throws (N)')
ax2.set_ylabel('Residual (estimated π − true π)')
ax2.set_title('Residual plot')

plt.show()

#Final value of pi estimated at N=5000
final_pi_estimate = estimate_pi(5000)
print(f"Final estimated value of π at N=5000: {final_pi_estimate}")

# End of code