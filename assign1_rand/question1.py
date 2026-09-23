# Problem: Use the iterative equation xi+1 = c xi (1 − xi) to generate 1,000 random numbers.
# Show the correlation among them by plotting xi vs xi+k. Use seed x0 = 0.1 but
# choose your own three (3) different c. Try various k, say 3, 5 and 10.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import matplotlib.pyplot as plt

random_numbers= []

#Function for the random number generation using the it. equation
def random_number_iterative(x0, c, n):
    global random_numbers
    random_numbers = []   # reset for each new c
    xi=x0
    for i in range(n):
        xi=c*xi*(1-xi)
        random_numbers.append(xi)

#Plot for corrolation among the random numbers generated (for different k)
def plot_correlation(x0, c, n, k, ax):
    #fill the list with RNs
    random_number_iterative(x0, c, n)
    x_vals=[]
    y_vals=[]
    for i in range(len(random_numbers)-k):
        x_vals.append(random_numbers[i])
        y_vals.append(random_numbers[i+k])
    ax.scatter(x_vals, y_vals, s=2)
    ax.set_xlabel('xi')
    ax.set_ylabel('xi+k')
    ax.set_title(f'c={c}, k={k}')


c_values = [3.85, 3.95, 4.0]
k_values = [3, 5, 10]

fig, axes = plt.subplots(3, 3, figsize=(10, 9))
fig.suptitle('Correlation Plot (N=1000)')

for i, c in enumerate(c_values):
    for j, k in enumerate(k_values):
        plot_correlation(0.1, c, 1000, k, axes[i][j])

plt.tight_layout()
plt.show()

# End of the code