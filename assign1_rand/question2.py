# Problem: Write your own LCG random generator with the following set of parameters a = 1103515245, c = 12345, m = 32768 
# and again check for correlation by plotting for k = 5. 
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from lcg import myrand
import matplotlib.pyplot as plt

def lcg(n, x0=1.2, a=1103515245, c=12345, m=32768):
    myrand(int(x0))
    lcg_random = [int(myrand() * m) for _ in range(n)]
    return lcg_random

def plot_correlation_lcg(x0, n, k):
    lcg_random=lcg(x0=x0, n=n)
    x_vals=[]
    y_vals=[]
    for i in range(len(lcg_random)-k):
        x_vals.append(lcg_random[i])
        y_vals.append(lcg_random[i+k])
    fig, ax = plt.subplots()
    ax.scatter(x_vals, y_vals, s=2)
    plt.xlabel('xi')
    plt.ylabel('xi+k')
    plt.title(f'Correlation Plot (k={k}) with (n={n}) of LCG Random Numbers')
    plt.savefig('correlation_Question2.png', dpi=120)
    plt.show()

plot_correlation_lcg(x0=1.2, n=2000, k=5)

# End of code