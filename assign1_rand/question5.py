# Generate pRNG having exponential distribution of the form exp(−x) from pRNG
# having uniform distribution in [0, 1). Generate at least 5,000 random numbers
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import matplotlib.pyplot as plt
from math import log
from lcg import myrand

if __name__ == '__main__':
    x_rn = []
    y_rn = []
    for i in range(5000):
        x = myrand()
        x_rn.append(x)
        y_rn.append(-log(x))

    plt.hist(y_rn, bins=40, color='skyblue', edgecolor='black')
    plt.xlabel('y')
    plt.ylabel('Counts')
    plt.title('Exponential distribution from LCG pRNG (N=5000)')
    plt.tight_layout()
    plt.savefig('exponential_distribution.png', dpi=120)
    plt.show()
