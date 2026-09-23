# lcg.py — pRNG library (LCG-based)
# Parameters (gcc standard): a=1103515245, c=12345, m=32768
# Author: Aryan Bandyopadhyay, Roll number: 2411014

def myrand(seed=0):
    if not hasattr(myrand, 'x0'):
        myrand.x0 = 0
    a = 1103515245
    c = 12345
    m = 32768
    if seed:
        myrand.x0 = seed
    x1 = (a * myrand.x0 + c) % m
    rn = x1 / m
    myrand.x0 = x1
    return rn

