# Code up library functions Bisection and Regula Falsi.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

# Bisection method for finding roots of a function f

def bisection(f, a, b, accuracy=1e-8, max_iterations=100):
    if f(a) * f(b) >= 0:
        while f(a) * f(b) >= 0:
            m=0.5
            if abs(f(a)) < abs(f(b)):
                a = a - m*(b-a)
            if abs(f(b)) < abs(f(a)):
                b = b + m*(b-a)
                m +=0.1                

    for i in range(max_iterations):
        c = (a + b) / 2
        if abs(f(c)) < accuracy and abs(b - a) < accuracy:
            print(f"The root is approximately {c}")
            print(f"The number of iterations taken is {i + 1}")
            return c
        if f(c) * f(a) < 0:
            b = c
        else:
            a = c
    return c, max_iterations

# Regula Falsi method

def regula_falsi(f, a, b, accuracy=1e-8, max_iterations=100):
    if f(a) * f(b) >= 0:
        while f(a) * f(b) >= 0:
            m=0.5
            if abs(f(a)) < abs(f(b)):
                a = a - m*(b-a)
            if abs(f(b)) < abs(f(a)):
                b = b + m*(b-a)
                m +=0.2

    for i in range(max_iterations):
        c = (a * f(b) - b * f(a)) / (f(b) - f(a))
        j=0
        if abs(f(c)) < accuracy and abs(b - a) < accuracy:
            print(f"The root is approximately {c}")
            print(f"The number of iterations taken is {i + 1}")
            return c
        if f(c) * f(a) < 0:
            b = c
            j+=1
        else:
            a = c
            j+=1

    return c, max_iterations 

# End of Code