# Use Forward (explicit) Euler and Backward (implicit) Euler and Predictor-Corrector
# methods to solve and the following differential equations and compare with analytical 
# form over range x ∈ [0, 2] and x ∈ [0, π/5] for the first and the second differential equation by plotting . 
# Use multiple step sizes 0.02, 0.1, 0.3, 0.5 in both the cases.
# 1. dy/dx=y-x^2 with y(0)=0 and y(x)=x^2 + 2x + 2 - 2e^x
# 2. dy/dx=(x+y)^2 with y(0)=1 tan^(-1)(x+y)= x + pi/4
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from mylib import newton_raphson
import numpy as np
import matplotlib.pyplot as plt

def euler(f, a, b, y0, h, method='forward'):
    if method=='forward':
        n = int((b - a) / h)
        x = [a + i * h for i in range(n + 1)]
        y = [y0]
        for i in range(n):
            y.append(y[i] + h * f(x[i], y[i]))
        return x, y
    
    if method=='backward':
        n = int((b - a) / h)
        x = [a + i * h for i in range(n + 1)]
        y = [y0]
        for i in range(n):
            # Use Newton-Raphson to solve for y[i+1]
            func = lambda y_next: y_next - y[i] - h * f(x[i + 1], y_next)
            guess = y[i] + h * f(x[i], y[i])
            res = newton_raphson(func, x0=guess)
            if isinstance(res, tuple):
                res = res[0]
            elif res is None:
                res = guess
            y.append(res)
        return x, y

def predictor_corrector(f, a, b, y0, h):
    n = int((b - a) / h)
    x = [a + i * h for i in range(n + 1)]
    y = [y0]
    for i in range(n):
        # Predictor step (Forward Euler)
        y_predictor = y[i] + h * f(x[i], y[i])
        # Corrector step 
        y_corrector = y[i] + (h / 2) * (f(x[i], y[i]) + f(x[i + 1], y_predictor))
        y.append(y_corrector)
    return x, y


if __name__=='__main__':

    # Define the differential equations
    def f1(x, y):
        return y - x**2

    def f2(x, y):
        return (x + y)**2

    # Derivative of functions
    def df1_dy(x, y):
        return 1.0

    def df2_dy(x, y):
        return 2.0 * (x + y)

    # Analytical solutions
    def analytical1(x):
        return x**2 + 2*x + 2 - 2*np.exp(x)

    def analytical2(x):
        return np.tan(x + np.pi/4) - x

    # Parameters
    a1, b1 = 0, 2
    a2, b2 = 0, np.pi/5
    y0_1 = 0
    y0_2 = 1
    step_sizes = [0.02, 0.1, 0.3, 0.5]

    # Solve and plot for the first differential equation with forward Euler method
    plt.figure(figsize=(12, 6))
    for h in step_sizes:
        x1, y1 = euler(f1, a1, b1, y0_1, h)
        plt.plot(x1, y1, label=f'Forward Euler h={h}')
    
    x_analytical1 = np.linspace(a1, b1, 100)
    plt.plot(x_analytical1, analytical1(x_analytical1), 'k--', label='Analytical Solution')
    plt.title('Differential Equation 1: dy/dx=y-x^2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()

    # Solve and plot for the second differential equation with foreward Euler method
    plt.figure(figsize=(12, 6))
    for h in step_sizes:
        x2, y2 = euler(f2, a2, b2, y0_2, h)
        plt.plot(x2, y2, label=f'Forward Euler h={h}')
    
    x_analytical2 = np.linspace(a2, b2, 100)
    plt.plot(x_analytical2, analytical2(x_analytical2), 'k--', label='Analytical Solution')
    plt.title('Differential Equation 2: dy/dx=(x+y)^2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()

    # Solve and plot for the first differential equation with backward Euler method
    plt.figure(figsize=(12, 6))
    for h in step_sizes:
        x1, y1 = euler(f1, a1, b1, y0_1, h, method='backward')
        plt.plot(x1, y1, label=f'Backward Euler h={h}')
    
    x_analytical1 = np.linspace(a1, b1, 100)
    plt.plot(x_analytical1, analytical1(x_analytical1), 'k--', label='Analytical Solution')
    plt.title('Differential Equation 1: dy/dx=y-x^2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()

    # Solve and plot for the second differential equation with backward Euler method
    plt.figure(figsize=(12, 6))
    for h in step_sizes:
        x2, y2 = euler(f2, a2, b2, y0_2, h, method='backward')
        plt.plot(x2, y2, label=f'Backward Euler h={h}')
    
    x_analytical2 = np.linspace(a2, b2, 100)
    plt.plot(x_analytical2, analytical2(x_analytical2), 'k--', label='Analytical Solution')
    plt.title('Differential Equation 2: dy/dx=(x+y)^2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()

    # Solve and plot for the first differential equation with predictor-corrector method
    plt.figure(figsize=(12, 6))
    for h in step_sizes:
        x1, y1 = predictor_corrector(f1, a1, b1, y0_1, h)
        plt.plot(x1, y1, label=f'Predictor-Corrector h={h}')
    
    x_analytical1 = np.linspace(a1, b1, 100)
    plt.plot(x_analytical1, analytical1(x_analytical1), 'k--', label='Analytical Solution')
    plt.title('Differential Equation 1: dy/dx=y-x^2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()

    # Solve and plot for the second differential equation with predictor-corrector method
    plt.figure(figsize=(12, 6))
    for h in step_sizes:
        x2, y2 = predictor_corrector(f2, a2, b2, y0_2, h)
        plt.plot(x2, y2, label=f'Predictor-Corrector h={h}')
    
    x_analytical2 = np.linspace(a2, b2, 100)
    plt.plot(x_analytical2, analytical2(x_analytical2), 'k--', label='Analytical Solution')
    plt.title('Differential Equation 2: dy/dx=(x+y)^2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()
