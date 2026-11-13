# Runge-Kutta 4th Order (RK4) routine for 1st and 2nd order ODEs
# Name: Aryan Bandyopadhyay, Roll number: 2411014

def rk4(f, a, b, y0, h=None, n=None):
    if h is not None:
        n = int((b - a) / h)
    elif n is not None:
        h = (b - a) / n
    else:
        raise ValueError("Either step size h or number of steps n must be provided.")

    x = [a + i * h for i in range(n + 1)]
    is_system = isinstance(y0, (list, tuple))

    if is_system:
        m = len(y0)
        y = [[float(v) for v in y0]]
        for i in range(n):
            xi = x[i]
            yi = y[i]

            # Forward Euler intermediate steps
            k1 = [float(v) for v in f(xi, yi)]
            y_mid1 = [yi[j] + 0.5 * h * k1[j] for j in range(m)]

            k2 = [float(v) for v in f(xi + 0.5 * h, y_mid1)]
            y_mid2 = [yi[j] + 0.5 * h * k2[j] for j in range(m)]

            k3 = [float(v) for v in f(xi + 0.5 * h, y_mid2)]
            y_end = [yi[j] + h * k3[j] for j in range(m)]

            k4 = [float(v) for v in f(xi + h, y_end)]

            y_next = [yi[j] + (h / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) for j in range(m)]
            y.append(y_next)

        return x, y
    else:
        y = [float(y0)]
        for i in range(n):
            xi = x[i]
            yi = y[i]

            # Forward Euler intermediate steps
            k1 = float(f(xi, yi))
            y_mid1 = yi + 0.5 * h * k1

            k2 = float(f(xi + 0.5 * h, y_mid1))
            y_mid2 = yi + 0.5 * h * k2

            k3 = float(f(xi + 0.5 * h, y_mid2))
            y_end = yi + h * k3

            k4 = float(f(xi + h, y_end))

            # RK4 combination
            y_next = yi + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            y.append(y_next)

        return x, y

