# 4. Solve the Rossler system for a = 0.2, b = 0.2, c = 5.7 and show a 3-dim (x, y, z) plot.
# dx/dt = -y - z, dy/dt = x + ay, dz/dt = b + z(x - c)
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import matplotlib.pyplot as plt
from mylib import rk4

# Parameters
a = 0.2
b = 0.2
c = 5.7

def rossler(t, r):
    x, y, z = r
    dxdt = -y - z
    dydt = x + a * y
    dzdt = b + z * (x - c)
    return [dxdt, dydt, dzdt]


if __name__ == '__main__':
    # Initial values and t range
    t0 = 0.0
    tf = 100.0
    h = 0.01
    init_vals = [0.1, 0.0, 0.0]

    t, sol = rk4(rossler, t0, tf, init_vals, h=h)

    x = [pt[0] for pt in sol]
    y = [pt[1] for pt in sol]
    z = [pt[2] for pt in sol]

    # 3D plot of the Rossler attractor
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, color='purple', lw=0.6)
    ax.set_title("Rossler Attractor (a=0.2, b=0.2, c=5.7)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    plt.tight_layout()
    plt.savefig("q4_rossler.png")
    plt.show()