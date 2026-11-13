"""
exact.py -- analytic benchmarks for the Burgers equation.

Three exact results are used to keep the numerical work honest:

1. The implicit solution along characteristics,
       u = u0(xi),   x = xi + u0(xi) t ,
   valid until characteristics first cross, i.e. until
       t_s = -1 / min_x u0'(x) .
   For u0 = sin x this gives t_s = 1, and the shock appears at x = pi.

2. The Cole-Hopf transformation,
       u = -2 nu (ln phi)_x ,     phi_t = nu phi_xx ,
   which turns the nonlinear equation into the heat equation and therefore
   gives an exact reference solution for any smooth u0, valid for all t.

3. The exact steady shock (a travelling wave with zero speed),
       u(x) = -U tanh( U (x - x0) / (2 nu) ),
   whose 10-90 % width is exactly 2 nu / U in the argument scale, and whose
   amplitude-to-slope ratio is 4 nu / U.

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import numpy as np

from fftlib import fourier_interpolate, fft, ifft

TWO_PI = 2.0 * np.pi


# ----------------------------------------------------------------------
# 1. characteristics
# ----------------------------------------------------------------------
def characteristic_shock_time(u0_func, x):
    """t_s = -1/min u0'(x), the first time characteristics cross."""
    xs = np.linspace(x[0], x[-1], 20001)[:-1]
    h = 1e-6
    du = (u0_func(xs + h) - u0_func(xs - h)) / (2.0 * h)
    slope = du.min()
    if slope >= 0.0:
        return np.inf
    return -1.0 / slope


def characteristic_solution(xi, t, u0_func=np.sin):
    """Exact implicit solution. Returns (x, u) points, no inversion needed.

    Plot u against x parametrically by xi. Valid only for t < t_s.
    """
    u = u0_func(xi)
    return xi + u * t, u


# ----------------------------------------------------------------------
# 2. Cole-Hopf
# ----------------------------------------------------------------------
def cole_hopf(u0_grid, x, t, nu, n_quad=4096, n_images=8):
    """Exact solution of u_t + u u_x = nu u_xx on a periodic domain.

    phi(x,t) = integral K(x-y,t) phi0(y) dy,  K = periodic heat kernel,
    phi0(y)  = exp(-(1/2nu) * integral_0^y u0(s) ds).

    The integral is done with the trapezoidal rule on a refined periodic grid
    (n_quad points). Because phi0 may span hundreds of orders of magnitude,
    it is rescaled by its maximum; the rescaling cancels in u = -2 nu d ln phi.
    """
    n = len(x)
    L = x[-1] - x[0] + (x[1] - x[0])
    y = np.linspace(0.0, L, n_quad, endpoint=False)
    dy = y[1] - y[0]

    # U0(y) = integral_0^y u0 ds. The initial data is interpolated onto the
    # quadrature grid *spectrally*, not linearly: a linear interpolant of a
    # coarse grid would inject an O(dx^2) error into the reference itself.
    if len(u0_grid) != n_quad and (len(u0_grid) & (len(u0_grid) - 1)) == 0:
        u0_on_y = fourier_interpolate(u0_grid, n_quad)
    else:
        u0_on_y = np.interp(y, x, u0_grid, period=L)
    # spectral integration: U0_hat = u_hat / (i k), exact for band-limited data
    kk = np.fft.fftfreq(n_quad, d=1.0 / n_quad).astype(np.int64)
    uh = fft(u0_on_y)
    U0_hat = np.zeros_like(uh)
    nz = kk != 0
    U0_hat[nz] = uh[nz] / (1j * kk[nz])
    U0 = np.real(ifft(U0_hat))
    U0 -= U0.mean()

    phi0 = np.exp(-(U0 - U0.min()) / (2.0 * nu))
    phi0 /= phi0.max()

    # heat kernel, summed over periodic images
    X, Y = np.meshgrid(x, y, indexing="ij")           # shape (n, n_quad)
    phi = np.zeros(n)
    phix = np.zeros(n)
    pref = 1.0 / np.sqrt(4.0 * np.pi * nu * t)
    for m in range(-n_images, n_images + 1):
        d = X - Y + m * L
        g = np.exp(-d * d / (4.0 * nu * t))
        G = pref * g
        phi += G @ phi0 * dy
        phix += (-d / (2.0 * nu * t)) * G @ phi0 * dy

    return -2.0 * nu * phix / phi


# ----------------------------------------------------------------------
# 3. exact steady shock
# ----------------------------------------------------------------------
def tanh_shock(x, U=1.0, nu=0.05, x0=0.0):
    """Exact steady solution u = -U tanh(U(x-x0)/(2 nu))."""
    return -U * np.tanh(U * (x - x0) / (2.0 * nu))


def tanh_width(nu, U=1.0):
    """Amplitude-to-maximum-slope width of the exact tanh shock: 4 nu / U."""
    return 4.0 * nu / U
