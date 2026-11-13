"""
fftlib.py -- Radix-2 Fast Fourier Transform, written from scratch.

Why this file exists
--------------------
The project needs three things that all rest on the same transform:

  1. a spectral solver for the 1-D Navier-Stokes (Burgers) equation,
  2. the energy spectrum E(k) of a shock train,
  3. the structure functions S_p(r) of Burgers turbulence.

The course rule is that no black-box numerical routine is used, so `np.fft`
is not called anywhere. This module implements the Cooley-Tukey radix-2
algorithm directly, plus a naive O(N^2) DFT used only to verify it.

Conventions
-----------
Forward:   U_k = (1/N) * sum_j u_j exp(-2*pi*i*j*k/N)        (normalised)
Inverse:   u_j = sum_k U_k exp(+2*pi*i*j*k/N)

With this normalisation Parseval reads  sum_j |u_j|^2 = N * sum_k |U_k|^2,
which is the form used in diagnostics.energy_spectrum.

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import numpy as np


# ----------------------------------------------------------------------
# bit reversal
# ----------------------------------------------------------------------
def bit_reverse_indices(n):
    """Return the bit-reversal permutation of range(n); n must be a power of 2."""
    if n & (n - 1) != 0:
        raise ValueError("fft length must be a power of 2")
    bits = n.bit_length() - 1
    idx = np.arange(n, dtype=np.int64)
    rev = np.zeros(n, dtype=np.int64)
    for _ in range(bits):
        rev = (rev << 1) | (idx & 1)
        idx >>= 1
    return rev


def fft(a):
    """In-place-free radix-2 forward FFT, normalised by 1/N.

    Parameters
    ----------
    a : array_like, length must be a power of 2

    Returns
    -------
    np.ndarray (complex128), same length as a
    """
    a = np.asarray(a, dtype=np.complex128)
    n = a.size
    out = a[bit_reverse_indices(n)].copy()

    m = 2
    while m <= n:
        half = m // 2
        # twiddle factors for this stage
        w = np.exp(-2j * np.pi * np.arange(half) / m)
        block = out.reshape(-1, m)
        even = block[:, :half].copy()
        odd = block[:, half:].copy() * w
        block[:, :half] = even + odd
        block[:, half:] = even - odd
        m *= 2

    return out / n


def ifft(a):
    """Inverse FFT, consistent with fft() being normalised by 1/N."""
    a = np.asarray(a, dtype=np.complex128)
    n = a.size
    return n * np.conjugate(fft(np.conjugate(a)))


def rfft(a):
    """Spectrum of a real signal: modes k = 0 ... N/2 (inclusive)."""
    n = np.asarray(a).size
    return fft(a)[: n // 2 + 1]


def dft(a):
    """Naive O(N^2) discrete Fourier transform. Used only to test fft()."""
    a = np.asarray(a, dtype=np.complex128)
    n = a.size
    j = np.arange(n)
    k = j.reshape(-1, 1)
    return (a.reshape(1, -1) * np.exp(-2j * np.pi * j * k / n)).sum(axis=1) / n


# ----------------------------------------------------------------------
# spectral derivatives on a periodic grid
# ----------------------------------------------------------------------
class PeriodicGrid:
    """Evenly spaced periodic grid on [0, L) with the integer wavenumbers.

    x_j = j*L/N,  j = 0 .. N-1
    k   = 0, 1, ..., N/2-1, -N/2, ..., -1   (standard FFT ordering)
    """

    def __init__(self, n, length=2.0 * np.pi):
        if n & (n - 1) != 0:
            raise ValueError("n must be a power of 2")
        self.n = n
        self.L = length
        self.dx = length / n
        self.x = np.arange(n) * self.dx
        self.k = np.fft.fftfreq(n, d=1.0 / n).astype(np.int64) * (2.0 * np.pi / length)

    def derivative(self, u):
        """Spectral first derivative du/dx of a periodic field."""
        uh = fft(u)
        return np.real(ifft(1j * self.k * uh))

    def wavenumbers(self):
        """Non-negative integer wavenumbers for spectra (0..N/2)."""
        kmax = self.n // 2
        return np.arange(kmax + 1)


def dealias_mask(n, fraction=2.0 / 3.0):
    """Orszag 2/3-rule mask: zero the top third of the Fourier modes."""
    k = np.fft.fftfreq(n, d=1.0 / n).astype(np.int64)
    return np.abs(k) <= fraction * (n // 2)


def fourier_interpolate(u, n_out):
    """Band-limited (spectral) interpolation of a periodic field onto n_out points.

    The Fourier coefficients of the n-point grid are zero-padded to n_out
    modes and transformed back, so the interpolant is the unique band-limited
    function through the samples. This matters when the result is used as a
    *reference*: linear interpolation of a coarse field is only O(dx^2)
    accurate and would contaminate a spectral benchmark.
    """
    u = np.asarray(u, dtype=float)
    n = u.size
    if n_out == n:
        return u.copy()
    if (n & (n - 1)) or (n_out & (n_out - 1)):
        raise ValueError("both lengths must be powers of 2")
    uh = fft(u)
    out = np.zeros(n_out, dtype=np.complex128)
    half = n // 2
    if n_out > n:
        out[:half] = uh[:half]
        out[n_out - half + 1:] = uh[half + 1:]
        out[half] = 0.5 * uh[half]          # split the Nyquist mode
        out[n_out - half] = 0.5 * uh[half]
    else:
        out[:n_out // 2] = uh[:n_out // 2]
        out[n_out // 2 + 1:] = uh[n_out // 2 + 1:]
    return np.real(ifft(out))
