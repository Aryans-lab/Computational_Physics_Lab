"""
question2.py
------------
Problem : Calculate the sum of N = 15 terms of a Geometric Progression (GP)
          and a Harmonic Progression (HP) for:
            common ratio  r  = 0.5
            common difference d = 1.5
            first term    t0 = 1.25
          Analytical closed-form formulae are NOT used; terms are summed
          iteratively.
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER
"""

# -- Parameters
N  = 15     # number of terms
t0 = 1.25   # first term (common to both series)
r  = 0.5    # common ratio (GP)
d  = 1.5    # common difference (underlying AP for HP)


def sum_gp(n: int, first: float, ratio: float) -> float:
    """
    Sum the first n terms of a GP iteratively.
    GP: t_k = first * ratio^(k-1)
    """
    total = 0.0
    term  = first
    for _ in range(n):
        total += term
        term  *= ratio
    return total


def sum_hp(n: int, first: float, diff: float) -> float:
    """
    Sum the first n terms of a HP iteratively.
    HP: t_k = 1 / a_k  where a_k is the k-th term of the underlying AP.
        a_k = first + (k-1)*diff
    """
    total   = 0.0
    ap_term = first          # a_k starts at t0 (the AP term)
    for _ in range(n):
        total   += 1.0 / ap_term
        ap_term += diff
    return total


# -- Driver
print("Sum of GP series:", sum_gp(N, t0, r))
print("Sum of HP series:", sum_hp(N, t0, d))

# -- Output
# Sum of GP series: 2.4999237060546875
# Sum of HP series: 2.4139570733659186
