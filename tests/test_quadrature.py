"""Unit tests for the quadrature routines.

Convergence orders are verified empirically (error ratios when N doubles)
and the Gaussian tables are checked against NumPy's orthogonal-polynomial
routines, which independently compute the same nodes and weights.
"""

from __future__ import annotations

import math

import pytest
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss

import mylib

# ---------------------------------------------------------------------------
# Newton-Cotes rules
# ---------------------------------------------------------------------------

def test_midpoint_convergence_order_2():
    f = lambda x: math.exp(x)  # noqa: E731
    exact = math.e - 1.0
    errors = [abs(mylib.midpoint(f, 0.0, 1.0, N) - exact) for N in (4, 8, 16, 32)]
    for old, new in zip(errors, errors[1:]):
        assert 3.0 < old / new < 5.0  # ~4x on doubling N => O(1/N^2)


def test_trapezoidal_convergence_order_2():
    f = lambda x: math.sin(x)  # noqa: E731
    exact = 2.0  # int_0^pi sin x dx = [-cos x]_0^pi = 2
    errors = [abs(mylib.trapezoidal(f, 0.0, math.pi, N) - exact) for N in (4, 8, 16, 32)]
    for old, new in zip(errors, errors[1:]):
        assert 3.0 < old / new < 5.0


def test_simpson_convergence_order_4():
    f = lambda x: math.cos(x)  # noqa: E731
    exact = math.sin(math.pi / 2) + math.sin(-math.pi / 2) * -1  # int cos on [-pi/2, pi/2] = 2
    exact = 2.0
    errors = [abs(mylib.simpson(f, -math.pi / 2, math.pi / 2, N)[0] - exact)
              for N in (4, 8, 16, 32)]
    for old, new in zip(errors, errors[1:]):
        assert 12.0 < old / new < 20.0  # ~16x on doubling N => O(1/N^4)


def test_simpson_exact_on_cubics():
    f = lambda x: x ** 3 + 2 * x ** 2 - x + 5  # noqa: E731
    exact = 0.25 + 2 / 3 - 0.5 + 5  # on [0, 1]
    for N in (2, 4, 100):
        val, _ = mylib.simpson(f, 0.0, 1.0, N)
        assert val == pytest.approx(exact, abs=1e-12)


def test_simpson_rejects_odd_N():
    with pytest.raises(ValueError, match="even"):
        mylib.simpson(math.sin, 0.0, 1.0, 5)


def test_simpson_reports_evaluation_count():
    val, evals = mylib.simpson(math.sin, 0.0, 1.0, 20)
    assert evals == 21


# ---------------------------------------------------------------------------
# Monte Carlo
# ---------------------------------------------------------------------------

def test_monte_carlo_is_seeded_and_reproducible():
    f = lambda x: x * x  # noqa: E731
    a = mylib.monte_carlo(f, 0.0, 1.0, 1000, seed=7)
    b = mylib.monte_carlo(f, 0.0, 1.0, 1000, seed=7)
    c = mylib.monte_carlo(f, 0.0, 1.0, 1000, seed=8)
    assert a == b
    assert a != c


def test_monte_carlo_within_sigma_of_exact():
    f = lambda x: math.sin(x) ** 2  # noqa: E731
    exact = 1.0 - math.sin(2.0) / 2.0
    value, _, sigma = mylib.monte_carlo(f, -1.0, 1.0, 20000, seed=3)
    assert abs(value - exact) < 5.0 * sigma  # a correct estimator: O(1) sigma


def test_monte_carlo_error_scaling():
    # sigma estimate should shrink like 1/sqrt(N)
    f = lambda x: x * x  # noqa: E731
    _, _, s1 = mylib.monte_carlo(f, 0.0, 1.0, 1000, seed=1)
    _, _, s2 = mylib.monte_carlo(f, 0.0, 1.0, 4000, seed=2)
    assert 1.4 < s1 / s2 < 2.1  # ~2x = sqrt(4)


# ---------------------------------------------------------------------------
# Error-bound helper
# ---------------------------------------------------------------------------

def test_error_bound_N_known_values():
    # 1/x on [1,2], M2 = 2: N = ceil(sqrt(2/24e-6)) = 289
    assert mylib.integration_error_bound_N("midpoint", 1, 2, 2.0, 1e-6) == 289
    # M4 = 24: N = ceil((24/180e-6)^(1/4)) = 20 (even)
    assert mylib.integration_error_bound_N("simpson", 1, 2, 24.0, 1e-6) == 20
    # trapezoidal: N = ceil(sqrt(L^3 M2 / (12 eps))) = ceil(sqrt(2/(12e-6))) = 409
    assert mylib.integration_error_bound_N("trapezoidal", 1, 2, 2.0, 1e-6) == 409


def test_error_bound_simpson_rounds_to_even():
    n = mylib.integration_error_bound_N("simpson", 0, 1, 1.0, 1e-6)
    assert n % 2 == 0


def test_error_bound_rejects_unknown_method():
    with pytest.raises(ValueError):
        mylib.integration_error_bound_N("gauss", 0, 1, 1.0)


# ---------------------------------------------------------------------------
# Gaussian quadrature
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", range(1, 7))
def test_legendre_tables_match_numpy(n):
    nodes, weights = leggauss(n)
    assert mylib.GAUSS_LEGENDRE[n][0] == pytest.approx(list(nodes), abs=1e-14)
    assert mylib.GAUSS_LEGENDRE[n][1] == pytest.approx(list(weights), abs=1e-14)


@pytest.mark.parametrize("n", range(1, 7))
def test_laguerre_tables_match_numpy(n):
    nodes, weights = laggauss(n)
    assert mylib.GAUSS_LAGUERRE[n][0] == pytest.approx(list(nodes), abs=1e-12)
    assert mylib.GAUSS_LAGUERRE[n][1] == pytest.approx(list(weights), abs=1e-11)


def test_gauss_legendre_exact_through_degree_2n_minus_1():
    for n in range(1, 7):
        for d in range(2 * n):
            f = lambda x, d=d: x ** d  # noqa: E731
            exact = 0.0 if d % 2 == 1 else 2.0 / (d + 1)
            val = mylib.gaussian_quadrature(f, -1.0, 1.0, n, "legendre")
            assert val == pytest.approx(exact, abs=1e-12), (n, d)


def test_gauss_laguerre_exact_for_moment_factorials():
    for n in range(1, 7):
        for d in range(2 * n):
            f = lambda x, d=d: x ** d  # noqa: E731
            val = mylib.gaussian_quadrature(f, 0.0, 0.0, n, "laguerre")
            assert val == pytest.approx(math.factorial(d), rel=1e-10), (n, d)


def test_gauss_laguerre_moment_11_with_n6():
    # 11! = 39916800 — exact for the 6-node rule (2n-1 = 11)
    val = mylib.gaussian_quadrature(lambda x: x ** 11, 0.0, 0.0, 6, "laguerre")
    assert val == pytest.approx(39916800.0, rel=1e-9)


def test_gauss_legendre_interval_transform():
    # int_0^pi sin(x) dx = 2, via the [-1,1] -> [0, pi] map
    val = mylib.gaussian_quadrature(math.sin, 0.0, math.pi, 5, "legendre")
    assert val == pytest.approx(2.0, abs=1e-6)


def test_lab_problem_gauss_legendre_4point():
    # week06c q2: the committed value, reproduced exactly
    f = lambda x: x * x / (1.0 + x ** 4)  # noqa: E731
    val = mylib.gaussian_quadrature(f, -1.0, 1.0, 4, "legendre")
    assert val == pytest.approx(0.4816354816354816, abs=1e-14)
    # and it sits close to the true value 0.4874954944 (the problem statement)
    assert abs(val - 0.4874954943993610) < 1e-2


def test_lab_problem_laguerre_5point():
    # week06c q3: int e^-x/(1+x) on [0, inf)
    val = mylib.gaussian_quadrature(lambda x: 1.0 / (1.0 + x), 0.0, 0.0, 5, "laguerre")
    # 5 nodes give ~1.3e-3 accuracy on this non-polynomial integrand
    assert val == pytest.approx(0.5963473623231941, abs=2e-3)


def test_gaussian_quadrature_rejects_bad_arguments():
    with pytest.raises(ValueError):
        mylib.gaussian_quadrature(math.sin, -1, 1, 7)
    with pytest.raises(ValueError):
        mylib.gaussian_quadrature(math.sin, -1, 1, 4, "chebyshev")
