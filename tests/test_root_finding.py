"""Unit tests for the root-finding routines."""

from __future__ import annotations

import math

import pytest

import mylib

# ---------------------------------------------------------------------------
# Bracketing
# ---------------------------------------------------------------------------

def test_bracket_root_expands_from_bad_interval():
    f = lambda x: -x - math.cos(x)  # noqa: E731
    a, b = mylib.bracket_root(f, 2.0, 3.0)
    assert f(a) * f(b) < 0.0
    # the (unique) root near -0.739 must lie inside the expanded bracket
    assert a < -0.739085 < b


def test_bracket_root_accepts_good_interval():
    f = lambda x: x * x - 2.0  # noqa: E731
    a, b = mylib.bracket_root(f, 1.0, 2.0)
    assert (a, b) == (1.0, 2.0)


def test_bisection_converges_to_sqrt2():
    f = lambda x: x * x - 2.0  # noqa: E731
    root, it = mylib.bisection(f, 1.0, 2.0, accuracy=1e-10)
    assert root == pytest.approx(math.sqrt(2.0), abs=1e-9)
    assert it <= 60


def test_bisection_iteration_count_matches_log2():
    # interval width 1, accuracy 1e-8 -> ~log2(1e8) ~ 27 halvings
    f = lambda x: x * x - 2.0  # noqa: E731
    _, it = mylib.bisection(f, 1.0, 2.0, accuracy=1e-8)
    assert 20 <= it <= 40


def test_regula_falsi_converges_to_sqrt2():
    f = lambda x: x * x - 2.0  # noqa: E731
    root, it = mylib.regula_falsi(f, 1.0, 2.0, accuracy=1e-10)
    assert root == pytest.approx(math.sqrt(2.0), abs=1e-8)
    assert it <= 100


def test_bisection_lab_problem():
    # week05a q2: log(x/2) - sin(5x/2) on [1.0, 1.5], root ~ 1.40193
    f = lambda x: math.log(x / 2.0) - math.sin(5.0 * x / 2.0)  # noqa: E731
    root, _ = mylib.bisection(f, 1.0, 1.5, accuracy=1e-12)
    assert root == pytest.approx(1.401929931614613, abs=1e-9)


# ---------------------------------------------------------------------------
# Fixed point
# ---------------------------------------------------------------------------

def test_fixed_point_converges():
    g = lambda x: math.sqrt(2.0 * x + 3.0)  # noqa: E731
    root, it = mylib.fixed_point(g, 1.0, accuracy=1e-10)
    assert root == pytest.approx(3.0, abs=1e-8)
    assert 0 < it <= 30


def test_fixed_point_reports_iteration_count():
    g = lambda x: 0.5 * (x + 2.0 / x)  # noqa: E731  # Newton for sqrt(2) in
    # fixed-point clothing; |g'| small near root
    root, it = mylib.fixed_point(g, 1.0, accuracy=1e-12)
    assert root == pytest.approx(math.sqrt(2.0), abs=1e-10)
    assert it <= 20


# ---------------------------------------------------------------------------
# Newton
# ---------------------------------------------------------------------------

def test_newton_analytic_vs_numeric_derivative():
    f = lambda x: 3.0 * x + math.sin(x) - math.exp(x)  # noqa: E731
    df = lambda x: 3.0 + math.cos(x) - math.exp(x)  # noqa: E731
    r_an, it_an = mylib.newton_raphson(f, df, x0=0.0, accuracy=1e-10)
    r_num, it_num = mylib.newton_raphson(f, None, x0=0.0, accuracy=1e-10)
    assert r_an == pytest.approx(0.3604217029603245, abs=1e-10)
    assert r_num == pytest.approx(r_an, abs=1e-7)
    assert it_an <= 10 and it_num <= 10


def test_newton_quadratic_convergence():
    # error should roughly square each step: count iterations from a
    # moderate start to 1e-14 and check it is tiny (quadratic signature)
    f = lambda x: x * x - 2.0  # noqa: E731
    df = lambda x: 2.0 * x  # noqa: E731
    _, it = mylib.newton_raphson(f, df, x0=2.0, accuracy=1e-14)
    assert it <= 8  # bisection would need ~47


def test_newton_zero_derivative_raises():
    # a supplied (singular) derivative must be caught by the guard
    f = math.sin
    with pytest.raises(ValueError, match="derivative"):
        mylib.newton_raphson(f, lambda x: 0.0, x0=2.0)


def test_newton_multivariate_system():
    F = lambda x, y: [3 * x * x - 2 * x * y - 3, 3 * y * y - 4 * x * y]  # noqa: E731
    J = lambda x, y: mylib.jacobian(F, [x, y])  # noqa: E731
    x, it = mylib.newton_raphson_system(F, J, [2.5, 3.1], accuracy=1e-12)
    assert x == pytest.approx([3.0, 4.0], abs=1e-8)
    assert it <= 15
    assert max(abs(v) for v in F(*x)) < 1e-10


# ---------------------------------------------------------------------------
# Finite differences
# ---------------------------------------------------------------------------

def test_central_difference_first_derivative():
    f = math.sin
    assert mylib.finite_difference_derivative(f, 0.7) == pytest.approx(math.cos(0.7), abs=1e-9)


def test_central_difference_second_derivative():
    f = lambda x: math.exp(x)  # noqa: E731
    assert mylib.finite_difference_derivative(f, 0.3, order=2) == \
        pytest.approx(math.exp(0.3), abs=1e-6)


def test_jacobian_matches_analytic():
    F = lambda x, y: [x * x + y, x * y]  # noqa: E731
    J = mylib.jacobian(F, [2.0, 3.0])
    # dF1/dx = 2x = 4, dF1/dy = 1, dF2/dx = y = 3, dF2/dy = x = 2
    expected = [[4.0, 1.0], [3.0, 2.0]]
    for i in range(2):
        for j in range(2):
            assert J[i][j] == pytest.approx(expected[i][j], abs=1e-6)
