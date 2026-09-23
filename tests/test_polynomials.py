"""Unit tests for polynomial evaluation, Laguerre's method and deflation."""

from __future__ import annotations

import numpy as np
import pytest

import mylib


def test_horner_evaluation_matches_numpy():
    rng = np.random.default_rng(5)
    coeffs = rng.normal(size=8).tolist()
    for x in (-3.7, -1.0, 0.0, 0.42, 2.9):
        assert mylib.polynomial_value(coeffs, x) == \
            pytest.approx(np.polyval(coeffs, x), rel=1e-12)


def test_derivative_coefficients_match_numpy():
    coeffs = [1.0, -6.0, 11.0, -6.0]
    d1 = mylib.polynomial_first_derivative(coeffs)
    d2 = mylib.polynomial_second_derivative(coeffs)
    assert d1 == pytest.approx([3.0, -12.0, 11.0])
    assert d2 == pytest.approx([6.0, -12.0])


def _real_roots_of(coeffs: list[float]) -> list[float]:
    roots = np.roots(coeffs)
    return sorted(r.real for r in roots if abs(r.imag) < 1e-8)


def test_laguerre_cubic_known_roots():
    coeffs = [1.0, -6.0, 11.0, -6.0]  # (x-1)(x-2)(x-3)
    roots, _ = mylib.laguerre_roots(coeffs, b0=2.5)
    assert sorted(roots) == pytest.approx([1.0, 2.0, 3.0], abs=1e-7)


@pytest.mark.parametrize("coeffs, b0", [
    ([1, -1, -7, 1, 6], 1.0),          # P1: roots -2, -1, 1, 3
    ([1, 0, -5, 0, 4], 1.0),           # P2: roots -2, -1, 1, 2
    ([2, 0, -19.5, 0.5, 13.5, -4.5], 1.0),  # P3: -3, -1, 0.5 (double), 3
])
def test_laguerre_lab_polynomials(coeffs, b0):
    roots, _ = mylib.laguerre_roots(coeffs, b0)
    ref = _real_roots_of(coeffs)
    mine = sorted(roots)
    # P3 has a double root: numpy splits it into a complex pair, so compare
    # with the deflated polynomial instead (exact roots known analytically).
    for r in mine:
        assert abs(mylib.polynomial_value(coeffs, r)) < 1e-5
    if coeffs != [2, 0, -19.5, 0.5, 13.5, -4.5]:
        assert mine == pytest.approx(ref, abs=1e-6)
    else:
        # double root 0.5 found as two nearby values
        near = [r for r in mine if abs(r - 0.5) < 1e-3]
        assert len(near) == 2
        assert sorted(mine) == pytest.approx(sorted([-3.0, -1.0, 0.5, 0.5, 3.0]), abs=1e-3)


def test_laguerre_stops_at_real_roots_only():
    # x^4 - x^2 + 1 has no real roots: the real-only method must fail loudly
    # rather than returning a bogus number.
    with pytest.raises(ValueError):
        mylib.laguerre_roots([1.0, 0.0, -1.0, 0.0, 1.0], b0=1.0)


def test_synthetic_division_quotient():
    # (x^3 - 6x^2 + 11x - 6) / (x - 2) = x^2 - 4x + 3
    q = mylib.synthetic_division([1.0, -6.0, 11.0, -6.0], 2.0)
    assert q == pytest.approx([1.0, -4.0, 3.0], abs=1e-12)


def test_synthetic_division_rejects_non_root():
    with pytest.raises(ValueError, match="remainder"):
        mylib.synthetic_division([1.0, -6.0, 11.0, -6.0], 5.0)
