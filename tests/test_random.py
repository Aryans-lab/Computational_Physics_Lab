"""Unit tests for the LCG, distributions, and MyComplex."""

from __future__ import annotations

import math

import pytest

import mylib

# ---------------------------------------------------------------------------
# LCG
# ---------------------------------------------------------------------------

def test_lcg_full_period():
    # Hull-Dobell full-period check: m draws after a reset cover all m states
    # with no repetition (the defining property of a full-period LCG).
    mylib.myrand_reset(1)
    period = mylib.LCG_M
    seen = set()
    for _ in range(period):
        seen.add(mylib.myrand())
    assert len(seen) == period


def test_lcg_hull_dobell_conditions():
    a, c, m = mylib.LCG_A, mylib.LCG_C, mylib.LCG_M
    assert m & (m - 1) == 0            # m is a power of two
    assert a % 4 == 1                  # a = 1 mod 4
    assert c % 2 == 1                  # c odd


def test_lcg_uniform_moments():
    mylib.myrand_reset(0)
    n = 100_000
    draws = [mylib.myrand() for _ in range(n)]
    mean = sum(draws) / n
    var = sum((d - mean) ** 2 for d in draws) / n
    assert abs(mean - 0.5) < 0.01          # exact mean 0.5
    assert abs(var - 1.0 / 12.0) < 0.01    # exact variance 1/12


def test_lcg_reseed_is_deterministic():
    mylib.myrand_reset(42)
    seq1 = [mylib.myrand() for _ in range(50)]
    mylib.myrand_reset(42)
    seq2 = [mylib.myrand() for _ in range(50)]
    assert seq1 == seq2
    mylib.myrand_reset(43)
    seq3 = [mylib.myrand() for _ in range(50)]
    assert seq1 != seq3


def test_random_uniform_range():
    mylib.myrand_reset(1)
    vals = [mylib.random_uniform(2.0, 5.0) for _ in range(10_000)]
    assert all(2.0 <= v < 5.0 for v in vals)


def test_random_exponential_moments():
    mylib.myrand_reset(1)
    lam = 2.0
    n = 50_000
    samples = [mylib.random_exponential(lam) for _ in range(n)]
    mean = sum(samples) / n
    assert abs(mean - 1.0 / lam) < 0.02  # within a couple of sigma of 0.5


# ---------------------------------------------------------------------------
# MyComplex
# ---------------------------------------------------------------------------

def test_mycomplex_arithmetic():
    c1 = mylib.MyComplex(1.3, -2.2)
    c2 = mylib.MyComplex(-0.8, 1.7)
    s = c1 + c2
    assert s.real == pytest.approx(0.5, abs=1e-12)
    assert s.imag == pytest.approx(-0.5, abs=1e-12)
    d = c1 - c2
    assert d.real == pytest.approx(2.1, abs=1e-12)
    assert d.imag == pytest.approx(-3.9, abs=1e-12)
    p = c1 * c2
    assert p.real == pytest.approx(2.7, abs=1e-12)
    assert p.imag == pytest.approx(3.97, abs=1e-12)
    assert abs(c1) == pytest.approx(math.hypot(1.3, -2.2))
    # equality works too (exact construction, no float drift)
    assert mylib.MyComplex(1.5, -2.5) == mylib.MyComplex(1.5, -2.5)
    assert mylib.MyComplex(1.0, 0.0) == 1.0


def test_mycomplex_scalar_multiplication():
    c = mylib.MyComplex(2.0, 3.0)
    assert (c * 2) == mylib.MyComplex(4.0, 6.0)
