"""Unit tests for the direct and iterative linear solvers.

Reference values come from NumPy (a test-only dependency) and from exact
integer solutions of the lab's own test systems.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import mylib

# ---------------------------------------------------------------------------
# Gauss-Jordan
# ---------------------------------------------------------------------------

def test_gauss_jordan_known_solution():
    A = [[0.0, 2.0, 5.0], [3.0, -1.0, 2.0], [1.0, -1.0, 3.0]]
    b = [1.0, -2.0, 3.0]
    aug = [A[i][:] + [b[i]] for i in range(3)]
    x = mylib.gauss_jordan_elimination_augmented(aug)
    assert x == pytest.approx([-2.0, -2.0, 1.0], abs=1e-12)
    assert mylib.matrix_residual(A, x, b) < 1e-12


def test_gauss_jordan_pivot_in_later_column():
    # Column 1 is fine, but column 2 has a zero on the diagonal without a
    # row exchange: non-pivoted codes fail here.
    A = [[1.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
    b = [1.0, 2.0, 3.0]
    aug = [A[i][:] + [b[i]] for i in range(3)]
    x = mylib.gauss_jordan_elimination_augmented(aug)
    assert x == pytest.approx([-2.0, 3.0, 2.0], abs=1e-12)


@pytest.mark.parametrize("augmented, message", [
    ([[1.0, 2.0, 3.0], [1.0, 2.0, 4.0]], "Inconsistent"),
    ([[1.0, 2.0, 5.0], [2.0, 4.0, 10.0]], "Dependent"),
    ([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], "Dependent"),
])
def test_gauss_jordan_singular_systems_raise(augmented, message):
    with pytest.raises(ValueError, match=message):
        mylib.gauss_jordan_elimination_augmented(augmented)


def test_gauss_jordan_rejects_wrong_shape():
    with pytest.raises(ValueError, match="shape"):
        mylib.gauss_jordan_elimination_augmented([[1.0, 1.0], [1.0, 1.0]])


def test_gauss_jordan_inverse():
    A = [[4.0, 7.0], [2.0, 6.0]]
    A_inv = mylib.gauss_jordan_inverse(A)
    product = mylib.matrix_multiply(A, A_inv)
    for i in range(2):
        for j in range(2):
            assert product[i][j] == pytest.approx(1.0 if i == j else 0.0, abs=1e-12)


def test_gauss_jordan_inverse_singular_raises():
    with pytest.raises(ValueError, match="singular"):
        mylib.gauss_jordan_inverse([[1.0, 2.0], [2.0, 4.0]])


# ---------------------------------------------------------------------------
# LU decomposition (partial pivoting)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(30))
def test_lu_factorisation_matches_numpy(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(2, 12))
    A = rng.normal(size=(n, n))
    L, U, perm = mylib.lu_decomposition(A.tolist(), return_perm=True)
    P = np.zeros((n, n))
    P[np.arange(n), perm] = 1.0
    assert np.allclose(np.array(L) @ np.array(U), P @ A, atol=1e-9)
    # L must be unit lower triangular
    assert np.allclose(np.array(L), np.tril(np.array(L), 0), atol=1e-15)
    assert np.allclose(np.diag(np.array(L)), 1.0, atol=1e-15)


@pytest.mark.parametrize("seed", range(20))
def test_lu_solve_matches_numpy(seed):
    rng = np.random.default_rng(1000 + seed)
    n = int(rng.integers(2, 12))
    A = rng.normal(size=(n, n))
    b = rng.normal(size=n)
    x = mylib.lu_forback(A.tolist(), b.tolist())
    assert np.allclose(A @ np.array(x), b, atol=1e-9)


def test_lu_zero_pivot_raises():
    # No pivot available in column 0 -> must raise, not divide by zero.
    with pytest.raises(ValueError, match="zero pivot"):
        mylib.lu_decomposition([[0.0, 1.0], [0.0, 1.0]])


# ---------------------------------------------------------------------------
# Determinants (sign-sensitive: catches permutation-sign bugs)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(30))
def test_determinant_matches_numpy(seed):
    rng = np.random.default_rng(7000 + seed)
    n = int(rng.integers(2, 13))
    A = rng.normal(size=(n, n))
    ref = np.linalg.det(A)
    mine = mylib.matrix_determinant_lu(A.tolist())
    assert math.copysign(1.0, mine) == math.copysign(1.0, ref)
    assert abs(mine - ref) <= 1e-6 * max(1.0, abs(ref))


def test_determinant_2x2_exact():
    assert mylib.matrix_determinant_lu([[4.0, 7.0], [2.0, 6.0]]) == pytest.approx(10.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Cholesky
# ---------------------------------------------------------------------------

def test_cholesky_spd_round_trip():
    rng = np.random.default_rng(42)
    n = 8
    M = rng.normal(size=(n, n))
    A = M @ M.T + n * np.eye(n)  # symmetric positive definite
    L = mylib.cholesky_decomposition(A.tolist())
    assert np.allclose(np.array(L) @ np.array(L).T, A, atol=1e-10)


def test_cholesky_known_factor():
    # The week03 lab system: A = L L^T with the committed L recovered.
    A = [[4.0, 1.0, 1.0, 1.0], [1.0, 3.0, -1.0, 1.0],
         [1.0, -1.0, 2.0, 0.0], [1.0, 1.0, 0.0, 2.0]]
    L = mylib.cholesky_decomposition(A)
    assert L[0][0] == pytest.approx(2.0)
    x = mylib.cholesky_forback(L, [3.0, 3.0, 1.0, 3.0])
    assert x == pytest.approx([0.0, 1.0, 1.0, 1.0], abs=1e-12)
    assert mylib.matrix_residual(A, x, [3.0, 3.0, 1.0, 3.0]) < 1e-12


def test_cholesky_rejects_non_symmetric():
    with pytest.raises(ValueError, match="symmetric"):
        mylib.cholesky_decomposition([[1.0, 2.0], [3.0, 4.0]])


def test_cholesky_rejects_non_positive_definite():
    with pytest.raises(ValueError, match="positive-definite"):
        mylib.cholesky_decomposition([[-1.0, 0.0], [0.0, 1.0]])
    # zero diagonal entry
    with pytest.raises(ValueError, match="positive-definite"):
        mylib.cholesky_decomposition([[0.0, 1.0], [1.0, 2.0]])


# ---------------------------------------------------------------------------
# Iterative solvers
# ---------------------------------------------------------------------------

A6 = [[4.0, -1.0, 0.0, -1.0, 0.0, 0.0],
      [-1.0, 4.0, -1.0, 0.0, -1.0, 0.0],
      [0.0, -1.0, 4.0, 0.0, 0.0, -1.0],
      [-1.0, 0.0, 0.0, 4.0, -1.0, 0.0],
      [0.0, -1.0, 0.0, -1.0, 4.0, -1.0],
      [0.0, 0.0, -1.0, 0.0, -1.0, 4.0]]
b6 = [2.0, 1.0, 2.0, 2.0, 1.0, 2.0]


def test_libraries_report_diagonal_dominance():
    assert mylib.is_diagonally_dominant(A6, strict=True)
    assert not mylib.is_diagonally_dominant([[1.0, 2.0], [1.0, 1.0]], strict=True)


@pytest.mark.parametrize("solver, name", [
    (mylib.jacobi_it, "jacobi"),
    (mylib.gauss_seidel, "gauss_seidel"),
])
def test_iterative_solvers_converge(solver, name):
    x, it = solver(A6, b6, tol=1e-10)
    assert 0 < it <= 1000
    assert mylib.matrix_residual(A6, x, b6) < 1e-8
    assert x == pytest.approx([1.0] * 6, abs=1e-5)


def test_sor_converges_faster_than_gauss_seidel():
    # Poisson model: weakly diagonally dominant -> SOR should win.
    n = 10
    A = [[2.0 if i == j else (-1.0 if abs(i - j) == 1 else 0.0)
          for j in range(n)] for i in range(n)]
    b = [1.0] * n
    _, it_gs = mylib.gauss_seidel(A, b, tol=1e-12, max_iterations=20000)
    x_sor, it_sor = mylib.sor_gauss_seidel(A, b, tol=1e-12,
                                           max_iterations=20000, omega=1.57)
    assert it_sor < it_gs / 2
    assert mylib.matrix_residual(A, x_sor, b) < 1e-10


def test_sor_rejects_bad_omega():
    with pytest.raises(ValueError, match="omega"):
        mylib.sor_gauss_seidel(A6, b6, omega=2.5)


def test_iterative_solver_zero_diagonal_raises():
    with pytest.raises(ValueError, match="non-zero diagonal"):
        mylib.jacobi_it([[0.0, 1.0], [1.0, 0.0]], [1.0, 1.0])


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def test_matrix_residual_detects_wrong_solution():
    A = [[4.0, 1.0], [1.0, 3.0]]
    good = mylib.lu_forback(A, [5.0, 4.0])
    assert mylib.matrix_residual(A, good, [5.0, 4.0]) < 1e-12
    assert mylib.matrix_residual(A, [2.0, 2.0], [5.0, 4.0]) > 1e-6


def test_read_write_matrix_round_trip(tmp_path):
    A = [[1.5, -2.0, 3.25], [0.0, 7.0, -1.0]]
    p = tmp_path / "mat.txt"
    mylib.write_matrix_to_file(A, str(p))
    back = mylib.read_matrix_from_file(str(p))
    for i in range(2):
        for j in range(3):
            assert back[i][j] == pytest.approx(A[i][j])
    vec = mylib.read_vector_from_file(str(p))
    for i, v in enumerate([1.5, 0.0]):
        assert vec[i] == pytest.approx(v)


def test_read_matrix_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        mylib.read_matrix_from_file("/nonexistent/path/mat.txt")


def test_matrix_multiply_shapes():
    assert mylib.matrix_multiply([[1.0, 2.0]], [[3.0], [4.0]]) == [[11.0]]
    with pytest.raises(ValueError, match="Dimension mismatch"):
        mylib.matrix_multiply([[1.0, 2.0]], [[1.0, 2.0]])
