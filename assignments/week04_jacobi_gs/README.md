# Week 04 — Stationary Iterative Methods: Jacobi, Gauss-Seidel, SOR

Direct $O(n^3)$ solvers hit a wall for the huge, sparse systems that emerge
from discretised field equations (Poisson, heat, Schrödinger). Iterative
methods trade one big factorisation for many cheap sweeps — and for sparse
matrices each sweep is $O(n)$.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Benchmarks Jacobi, Gauss-Seidel, and SOR ($\omega = 1.3$) on the diagonally-dominant 6×6 system `data/asgn4_mat1` | iteration counts + residuals for all three; G-S ≈ ½ the Jacobi sweeps |
| 2 | `question2.py` | The same 6×6 system, Jacobi vs Gauss-Seidel head-to-head | speed-up ≈ 2× in iteration count |
| 3 | `question3.py` | The 10×10 tridiagonal Poisson model, Gauss-Seidel vs SOR with $\omega = 1.57$ | **339 → 58 iterations** (≈6× faster); both hit the exact solution $(5,9,12,14,15,15,14,12,9,5)$ |

## What's worth reading

- **Why the Poisson model crawls.** The tridiagonal matrix
  $(-1, 2, -1)$ is only *weakly* diagonally dominant: the spectral radius of
  the G-S iteration matrix sits at $1 - O(n^{-2})$, so convergence is
  glacial for large $n$.
- **Why 1.57.** For this model the optimal SOR parameter has the analytic
  form $\omega_{\text{opt}} = \dfrac{2}{1 + \sin(\pi/(n+1))} \approx 1.56$
  at $n = 10$ — the driver's $\omega = 1.57$ is essentially on target, which
  is exactly why the run is 6× shorter.
- **Convergence is a theorem, not a hope.** Both G-S and Jacobi converge for
  strictly diagonally dominant matrices (and more generally SPD ones); the
  library's `is_diagonally_dominant` makes the hypothesis checkable.

## Running

```bash
python question1.py   # zero arguments: writes output/q1_output.txt
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
