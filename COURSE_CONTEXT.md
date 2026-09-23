# Computational Physics Laboratory (PHY341 / PHY745) — Course Context

**Course:** PHY341 / PHY745 — Computational Physics
**Program:** Integrated M.Sc. Physics, Semester V
**Institution:** School of Physical Sciences, National Institute of Science Education and Research (NISER), Bhubaneswar
**Semester:** 2025–26

This document records the rules the coursework was written under, the
layout of the repository, the mathematical reference for the library, and
the standard code pattern every driver follows. It is the provenance layer:
each `assignments/weekNN_*/` folder is a complete, self-contained submission
that obeys these rules exactly.

---

## 1. Course Philosophy & Mandatory Ground Rules

The course teaches foundational numerical algorithms from scratch, without
black-box libraries. Submissions are audited under the following
non-negotiable constraints:

1. **Zero high-level library policy**
   - Functions from `scipy`, `numpy.linalg`, or the standard-library
     `random` module must **never** be used for core numerical algorithms
     (matrix inversion, linear solves, random number generation, root
     finding, numerical integration).
   - Only elementary math (`math.sqrt`, `math.sin`, `math.log`, …) and
     basic plotting (`matplotlib.pyplot`) are permitted; NumPy/SciPy may
     not touch the numerics at all.
   - Standard data structures are native Python lists and nested lists.
2. **Single central library rule (`mylib.py`)**
   - All numerical routines (LCG, LU, Cholesky, Jacobi, Gauss-Seidel,
     Newton-Raphson, Laguerre, Simpson, Gaussian quadrature, …) live in
     `mylib.py` in the same directory as the front code.
   - Front-end scripts (`question1.py`, …) stay lightweight: read input,
     call the library, write outputs/plots. Nothing else.
3. **External file I/O (non-interactive)**
   - Inputs are never hardcoded or taken via `input()`; they come from
     external ASCII files in `data/`.
   - Outputs are written to `output/qN_output.txt`, mirroring the console
     run; explanatory physics notes are appended below the raw output.
   - Plots go to `figures/` (saved, never displayed interactively).
4. **Program headers & attribution**
   - Every script opens with a comment block: problem statement, usage,
     author.
5. **Style & academic integrity**
   - Code reads like that of a rigorous physics student: clean, commented,
     direct, unbloated.

> **Repository convention on top of the course rules:** this repository
> keeps one canonical `mylib.py` at the root (the single source of truth,
> type-annotated and unit-tested), and `scripts/bundle.py` copies it
> verbatim into each week folder. Each folder therefore stays
> independently submittable — exactly as the course prescribes — while the
> root copy is what the test suite and CI validate.

---

## 2. Directory Structure & Assignments Map

```
.
├── mylib.py                     # canonical library (pure stdlib)
├── COURSE_CONTEXT.md            # this document
├── scripts/
│   ├── bundle.py                # mylib.py -> every week folder (idempotent)
│   └── run_all.py               # run all 36 drivers, report ok/FAIL
├── tests/                       # unit tests + golden-output regression
├── Makefile                     # make all = bundle + run + lint + test
├── .github/workflows/ci.yml     # Python 3.10 & 3.12
└── assignments/
    ├── week00_warmup/           # loops, series (GP/HP), complex numbers, matrix I/O
    │   ├── question1.py         (sum of first N even numbers, factorial of M)
    │   ├── question2.py         (GP & HP series sums, no closed form)
    │   ├── question3.py         (matrix multiplication from files: AB, BC, D·C)
    │   └── question5.py         (MyComplex: +, -, *, |z|)
    ├── week01_random/           # PRNG (LCG), distributions, simulations
    │   ├── question1.py         (chaos map vs LCG, lag-correlation plots)
    │   ├── question2.py         (LCG, full-period audit, k-lag correlation)
    │   ├── question3.py         (Monte-Carlo pi, residual analysis)
    │   ├── question4.py         (radioactive decay A -> B -> C)
    │   └── question5.py         (inverse-transform: exponential deviates)
    ├── week02_gauss_jordan_lu/  # direct solvers: Gauss-Jordan & LU
    │   ├── question1.py         (pivot-stress GJ system)
    │   ├── question2.py         (3-variable system via Gauss-Jordan)
    │   └── question3.py         (LU verification: P A = L U)
    ├── week03_lu_cholesky/      # forward-backward & Cholesky
    │   ├── question1.py         (Cholesky routine with symmetry check)
    │   ├── question2.py         (6x6 system via LU forward-backward)
    │   └── question3.py         (Cholesky solve: L y = b, L^T x = y)
    ├── week04_jacobi_gs/        # iterative solvers
    │   ├── question1.py         (benchmark: Jacobi vs G-S vs SOR)
    │   ├── question2.py         (Jacobi vs G-S on the 6x6 system)
    │   └── question3.py         (Poisson model: G-S vs SOR, omega = 1.57)
    ├── week05a_bisection_regf/  # root finding: bracketing
    │   ├── question1.py         (validation + auto-bracketing demo)
    │   ├── question2.py         (root of ln(x/2) - sin(5x/2) = 0)
    │   └── question3.py         (root of -x - cos(x) = 0 from [2, 3])
    ├── week05b_fixedpoint_newton/  # fixed-point & Newton
    │   ├── question1.py         (validation: 1-D & multivariate Newton)
    │   ├── question2.py         (rate comparison: bi. / rf / NR, analytic vs numeric f')
    │   ├── question3.py         (fixed-point for both roots of x^2 - 2x - 3)
    │   └── question4.py         (coupled nonlinear system, numerical Jacobian)
    ├── week05c_polynomial_roots/   # Horner, Laguerre, deflation
    │   ├── question1.py         (library validation on a known cubic)
    │   └── question2.py         (all real roots of the lab's quartic & quintic)
    ├── week06a_midpoint_trapezoid/  # closed Newton-Cotes
    │   ├── question1.py         (midpoint & trapezoidal validation, O(h^2) table)
    │   └── question2.py         (1/x, x cos x, x arctan x at N = 4,8,15,20)
    ├── week06b_simpson_montecarlo/  # Simpson & Monte Carlo
    │   ├── question1.py         (Simpson exactness on cubics; MC sigma check)
    │   ├── question2.py         (6 d.p. integration, N from error bounds)
    │   └── question3.py         (MC of sin^2 on [-1,1], value & deviation vs N)
    └── week06c_gaussian_quadrature/  # Gauss-Legendre & Gauss-Laguerre
        ├── question1.py         (unified routine, tables n = 1..6, order-of-precision proof)
        ├── question2.py         (4-pt Gauss-Legendre on x^2/(1+x^4) vs Simpson)
        ├── question3.py         (5-pt Gauss-Laguerre on int_0^inf e^-x/(1+x))
        └── question4.py         (3-pt Gauss-Laguerre on int_0^inf e^-x x^4 = 4!)
```

Every week folder additionally contains: `README.md`, `data/` (inputs),
`output/qN_output.txt` (committed golden outputs), `figures/` (PNGs), and
`mylib.py` (bundled copy, gitignored).

---

## 3. Mathematical & Algorithmic Reference

### 3.1 Pseudo-Random Numbers (`myrand`, LCG)
- **Recurrence:** $x_{i+1} = (a \, x_i + c) \pmod m$
- **Parameters (full-period, Hull-Dobell):** $a = 1103515245$,
  $c = 12345$, $m = 32768 = 2^{15}$.
- **Uniform $[A, B)$:** $X = A + (B - A)\,\xi$, $\xi = x_i / m$.
- **Exponential $\lambda e^{-\lambda y}$:** inverse-CDF
  $y = -\lambda^{-1}\ln u$, $u \in (0, 1]$.

### 3.2 Direct Linear Solvers ($A x = b$)
- **Gauss-Jordan:** $[A \mid b] \to [I \mid x]$ via RREF, with
  partial row pivoting (scale-relative tolerance) and explicit
  classification of inconsistent / dependent / non-square systems.
- **Matrix inversion:** $[A \mid I] \to [I \mid A^{-1}]$.
- **LU (Doolittle, partial pivoting):** $P A = L U$, $L_{ii} = 1$.
  $$U_{ij} = A_{ij} - \sum_{k < i} L_{ik} U_{kj}, \qquad
    L_{ji} = \frac{1}{U_{ii}} \Big( A_{ji} - \sum_{k < i} L_{jk} U_{ki} \Big)$$
  At each pivot, the already-computed multipliers $L[i][k], k < i$ move with
  the row swap; swapping whole $L$ rows would dislocate the unit diagonal.
- **Forward-backward:** solve $L y = b$, then $U x = y$.
- **Free determinant:** $\det A = \operatorname{sign}(P)\prod_i U_{ii}$,
  with $\operatorname{sign}(P)$ from the permutation's cycle decomposition.

### 3.3 Cholesky ($A = L L^T$)
- **Applicability:** symmetric ($A = A^T$) and positive-definite.
  Both conditions are checked *before* the first square root.
- **Algorithm:**
  $$L_{ii} = \sqrt{A_{ii} - \sum_{k < i} L_{ik}^2}, \qquad
    L_{ji} = \frac{1}{L_{ii}} \Big( A_{ji} - \sum_{k < i} L_{jk} L_{ik} \Big)$$
- **Cost:** $\tfrac{1}{3}n^3$ flops — about half of a full LU.

### 3.4 Iterative Solvers (Jacobi, Gauss-Seidel, SOR)
- **Convergence:** guaranteed for strictly diagonally dominant or
  symmetric positive-definite $A$ (Gauss-Seidel).
- **Jacobi:** $x_i^{(k+1)} = \frac{1}{A_{ii}}\big( b_i - \sum_{j \ne i} A_{ij} x_j^{(k)} \big)$
- **Gauss-Seidel:** same, but reuses $x_j^{(k+1)}$ for $j < i$ in place.
- **SOR:** $x_i^{(k+1)} = (1-\omega) x_i^{(k)} + \frac{\omega}{A_{ii}}
  \big( b_i - \sum_{j < i} A_{ij} x_j^{(k+1)} - \sum_{j > i} A_{ij} x_j^{(k)} \big)$,
  with $\omega = 1$ recovering Gauss-Seidel; $1 < \omega < 2$ accelerates.
  For the tridiagonal Poisson model $(-1, 2, -1)$ the optimal value is
  $\omega_{\text{opt}} = 2/(1 + \sin(\pi/(n+1)))$.

### 3.5 Root Finding
- **Auto-bracketing:** if $f(a) f(b) > 0$, push out the side with smaller
  $|f|$ (step growing 10% per push) until a sign change appears.
- **Bisection:** $c = (a+b)/2$; linear, error halves each step.
- **Regula falsi:** secant
  $c = b - \dfrac{f(b)(b-a)}{f(b) - f(a)}$; superlinear, with the known
  "stuck endpoint" pathology.
- **Fixed point:** $x_{n+1} = g(x_n)$; converges iff $|g'(x^*)| < 1$.
- **Newton (1-D):** $x_{n+1} = x_n - f(x_n)/f'(x_n)$; $f'$ analytic or
  central-difference; quadratic near simple roots.
- **Newton (multivariate):** $x^{(k+1)} = x^{(k)} + \delta x$ where
  $J(x^{(k)})\,\delta x = -F(x^{(k)})$ is *solved* (Gauss-Jordan) rather
  than inverted; $J_{ij} = \partial f_i / \partial x_j$ by central
  differences.

### 3.6 Polynomial Roots
- **Horner:** $b_n = a_n$, $b_k = a_k + x\, b_{k+1}$ — $O(n)$ evaluation;
  the same scheme yields $P'(x), P''(x)$ by coefficient differentiation.
- **Laguerre:** with $G = P'/P$, $H = G^2 - P''/P$,
  $$a = \frac{n}{G \pm \sqrt{(n-1)(nH - G^2)}}$$
  (sign chosen to avoid cancellation); update $x \leftarrow x - a$.
  Cubic near simple roots; **linear** near multiple roots.
- **Deflation:** synthetic division by $(x - r)$ (remainder must vanish);
  repeat until linear, then solve the last factor exactly.
- **Contract:** real roots only — a step with no real solution raises
  `ValueError` instead of returning a phantom.

### 3.7 Quadrature
- **Composite midpoint:** $M_N = h \sum_{i=0}^{N-1} f(a + (i+\tfrac12)h)$,
  $\lVert E \rVert \le \frac{(b-a)^3}{24 N^2} \max |f''|$.
- **Composite trapezoid:** $T_N = \frac{h}{2}\big[ f(a) + 2\sum_{i=1}^{N-1} f(a+ih) + f(b) \big]$,
  $\lVert E \rVert \le \frac{(b-a)^3}{12 N^2} \max |f''|$.
- **Composite Simpson 1/3** ($N$ even): $S_N = \frac{h}{3}\big[ f(a) + 4\sum_{odd} f + 2\sum_{even} f + f(b) \big]$,
  $\lVert E \rVert \le \frac{(b-a)^5}{180 N^4} \max |f^{(4)}|$;
  identity $S_{2N} = \tfrac23 M_N + \tfrac13 T_N$.
- **Uniform Monte Carlo:** $F_N = (b-a)\,\bar f$,
  $\sigma_I = \dfrac{(b-a)\,\sigma_f}{\sqrt N}$; rate $O(N^{-1/2})$,
  dimension-independent.
- **Gaussian quadrature:** $n$ nodes exact for $\deg \le 2n-1$.
  - *Legendre:* $[-1,1]$ mapped by $x = \frac{b-a}{2}t + \frac{b+a}{2}$.
  - *Laguerre:* $\int_0^\infty e^{-t} f(t)\,dt \approx \sum_i w_i f(t_i)$ —
    the weight is the quadrature's own, so the supplied $f$ excludes
    $e^{-t}$.
  - Tables for $n = 1 \dots 6$ are hardwired in `mylib.py` and validated
    against NumPy's `leggauss` / `laggauss` in the test suite.

---

## 4. Library Reference

The authoritative, typed reference is the module docstring of
[`mylib.py`](mylib.py) itself (kept in sync by the test suite and CI). The
surface, in module order:

1. **RNG & distributions** — `myrand`, `myrand_reset`, `random_uniform`,
   `random_exponential`, plus the `LCG_A / LCG_C / LCG_M` constants.
2. **Numbers & vectors** — `MyComplex` (with `+`, `-`, `*`, `abs`, `==`),
   `dot_product_vector`, `vector_norm`.
3. **Matrix utilities & I/O** — `read_matrix_from_file`,
   `read_vector_from_file`, `write_matrix_to_file`, `print_matrix`,
   `matrix_multiply`, `matrix_transpose`, `matrix_add`, `matrix_sub`,
   `matrix_vector_multiply`, `matrix_residual`, `is_symmetric`,
   `is_diagonally_dominant`.
4. **Direct solvers** — `gauss_jordan_elimination_augmented`,
   `gauss_jordan_inverse`, `lu_decomposition` (returns `L, U, perm`),
   `lu_forback`, `matrix_determinant_lu`, `cholesky_decomposition`,
   `cholesky_forback`.
5. **Iterative solvers** — `jacobi_it`, `gauss_seidel`,
   `sor_gauss_seidel` (each returns `(solution, iterations)`).
6. **Root finding** — `bracket_root`, `bisection`, `regula_falsi`,
   `fixed_point`, `finite_difference_derivative`, `partial_derivative`,
   `jacobian`, `newton_raphson`, `newton_raphson_system`.
7. **Polynomials** — `polynomial_value`, `polynomial_first_derivative`,
   `polynomial_second_derivative`, `laguerre`, `synthetic_division`,
   `laguerre_roots`.
8. **Quadrature** — `midpoint`, `trapezoidal`, `simpson` (returns
   `(integral, evaluations)`), `monte_carlo` (returns
   `(integral, sigma_f, sigma_I)`, seeded),
   `integration_error_bound_N`, `gaussian_quadrature` (one routine for
   `method='legendre' | 'laguerre'`), plus the `GAUSS_LEGENDRE` and
   `GAUSS_LAGUERRE` tables.

---

## 5. Standard Code Template (front code)

Every driver follows this shape — I/O and presentation only, numerics in
the library:

```python
"""
questionN.py — weekNN_topic
---------------------------
Problem : <one-line statement from the assignment sheet>
Usage   : python questionN.py [output_file]
          Defaults: output/qN_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import some_routine, matrix_residual   # numerics: library only

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "qN_output.txt"

    # 1. read inputs (non-interactive, from data/ relative to this script)
    A = mylib.read_matrix_from_file(BASE_DIR / "data" / "asgn3_mat1")
    b = mylib.read_vector_from_file(BASE_DIR / "data" / "asgn3_vec1")

    # 2. numerics: call the library
    x = lu_forback(A, b)
    res = matrix_residual(A, x, b)

    # 3. report: verified output to file (and console)
    lines = [f"  x[{i}] = {v:.10f}" for i, v in enumerate(x)]
    lines.append(f"  residual ||Ax - b||_2 = {res:.3e}")
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
```

Plotting drivers additionally set `matplotlib.use("Agg")` before importing
`pyplot` and save to `figures/` — never `plt.show()`.

---

## 6. Post-Midsem Roadmap

Topics scheduled for the second half of the semester (not yet part of this
repository):

1. **ODEs** — Euler (forward/backward), Heun/Milne predictor-corrector,
   classical RK4 for coupled systems, IVPs and BVPs (shooting, finite
   differences).
2. **Least-squares fitting** — linear/polynomial regression via normal
   equations, Gauss-Newton, Levenberg-Marquardt.
3. **PDEs & eigenvalues (if time permits)** — FTCS / Crank-Nicolson for
   diffusion, the wave equation, power method and QR iteration.
