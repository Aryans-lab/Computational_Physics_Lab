# Computational Physics Laboratory (PHY341 / PHY745) - Course Context & System Architecture

**Student:** Aryan Bandyopadhyay  
**Roll Number:** 2411014  
**Program:** Integrated MSc (3rd Year, Semester V), Physics  
**Institution:** School of Physical Sciences, National Institute of Science Education and Research (NISER), Bhubaneswar  
**Instructor:** Prof. Subhasish Basak (`sbasak@niser.ac.in`)  
**Teaching Assistants:** Hemant Lohumi, Chinmoy Samanta  
**Primary Repository Path:** `D:\NISER\Computational Physics Lab\`

---

## 1. Course Philosophy & Mandatory Ground Rules

This course is designed to teach foundational numerical algorithms from scratch without reliance on "black-box" libraries. Submissions are strictly audited by the instructor and TAs under the following non-negotiable constraints:

1. **Zero High-Level Library Policy:**
   - Functions from `scipy`, `numpy.linalg`, or standard library `random` must **never** be used for core numerical algorithms (e.g., matrix inversion, solving linear systems, random number generation, root finding, numerical integration).
   - Only elementary mathematical operations (`math.sqrt`, `math.sin`, `math.cos`, `math.log`, `math.exp`, `math.pi`, `math.ceil`) and basic plotting (`matplotlib.pyplot`) are permitted.
   - Standard data structures are native Python lists and nested lists.

2. **Single Central Library Rule (`mylib.py`):**
   - All numerical routines (LCG, LU decomposition, Cholesky, Jacobi, Gauss-Seidel, Newton-Raphson, Laguerre, Simpson, Gaussian quadrature, etc.) **must** live inside `mylib.py` at the root directory.
   - Front-end assignment scripts (`question1.py`, `question2.py`, etc.) must remain lightweight. They are strictly restricted to:
     - Reading input data / command-line arguments.
     - Invoking routines from `mylib.py`.
     - Outputting results to console or external text files.
     - Producing plots where requested.

3. **External File I/O (Non-Interactive Requirement):**
   - Inputs must **never** be hardcoded inside scripts or taken via interactive terminal prompts (`input()`).
   - All input data (matrices, vectors, parameters) must be read from external ASCII files (e.g., `asgn0_matA`, `asgn2_mat1`, `asgn3_vec1`).
   - Outputs must be redirected or saved to separate files (`output2.txt`, etc.) exactly mirroring program execution. Any explanatory physics comments must be appended below the raw output.

4. **Program Headers & Attribution:**
   - Every script must begin with a clear comment block containing:
     - Problem statement.
     - Name: Aryan Bandyopadhyay.
     - Roll Number: 2411014.

5. **Style & Academic Integrity:**
   - Code must read like that of a rigorous, mathematically sound physics student—clean, well-commented, direct, and unbloated.
   - Avoid robotic AI idioms, generic corporate boilerplate, and unnecessary object-oriented wrappers where simple, modular functions are expected.

---

## 2. Directory Structure & Assignments Map

```
D:\NISER\Computational Physics Lab\
│
├── mylib.py                 <-- Master library containing all algorithms
├── mylib.py.bak             <-- Safe backup of original library
├── COURSE_CONTEXT.md        <-- This architectural and context document
│
├── Assgn0\                  <-- Warmup: basic loops, series (GP/HP), complex numbers, matrix file I/O
│   ├── question1.py         (Sum of first 20 even numbers, factorial of 8)
│   ├── question2.py         (GP and HP series sums without analytical formula)
│   ├── question3.py         (Matrix multiplication from external files: AB, BC, D.C)
│   └── question5.py         (MyComplex class: addition, subtraction, multiplication, modulus)
│
├── assign1_rand\            <-- Pseudo-Random Number Generation (LCG, distributions, simulations)
│   ├── question1.py         (Non-linear chaotic map x_{i+1} = c*x_i*(1-x_i), correlation plots)
│   ├── question2.py         (LCG implementation, k-lag correlation check)
│   ├── question3.py         (Monte Carlo determination of pi, throwing method & residual analysis)
│   ├── question4.py         (Radioactive decay simulation A -> B -> C with kinetics overlay)
│   └── question5.py         (Inverse-transform sampling: exponential distribution exp(-x))
│
├── assign2_gj_lu\           <-- Direct Linear Solvers: Gauss-Jordan & LU (Doolittle)
│   ├── question1.py         (Library implementations with partial row pivoting)
│   ├── question2.py         (3-variable linear system solution via Gauss-Jordan)
│   └── question3.py         (LU decomposition verification: A = L * U)
│
├── assign3_lu_chsk\         <-- LU Forward-Backward & Cholesky Decomposition
│   ├── question1.py         (Cholesky factorization routine with symmetry check)
│   ├── question2.py         (Solving 6x6 linear system via LU forward-backward substitution)
│   └── question3.py         (Solving Ax = b using Cholesky factor L: L y = b, L^T x = y)
│
├── assign4_jac_gs\          <-- Iterative Solvers: Jacobi, Gauss-Seidel, and SOR
│   ├── question1.py         (Library routines with diagonal dominance row swapping)
│   ├── question2.py         (Solving 4x4 system: comparing Jacobi vs Gauss-Seidel convergence)
│   └── question3.py         (Tri-diagonal / sparse system: Gauss-Seidel vs SOR with omega = 1.57)
│
├── assign5_bi_rf\           <-- Root Finding: Bisection & Regula Falsi
│   ├── question1.py         (Library routines with interval auto-bracketing)
│   ├── question2.py         (Root of transcendental equation: ln(x/2) - sin(5x/2) = 0 in [1.5, 3.0])
│   └── question3.py         (Root of -x - cos(x) = 0, bracket searching from [2, 3])
│
├── assign5_fx_nr\           <-- Root Finding: Fixed-Point & Newton-Raphson (1D & Multivariate)
│   ├── question1.py         (Fixed-point and 1D Newton-Raphson library routines)
│   ├── question2.py         (Convergence comparison: Bisection vs Regula Falsi vs Newton-Raphson)
│   ├── question3.py         (Fixed-point iteration x = g(x) for x^2 - 2x - 3 = 0)
│   └── question4.py         (Multivariable non-linear system: Jacobian via central differences)
│
├── assign5_poly\            <-- Polynomial Roots: Horner, Laguerre, & Synthetic Division
│   ├── question1.py         (Horner value, 1st & 2nd derivatives, Laguerre iteration, deflation)
│   └── question2.py         (Extracting all real roots of quartic and quintic polynomials)
│
├── assign6_midtrap\         <-- Closed Newton-Cotes Integration: Midpoint & Trapezoidal
│   ├── question1.py         (Library implementations of Midpoint and Trapezoidal rules)
│   └── question2.py         (Integrating 1/x, x*cos(x), x*arctan(x) for N = 4, 8, 15, 20 with analytical checks)
│
├── assign6b_simmc\          <-- Simpson's 1/3-Rule & Uniform Monte Carlo Integration
│   ├── question1.py         (Simpson's rule with even-N validation, Monte Carlo with variance)
│   ├── question2.py         (Integration to 10^-6 precision with theoretical N from error bounds)
│   └── question3.py         (Monte Carlo integration of sin^2(x) on [-1, 1], tracking sigma / sqrt(N))
│
├── assign6c_gq\             <-- Gaussian Quadrature: Gauss-Legendre & Gauss-Laguerre
│   ├── question1.py         (Unified routine for orders N = 1 to 6 using hardwired weights and nodes)
│   ├── question2.py         (Gauss-Legendre 4-point vs Simpson 1/3 for x^2 / (1 + x^4))
│   ├── question3.py         (Gauss-Laguerre 5-point for integral_0^inf e^(-x)/(1 + x) dx)
│   └── question4.py         (Gauss-Laguerre integration of x^4 * e^(-x) comparing to 4! = 24)
│
└── Slides\                  <-- Lecture notes (PDFs from SPS NISER)
```

---

## 3. Mathematical & Algorithmic Blueprint

### 3.1 Pseudo-Random Numbers (`myrand`, LCG)
- **Recurrence:** $x_{i+1} = (a \cdot x_i + c) \pmod m$
- **Parameters (GCC Runtime / Slide Standard):** $a = 1103515245$, $c = 12345$, $m = 32768$ ($2^{15}$).
- **Transformation to uniform $[a, b)$:** $X = a + (b - a) \cdot \xi$, where $\xi = x_i / m \in [0, 1)$.
- **Transformation to exponential $\lambda e^{-\lambda y}$:** By inverse CDF: $y = -\frac{1}{\lambda} \ln(u)$, where $u \in (0, 1]$.

### 3.2 Direct Linear Solvers ($A x = b$)
- **Gauss-Jordan Elimination:** Augmented matrix $[A \mid b] \to [I \mid x]$. Uses partial row pivoting (choosing pivot row $p \ge i$ maximizing $|A_{pi}|$) to eliminate roundoff explosion and prevent division by zero.
- **Matrix Inversion:** Gauss-Jordan on $[A \mid I] \to [I \mid A^{-1}]$.
- **LU Decomposition (Doolittle):** $A = L \cdot U$, where $L_{ii} = 1$ (unit lower triangular) and $U$ is upper triangular.
  $$U_{ij} = A_{ij} - \sum_{k=0}^{i-1} L_{ik} U_{kj} \quad (j \ge i)$$
  $$L_{ji} = \frac{1}{U_{ii}} \left( A_{ji} - \sum_{k=0}^{i-1} L_{jk} U_{ki} \right) \quad (j > i)$$
- **Forward-Backward Substitution:**
  1. Solve $L y = b$ (forward: $y_0 = b_0$, $y_i = b_i - \sum_{j<i} L_{ij} y_j$).
  2. Solve $U x = y$ (backward: $x_{n-1} = y_{n-1}/U_{n-1,n-1}$, $x_i = \frac{1}{U_{ii}} (y_i - \sum_{j>i} U_{ij} x_j)$).
- **Free Determinant:** $\det(A) = \det(L) \det(U) = \prod_{i=0}^{n-1} U_{ii}$.

### 3.3 Cholesky Decomposition ($A = L L^T$)
- **Applicability:** Strictly symmetric ($A = A^T$) and positive-definite ($x^T A x > 0 \iff$ all leading principal minors $> 0$).
- **Algorithm:**
  $$L_{ii} = \sqrt{A_{ii} - \sum_{k=0}^{i-1} L_{ik}^2}$$
  $$L_{ji} = \frac{1}{L_{ii}} \left( A_{ji} - \sum_{k=0}^{i-1} L_{jk} L_{ik} \right) \quad (j > i)$$
- **Efficiency:** Twice as fast as LU ($\sim \frac{1}{3} n^3$ flops vs $\frac{2}{3} n^3$).
- **Forward-Backward Solver:** Solve $L y = b$, then $L^T x = y$.

### 3.4 Iterative Linear Solvers (Jacobi, Gauss-Seidel, SOR)
- **Convergence Requirement:** Guaranteed if $A$ is strictly diagonally dominant ($|A_{ii}| > \sum_{j \ne i} |A_{ij}|$) or symmetric positive-definite (for Gauss-Seidel).
- **Jacobi:** Updates evaluated strictly using prior step $x^{(k)}$:
  $$x_i^{(k+1)} = \frac{1}{A_{ii}} \left( b_i - \sum_{j \ne i} A_{ij} x_j^{(k)} \right)$$
- **Gauss-Seidel:** Uses newly computed values immediately in-place:
  $$x_i^{(k+1)} = \frac{1}{A_{ii}} \left( b_i - \sum_{j < i} A_{ij} x_j^{(k+1)} - \sum_{j > i} A_{ij} x_j^{(k)} \right)$$
- **Successive Over-Relaxation (SOR):**
  $$x_i^{(k+1)} = (1 - \omega) x_i^{(k)} + \frac{\omega}{A_{ii}} \left( b_i - \sum_{j < i} A_{ij} x_j^{(k+1)} - \sum_{j > i} A_{ij} x_j^{(k)} \right)$$
  - $\omega = 1$: Standard Gauss-Seidel.
  - $1 < \omega < 2$: Over-relaxation (accelerates slow convergence).
  - $0 < \omega < 1$: Under-relaxation (stabilizes divergence / oscillations).

### 3.5 Root Finding for Non-Linear Equations
- **Auto-Bracketing (`bracket_root`):**
  If $f(a) f(b) > 0$, expand outwards:
  - If $|f(a)| < |f(b)| \implies a \leftarrow a - \beta(b - a)$
  - If $|f(b)| < |f(a)| \implies b \leftarrow b + \beta(b - a)$
- **Bisection:** $c = (a + b)/2$. Retains subinterval with opposite sign. Linear convergence: error halves every iteration ($|b_n - a_n| = |b_0 - a_0| / 2^n$).
- **Regula Falsi (False Position):** Secant chord:
  $$c = \frac{a f(b) - b f(a)}{f(b) - f(a)} = b - \frac{f(b)(b - a)}{f(b) - f(a)}$$
- **Fixed-Point Iteration:** $x = g(x) \implies x_{n+1} = g(x_n)$. Converges iff $|g'(x^*)| < 1$.
- **Newton-Raphson (1D):**
  $$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$
  Derivative $f'(x)$ can be supplied analytically or computed numerically via central difference: $f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}$. Quadratic convergence near simple roots.
- **Multivariate Newton-Raphson System:** For $F(x) = [f_1(x), \dots, f_n(x)]^T = 0$:
  $$x^{(k+1)} = x^{(k)} - [J(x^{(k)})]^{-1} F(x^{(k)})$$
  where $J_{ij} = \partial f_i / \partial x_j$ is evaluated numerically via central differences and inverted via Gauss-Jordan elimination.

### 3.6 Polynomial Roots (Horner, Laguerre, Deflation)
- **Horner's Scheme (`polynomial_value`):** Evaluates $P(x) = a_n x^n + \dots + a_0$ in $O(n)$ additions and multiplications:
  $$b_n = a_n, \quad b_k = a_k + x \cdot b_{k+1}$$
- **Laguerre's Method:**
  $$G = \frac{P'(x)}{P(x)}, \quad H = G^2 - \frac{P''(x)}{P(x)}$$
  $$a = \frac{n}{G \pm \sqrt{(n - 1)(n H - G^2)}}$$
  Sign in denominator chosen to maximize $|G \pm \dots|$ (avoiding round-off cancellation). Update: $x_{new} = x - a$.
- **Deflation via Synthetic Division:** Divide $P(x)$ by $(x - r)$ to obtain $Q(x)$ of degree $n - 1$. Deflate iteratively to isolate all real roots.

### 3.7 Numerical Integration & Quadrature
- **Composite Midpoint Rule:**
  $$M_N = h \sum_{i=0}^{N-1} f(a + (i + 0.5)h), \quad h = \frac{b - a}{N}$$
  Theoretical error bound: $E_M \le \frac{(b - a)^3}{24 N^2} \max_{x \in [a, b]} |f''(x)|$.
- **Composite Trapezoidal Rule:**
  $$T_N = \frac{h}{2} \left[ f(a) + 2 \sum_{i=1}^{N-1} f(a + i h) + f(b) \right]$$
  Theoretical error bound: $E_T \le \frac{(b - a)^3}{12 N^2} \max_{x \in [a, b]} |f''(x)|$.
- **Composite Simpson's 1/3-Rule:**
  $$S_N = \frac{h}{3} \left[ f(a) + 4 \sum_{odd} f(x_i) + 2 \sum_{even} f(x_i) + f(b) \right]$$
  Requires $N$ to be an **even integer** (odd number of grid points).  
  Theoretical error bound: $E_S \le \frac{(b - a)^5}{180 N^4} \max_{x \in [a, b]} |f''''(x)|$.  
  Relation: $S_{2N} = \frac{2}{3} M_N + \frac{1}{3} T_N$.
- **Uniform Monte Carlo Integration:**
  $$F_N = (b - a) \frac{1}{N} \sum_{i=1}^N f(X_i), \quad \sigma_{integral} = \frac{(b - a) \sigma_f}{\sqrt{N}}$$
  Convergence rate scales as $O(N^{-1/2})$, independent of dimension.
- **Gaussian Quadrature:**
  - **Gauss-Legendre:** Over $[-1, 1]$, mapped to $[a, b]$ via $x = \frac{b - a}{2} t + \frac{b + a}{2}$:
    $$\int_a^b f(x) dx \approx \frac{b - a}{2} \sum_{i=1}^N w_i f\left( \frac{b - a}{2} t_i + \frac{b + a}{2} \right)$$
    Exact for all polynomials up to degree $2N - 1$.
  - **Gauss-Laguerre:** Over $[0, \infty)$ with weight function $e^{-t}$:
    $$\int_0^\infty e^{-t} f(t) dt \approx \sum_{i=1}^N w_i f(t_i)$$
    Note: When calling `gaussian_quadrature(f, 0, inf, N, method='laguerre')`, $f(t)$ represents the factor *without* $e^{-t}$.

---

## 4. Master Library Functions Reference (`mylib.py`)

| Function | Signature | Description |
|---|---|---|
| `myrand` | `(seed=0)` | LCG PRNG generator returning float in $[0, 1)$ |
| `random_uniform` | `(a=0.0, b=1.0)` | Uniform random float in $[a, b)$ |
| `random_exponential` | `(lam=1.0)` | Exponential random deviate $\sim \lambda e^{-\lambda x}$ |
| `MyComplex` | `(real, image=0.0)` | Complex number class with arithmetic methods & operator overloading |
| `read_matrix_from_file` | `(filename)` | Robust ASCII matrix loader (handles full or relative paths) |
| `write_matrix_to_file` | `(matrix, filename, precision=6)` | Writes 2D/1D list to external file in aligned columns |
| `print_matrix` | `(matrix, label=None, precision=4)` | Formatted column-aligned printer |
| `matrix_multiply` | `(X, Y)` | Standard matrix multiplication $X \cdot Y$ |
| `dot_product_vector` | `(X, Y)` | Vector dot product (supports 1D lists or column vectors) |
| `matrix_transpose` | `(A)` | Computes transpose $A^T$ |
| `matrix_add` / `matrix_sub` | `(A, B)` | Element-wise matrix addition / subtraction |
| `matrix_vector_multiply` | `(A, x)` | Matrix-vector product $A x$ |
| `vector_norm` | `(v, p=2)` | Vector $L_2$, $L_1$, or $L_\infty$ norm |
| `matrix_residual` | `(A, x, b)` | Verification helper: $\|A x - b\|_2$ |
| `is_symmetric` | `(A, tol=1e-9)` | Verification helper: checks if $A = A^T$ |
| `is_diagonally_dominant` | `(A, strict=False)` | Verification helper: checks diagonal dominance |
| `matrix_determinant_lu` | `(A)` | Evaluates $\det(A) = \prod U_{ii}$ via Doolittle LU |
| `gauss_jordan_elimination_augmented`| `(augmented_matrix)` | Solves $[A \mid b]$ with partial row pivoting |
| `gauss_jordan_inverse` | `(A)` | Inverts matrix via $[A \mid I] \to [I \mid A^{-1}]$ |
| `lu_decomposition` | `(A)` | Doolittle LU factorization ($L_{ii} = 1$) |
| `lu_forback` | `(A, b, method='lu')` | Solves $A x = b$ via LU/Cholesky forward-backward substitution |
| `cholesky_decomposition`| `(A)` | Computes lower factor $L$ such that $A = L L^T$ |
| `cholesky_forback` | `(L, b)` | Solves $L y = b$, then $L^T x = y$ |
| `jacobi_it` | `(A, b, tol=1e-8, max_iterations=1000)` | Jacobi iterative solver with diagonal dominance pivoting |
| `gauss_seidel` | `(A, b, tol=1e-8, max_iterations=500)` | Gauss-Seidel solver with in-place sequential updates |
| `sor_gauss_seidel` | `(A, b, tol=1e-8, max_iterations=500, omega=1.0)` | Successive Over-Relaxation solver |
| `bracket_root` | `(f, a, b, beta=0.1, max_iter=100)` | Expands interval until $f(a) f(b) < 0$ |
| `bisection` | `(f, a, b, accuracy=1e-8, max_iterations=100)` | Bisection root finder with auto-bracketing |
| `regula_falsi` | `(f, a, b, accuracy=1e-8, max_iterations=100)` | False position root finder with auto-bracketing |
| `fixed_point` | `(g, x0, accuracy=1e-6, max_iterations=30)` | Picard iteration $x_{n+1} = g(x_n)$ |
| `finite_difference_derivative` | `(f, x, h=1e-5, order=1)` | Central difference 1st ($order=1$) or 2nd ($order=2$) derivative |
| `partial_derivative` | `(f, var_index, point, h=1e-5)` | Multivariable partial derivative $\partial f / \partial x_i$ |
| `jacobian` | `(f, point, h=1e-6)` | Numerical Jacobian matrix for vector functions |
| `newton_raphson` | `(f, df=None, x0=0.0, accuracy=1e-6, max_iterations=30, h=1e-5)` | 1D Newton-Raphson (analytical or numerical derivative) |
| `newton_raphson_system` | `(f, J, initial_guess, accuracy=1e-6, max_iterations=30)` | Multivariate non-linear system solver via Newton-Raphson |
| `polynomial_value` | `(coefficients, x)` | Horner's method for polynomial evaluation |
| `polynomial_first_derivative` | `(coefficients)` | Analytical coefficients of $P'(x)$ |
| `polynomial_second_derivative`| `(coefficients)` | Analytical coefficients of $P''(x)$ |
| `laguerre` | `(coefficients, b0, accuracy1=1e-8, accuracy2=1e-6, max_iterations=30)` | Single real root finder via Laguerre formula |
| `synthetic_division` | `(coefficients, root, accuracy=1e-6)` | Polynomial deflation by $(x - r)$ |
| `laguerre_roots` | `(coefficients, b0=0.0)` | Finds all real roots via Laguerre + deflation |
| `midpoint` | `(f, a, b, N)` | Composite midpoint numerical integration |
| `trapezoidal` | `(f, a, b, N)` | Composite trapezoidal numerical integration |
| `simpson` | `(f, a, b, N)` | Simpson's 1/3-rule (returns `(integral, evaluations)`) |
| `monte_carlo` | `(f, a, b, N, seed=1)` | Uniform Monte Carlo (returns `(integral, sigma_f, sigma_int)`) |
| `integration_error_bound_N` | `(method, a, b, max_deriv, target_error)` | Evaluates theoretical minimum $N$ from error bounds |
| `gaussian_quadrature` | `(f, a, b, N, method='legendre')` | Unified Gauss-Legendre and Gauss-Laguerre quadrature ($N \in \{1..6\}$) |

---

## 5. Standard Code Template for Assignments

When generating new scripts for class assignments, adhere strictly to this pattern:

```python
# Problem: [Exact statement from assignment sheet]
# Author: Aryan Bandyopadhyay | Roll Number: 2411014
# Computational Physics Lab (PHY341 / PHY745), NISER

import math
import sys
import os

# Import required routines exclusively from central library
from mylib import (
    read_matrix_from_file,
    write_matrix_to_file,
    lu_forback,
    matrix_residual,
    print_matrix
)

if __name__ == '__main__':
    # 1. Read input parameters / matrices non-interactively
    matrix_file = 'asgn3_mat1'
    vector_file = 'asgn3_vec1'
    
    A = read_matrix_from_file(matrix_file)
    b = read_matrix_from_file(vector_file)
    
    # 2. Invoke library routine
    x_sol = lu_forback(A, b, method='lu')
    
    # 3. Verify solution integrity via residual
    res = matrix_residual(A, x_sol, b)
    
    # 4. Display formatted results
    print("Solution vector x:")
    for i, val in enumerate(x_sol):
        print(f"  x[{i}] = {val:.6f}")
    print(f"\nVerification Residual ||Ax - b|| = {res:.4e}")
    
    # 5. Save output to file as required by lab protocol
    output_filename = 'output2.txt'
    with open(output_filename, 'w') as f:
        f.write("# Computed Solution Vector:\n")
        for i, val in enumerate(x_sol):
            f.write(f"x[{i}] = {val:.6f}\n")
        f.write(f"\n# Verification Residual: {res:.4e}\n")
```

---

## 6. Post-Midsem Roadmap

Topics scheduled for the second half of the semester:
1. **Ordinary Differential Equations (ODEs):**
   - Forward Euler, Backward Euler, Predictor-Corrector (Heun / Milne).
   - Classical 4th-Order Runge-Kutta (RK4) for coupled systems.
   - Initial Value Problems (IVP) and Boundary Value Problems (Shooting Method, Finite Difference).
2. **Least Squares Curve Fitting:**
   - Linear regression, polynomial regression via normal equations.
   - Non-linear least squares (Gauss-Newton, Levenberg-Marquardt).
3. **Partial Differential Equations (PDEs) & Eigenvalues (if time permits):**
   - Heat/Diffusion equation (FTCS, Crank-Nicolson).
   - Wave equation.
   - Power method and QR algorithm for matrix eigenvalues.
