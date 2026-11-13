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
├── Midsem\                  <-- Mid-Semester Examination Suite
│   ├── question1.py         (Ehrenfest particle diffusion simulation between two compartments)
│   ├── question2.py         (Matrix inversion via Doolittle LU decomposition with det check)
│   ├── question3.py         (Wien's displacement law root finding via Newton-Raphson)
│   ├── question4.py         (Center of mass of non-uniform beam via composite midpoint rule)
│   └── question5.py         (Ground state hydrogen <r> expectation value via Gauss-Legendre)
│
├── assign7a_fwbk\           <-- ODEs: Forward Euler, Backward Euler & Predictor-Corrector
│   └── question1.py         (Comparative IVP solutions and step-size stability study)
│
├── assign7b_rk4\            <-- ODEs: Runge-Kutta 4th Order (RK4) Engine
│   ├── q1.py                (Core RK4 library function verification for 1st/2nd order systems)
│   ├── q2.py                (Solving dy/dx = (x+y)^2 against analytical tan(x + pi/4) - x)
│   ├── q3.py                (Damped harmonic oscillator: trajectory, phase portrait, energy decay)
│   └── q4.py                (3D chaotic Rössler attractor phase space evolution)
│
├── assign7c_sht\            <-- Two-Point Boundary Value Problems & FEM
│   ├── question1.py         (Cantilever beam depression via Shooting Method with RK4)
│   ├── question2.py         (Uninsulated rod heat conduction with convective ambient cooling)
│   └── question3.py         (1D Poisson equation BVP via Finite Element Method & Cholesky)
│
├── assign7d_felem\          <-- Partial Differential Equations (PDEs)
│   ├── question1.py         (1D Heat equation via explicit FTCS scheme with stability analysis)
│   └── question2.py         (Generalized parabolic convection-diffusion via Crank-Nicolson)
│
├── assign8_eval\            <-- Eigenvalue Problems: Power Iteration
│   ├── question1.py         (Core Power Iteration library routine with Rayleigh quotient)
│   ├── question2.py         (Dominant eigenvalue/eigenvector of 4x4 matrix)
│   └── question3.py         (First two eigenvalues of 5x5 tridiagonal matrix via deflation)
│
├── assign9a_fit\            <-- Interpolation & Comparative Data Fitting
│   ├── question1.py         (Lagrange polynomial interpolation library function)
│   ├── question2.py         (Weighted linear regression with Pearson r^2)
│   ├── question3.py         (Lagrange interpolation at x = 6.7 with ascending polynomial coefficients)
│   └── question4.py         (Comparative fitting: Power law vs Exponential model via Pearson r^2)
│
├── assign9b_fit\            <-- Polynomial & Chebyshev Least Squares Regression
│   ├── question1.py         (Monomial polynomial fit of order k with Hankel condition number)
│   ├── question2.py         (Shifted Chebyshev basis regression on [0, 1])
│   ├── question3.py         (Cubic least squares fit on experimental data assign9fit.txt)
│   └── question4.py         (Condition number and stability comparison: Monomial vs Chebyshev)
│
├── Slides\                  <-- Complete collection of course lecture slide decks (PDFs)
└── helper\                  <-- Comprehensive study notes (Computational_Physics_Complete_Notes.md)
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

### 3.8 Ordinary Differential Equations (ODEs)
- **Initial Value Problems (IVPs):** $\frac{dy}{dx} = f(x, y)$ with $y(x_0) = y_0$.
- **Explicit Forward Euler:**
  $$y_{n+1} = y_n + h f(x_n, y_n)$$
  Local truncation error $O(h^2)$, global error $O(h)$. Simple but conditionally stable ($|1 + h\lambda| \le 1$).
- **Implicit Backward Euler:**
  $$y_{n+1} = y_n + h f(x_{n+1}, y_{n+1})$$
  Solved via Newton-Raphson: $g(y) = y - y_n - h f(x_{n+1}, y) = 0$. Unconditionally stable for stiff equations ($\text{Re}(\lambda) < 0$).
- **Heun's Predictor-Corrector Method:**
  $$y_{n+1}^p = y_n + h f(x_n, y_n) \quad \text{(Predictor)}$$
  $$y_{n+1} = y_n + \frac{h}{2} \left[ f(x_n, y_n) + f(x_{n+1}, y_{n+1}^p) \right] \quad \text{(Corrector)}$$
  Global accuracy $O(h^2)$.
- **Runge-Kutta 4th-Order (RK4):**
  $$k_1 = f(x_n, y_n)$$
  $$k_2 = f\left(x_n + \frac{h}{2}, y_n + \frac{h}{2} k_1\right)$$
  $$k_3 = f\left(x_n + \frac{h}{2}, y_n + \frac{h}{2} k_2\right)$$
  $$k_4 = f(x_n + h, y_n + h k_3)$$
  $$y_{n+1} = y_n + \frac{h}{6} (k_1 + 2k_2 + 2k_3 + k_4)$$
  Fourth-order accuracy $O(h^4)$, derived from Simpson's 1/3 quadrature. Extends naturally to coupled systems $\mathbf{y}' = \mathbf{f}(x, \mathbf{y})$:
  - Second-order ODEs $y'' = g(x, y, y')$ decomposed into $\mathbf{y} = [y, v]^T$, $\mathbf{y}' = [v, g(x, y, v)]^T$.
  - Applied to damped harmonic oscillators ($E(t)$ tracking) and 3D chaotic attractors (Rössler system).
- **Two-Point Boundary Value Problems (BVPs):**
  - **Shooting Method:** Converts BVP $y'' = f(x, y, y'), y(a) = \alpha, y(b) = \beta$ into IVP with guessed initial slope $z_0 = y'(a)$. Shoots via RK4 to $x = b$. Uses linear interpolation between two trial slopes $z_1, z_2$ to find the exact initial condition $z^*$:
    $$z^* = z_1 + (z_2 - z_1) \frac{\beta - y(b; z_1)}{y(b; z_2) - y(b; z_1)}$$
  - **1D Finite Element Method (FEM):** Weak formulation $\int (y' v' + f v) dx = 0$ with piecewise linear hat basis functions. Assembles a symmetric positive-definite tridiagonal stiffness matrix, solved via Cholesky decomposition.

### 3.9 Partial Differential Equations (PDEs)
- **1D Heat / Diffusion Equation:** $\frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2}$ on $[0, L] \times [0, T]$ with parameter $\alpha = \frac{\Delta t}{\Delta x^2}$.
  - **Explicit FTCS:** $u_i^{j+1} = \alpha (u_{i+1}^j + u_{i-1}^j) + (1 - 2\alpha) u_i^j$. Stable only for $\alpha \le 0.5$.
  - **Implicit BTCS:** $-\alpha u_{i-1}^{j+1} + (1 + 2\alpha) u_i^{j+1} - \alpha u_{i+1}^{j+1} = u_i^j$. Unconditionally stable; tridiagonal system solved at each step via Thomas Algorithm (TDMA) in $O(N_x)$ time.
  - **Crank-Nicolson Scheme:** Averages explicit and implicit spatial differences:
    $$(1 + \alpha) u_i^{j+1} - \frac{\alpha}{2} (u_{i+1}^{j+1} + u_{i-1}^{j+1}) = (1 - \alpha) u_i^j + \frac{\alpha}{2} (u_{i+1}^j + u_{i-1}^j)$$
    Second-order accurate in both space and time: $O(\Delta t^2 + \Delta x^2)$. Unconditionally stable.
- **1D Wave Equation (Hyperbolic):** $\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$ with Courant number $r = \frac{c \Delta t}{\Delta x}$.
  - Discretization: $u_i^{j+1} = 2(1 - r^2) u_i^j + r^2(u_{i+1}^j + u_{i-1}^j) - u_i^{j-1}$.
  - Courant-Friedrichs-Lewy (CFL) stability criterion: $r \le 1.0$.
- **2D Elliptic Equations (Laplace / Poisson):** $\nabla^2 u = \rho(x, y)$.
  - 5-point discrete Laplacian stencil with Successive Over-Relaxation (SOR):
    $$u_{i,j}^{(k+1)} = (1 - \omega) u_{i,j}^{(k)} + \frac{\omega}{4} \left( u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - h^2 \rho_{i,j} \right)$$
    Tuning relaxation parameter $\omega \in (1, 2)$ accelerates spectral convergence.

### 3.10 Eigenvalue Problems (Power Iteration & Deflation)
- **Problem Formulation:** $A \mathbf{x} = \lambda \mathbf{x}$ for $n \times n$ matrix $A$.
- **Power Iteration Method:**
  Repeated matrix-vector multiplication $\mathbf{x}_{k+1} = A \mathbf{x}_k$. If eigenvalues satisfy $|\lambda_1| > |\lambda_2| \ge \dots \ge |\lambda_n|$:
  $$\mathbf{x}_k = c_1 \lambda_1^k \mathbf{v}_1 + \sum_{i=2}^n c_i \lambda_i^k \mathbf{v}_i \approx c_1 \lambda_1^k \mathbf{v}_1$$
  Rayleigh quotient for dominant eigenvalue:
  $$\lambda_1^{(k)} = \frac{A \mathbf{x}_k \cdot \mathbf{y}}{\mathbf{x}_k \cdot \mathbf{y}} \to \lambda_1$$
  Normalized at every iteration ($\mathbf{x}_k \leftarrow \mathbf{x}_k / \|\mathbf{x}_k\|_2$). Convergence rate scales as $|\lambda_2 / \lambda_1|^k$.
- **Matrix Deflation for Subdominant Eigenvalues:**
  For symmetric $A$, eliminate the dominant eigenspace via outer-product subtraction:
  $$A' = A - \lambda_1 \mathbf{u}_1 \mathbf{u}_1^T, \quad \mathbf{u}_1 = \frac{\mathbf{v}_1}{\|\mathbf{v}_1\|}$$
  Power iteration on $A'$ yields $\lambda_2$ and $\mathbf{v}_2$.
- **Analytical Tridiagonal Spectrum:**
  For a symmetric tridiagonal Toeplitz matrix with subdiagonal $a$, main diagonal $b$, superdiagonal $c$:
  $$\lambda_k = b + 2\sqrt{ac} \cos\left(\frac{k\pi}{n+1}\right), \quad v_{k,i} = \sin\left(\frac{ik\pi}{n+1}\right), \quad k=1,\dots,n$$

### 3.11 Interpolation & Least Squares Data Fitting
- **Lagrange Interpolation Formula:**
  $$P_N(x) = \sum_{i=0}^N y_i \prod_{k \ne i} \frac{x - x_k}{x_i - x_k}$$
  Polynomial coefficients extracted in ascending order $[c_0, c_1, \dots, c_N]$ using polynomial basis expansion.
- **Weighted Linear Regression ($y = m x + c$):**
  Minimizes $\chi^2 = \sum_i w_i (y_i - m x_i - c)^2$ with weights $w_i = 1 / \sigma_i^2$.
  $$S = \sum w_i, \quad S_x = \sum w_i x_i, \quad S_y = \sum w_i y_i, \quad S_{xx} = \sum w_i x_i^2, \quad S_{xy} = \sum w_i x_i y_i, \quad S_{yy} = \sum w_i y_i^2$$
  $$\Delta = S S_{xx} - S_x^2, \quad m = \frac{S S_{xy} - S_x S_y}{\Delta}, \quad c = \frac{S_{xx} S_y - S_x S_{xy}}{\Delta}$$
  Pearson's correlation coefficient:
  $$r^2 = \frac{(S S_{xy} - S_x S_y)^2}{\Delta (S S_{yy} - S_y^2)}$$
- **Polynomial Least Squares via Normal Equations:**
  Fitting $y(x) = \sum_{j=0}^k a_j x^j$ of order $k$ solves $(A^T W A) \mathbf{a} = A^T W \mathbf{y}$, where $A_{ij} = x_i^j$.
  The normal matrix $M = A^T W A$ is a Hankel matrix: $M_{jl} = \sum_i w_i x_i^{j+l}$. Solved via Doolittle LU decomposition.
- **Condition Number & Shifted Chebyshev Basis:**
  Monomial powers $x^j$ produce severely ill-conditioned Hankel matrices as degree $k$ increases ($\kappa(M) = \|M\|_F \|M^{-1}\|_F \gg 10^4$).
  Replacing monomials with orthogonal polynomials stabilizes the regression. On $[0, 1]$, shifted Chebyshev polynomials $\phi_m(x) = T_m(2x - 1)$ satisfy:
  $$\phi_0(x) = 1, \quad \phi_1(x) = 2x - 1, \quad \phi_2(x) = 8x^2 - 8x + 1, \quad \phi_3(x) = 32x^3 - 48x^2 + 18x - 1$$
  Recurrence relation for arbitrary $m \ge 1$:
  $$\phi_{m+1}(x) = 2(2x - 1) \phi_m(x) - \phi_{m-1}(x)$$
  Provides orders-of-magnitude reduction in condition number and eliminates numerical instability.

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
| `euler` | `(f, a, b, y0, h, method='forward')` | Explicit Forward or implicit Backward Euler for IVPs |
| `predictor_corrector` | `(f, a, b, y0, h)` | Heun's predictor-corrector 2nd-order method |
| `rk4` | `(f, a, b, y0, h=None, n=None)` | 4th-Order Runge-Kutta for 1st-order and coupled/2nd-order ODEs |
| `lagrange_linear_interpolation` | `(x1, y1, x2, y2, y_target)` | Linear interpolation for boundary condition matching |
| `fem_1d_poisson` | `(f_rhs, a, b, ya, yb, n_elements)` | 1D Finite Element Method with linear hat functions & Cholesky |
| `thomas_algorithm` | `(a, b, c, d)` | $O(n)$ Tridiagonal Matrix Algorithm (TDMA) |
| `pde_heat_1d_explicit` | `(L, T, nx, nt, u0_func, ...)` | Explicit FTCS scheme for 1D diffusion ($\alpha \le 0.5$) |
| `pde_heat_1d_implicit` | `(L, T, nx, nt, u0_func, ...)` | Implicit BTCS scheme via Thomas algorithm (unconditionally stable) |
| `pde_heat_1d_crank_nicolson` | `(L, T, nx, nt, u0_func, ...)` | Crank-Nicolson $O(\Delta t^2 + \Delta x^2)$ scheme |
| `pde_wave_1d_explicit` | `(L, T, nx, nt, c, u0_func, ...)` | Explicit 2nd-order wave equation solver ($r = c \Delta t / \Delta x \le 1$) |
| `pde_laplace_poisson_2d`| `(Lx, Ly, nx, ny, rho_func, ...)` | 2D elliptic boundary value solver via SOR relaxation |
| `power_iteration` | `(A, tolerance=1e-8, max_iterations=1000)` | Dominant eigenvalue and eigenvector via Rayleigh quotient |
| `lagrange_interpolation` | `(X, x)` | Lagrange polynomial interpolation with ascending power coefficients |
| `linear_regression` | `(X, sigma=None)` | Weighted linear regression with Pearson $r^2$ goodness-of-fit |
| `poly_eval` | `(coefficients, x)` | Horner's scheme for polynomial evaluation |
| `poly_fit` | `(X, k, sigma=None)` | Monomial polynomial least squares of degree $k$ with condition number |
| `modified_chebyshev_basis` | `(m, x)` | Shifted Chebyshev polynomial $\phi_m(x) = T_m(2x - 1)$ on $[0, 1]$ |
| `chebyshev_fit` | `(X, k=3, sigma=None)` | Least-squares fitting using shifted Chebyshev orthogonal basis |

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

## 6. Post-Midsem Status

All scheduled topics have been implemented, verified, and integrated into `mylib.py`:
1. **Ordinary Differential Equations (ODEs) — Complete (`assign7a_fwbk`, `assign7b_rk4`, `assign7c_sht`):**
   - Explicit Forward Euler, Implicit Backward Euler (Newton solver), Heun's Predictor-Corrector.
   - Classical 4th-Order Runge-Kutta (RK4) for coupled 1st- and 2nd-order ODEs (Damped SHO, Rössler attractor).
   - Boundary Value Problems via Shooting Method with RK4 and Lagrange interpolation.
   - 1D Finite Element Method (FEM) with piecewise linear hat functions and Cholesky decomposition.
2. **Partial Differential Equations (PDEs) — Complete (`assign7d_felem`):**
   - 1D Diffusion/Heat equation via explicit FTCS and unconditionally stable Crank-Nicolson schemes.
   - 1D hyperbolic Wave equation with Courant-Friedrichs-Lewy (CFL) stability enforcement.
   - 2D elliptic Laplace/Poisson solver via Successive Over-Relaxation (SOR).
3. **Eigenvalues & Eigenvectors — Complete (`assign8_eval`):**
   - Power Iteration method with Rayleigh quotient.
   - Matrix deflation via outer products for subdominant eigenvalues.
   - Analytical verification on tridiagonal Toeplitz spectrum.
4. **Interpolation & Least Squares Curve Fitting — Complete (`assign9a_fit`, `assign9b_fit`):**
   - Lagrange polynomial interpolation with ascending power coefficients.
   - Weighted linear regression with Pearson $r^2$ applied to power-law vs. exponential model comparisons.
   - Standard polynomial least squares of arbitrary degree $k$ with Frobenius condition number of the Hankel matrix.
   - Shifted Chebyshev orthogonal polynomial fitting $\phi_m(x) = T_m(2x - 1)$ on $[0, 1]$.
