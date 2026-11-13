# PHY341/745 — Physics Computer Lab: Complete Master Notes

> Compiled from **all 10 slide decks**: `intro`, `random`, `gaussj`, `ludecomp`, `itermat`, `root`, `integ`, `diffeq`, `eigv`, `datafit`.
> Everything your prof taught — every method, equation, algorithm, worked example and exercise — laid out in one place.

**Contents**

- [0. Course Map, Syllabus & Programming Rules](#0-course-map-syllabus--programming-rules)
- [1. Random Numbers](#1-random-numbers)
- [2. Linear Equations I — Gauss–Jordan Elimination](#2-linear-equations-i--gaussjordan-elimination)
- [3. Linear Equations II — LU Decomposition & Cholesky](#3-linear-equations-ii--lu-decomposition--cholesky)
- [4. Linear Equations III — Iterative Methods: Jacobi, Gauss–Seidel, SOR](#4-linear-equations-iii--iterative-methods-jacobi-gaussseidel-sor)
- [5. Roots of Nonlinear Equations](#5-roots-of-nonlinear-equations)
- [6. Numerical Integration](#6-numerical-integration)
- [7. Ordinary Differential Equations (ODEs)](#7-ordinary-differential-equations-odes)
- [8. Partial Differential Equations (PDEs)](#8-partial-differential-equations-pdes)
- [9. Eigenvalues & Eigenvectors](#9-eigenvalues--eigenvectors)
- [10. Interpolation & Data Fitting](#10-interpolation--data-fitting)
- [Appendix A — Formula Cheat Sheet](#appendix-a--formula-cheat-sheet)
- [Appendix B — All Exercises / Tryouts in One List](#appendix-b--all-exercises--tryouts-in-one-list)
- [Appendix C — Assignments (Repo): Questions & How They're Solved](#appendix-c--assignments-repo-questions--how-theyre-solved)

---

## 0. Course Map, Syllabus & Programming Rules

### 0.1 Syllabus (from intro slides)

1. Random number generation
2. Gauss–Jordan elimination — solving linear equations, inverse, determinant
3. LU decomposition — solving linear equations, inverse
4. Jacobi / Gauss–Seidel — iterative solution of the above
5. Root finding — Regula falsi, Newton–Raphson, Laguerre (+ Bisection, Fixed-point)
6. Integration — Midpoint, Trapezoidal, Simpson, Monte Carlo (+ Gaussian Quadrature)
7. ODE — Euler, Predictor–Corrector, Runge–Kutta (IVP + BVP) *(not in these slides)*
8. Least square fitting *(not in these slides)*
9. If time permits: PDE & Eigenvalue problem; ML regression and classification

### 0.2 Programming practices (prof's rules — follow for every assignment/DIY/exam)

1. Every code has a **comment header**: purpose, name, roll number.
2. **Input is NEVER interactive or hardwired.** Read everything (dimensions, precision, filenames, numbers) from **external file(s)** / command-line args.
3. All output goes to **separate file(s)**, exactly as generated. Without output the problem is incomplete. Extra comments appended after output, clearly separated.
4. Decent amount of comments explaining computing steps (not too many, not too few).
5. **Initialization in loops** where needed (reading matrix, coefficients, …).
6. Every routine you write (`myrand`, `lu_decomp`, `rk4`, …) lives in **one single growing library file**, NOT in the main front code.
7. Main/front code does only I/O, memory allocation, file opening, and **calls** to library routines.
8. **No system-built routines for major algorithms** — no `numpy`, `scipy`, `random`, `inverse`, etc. for the algorithm itself. Basic `sqrt()`, `sin()`, `log()` etc. are fine.
9. Endsem: full syllabus, must use **your own library functions**.

### 0.3 Git basics (from intro slides)

```bash
# create + init
mkdir "Computational physics Examples" && cd "Computational physics Examples"
git init

# work, stage, commit
echo "an example of git workflow" >> chapter1.txt
git status
git add chapter1.txt
git commit -m "create chapter1"
git log

# second file
echo "creating chapter2" >> chapter2.txt
git add chapter2.txt
git commit -m "add chapter2"

# link remote + push
git remote add origin git@github.com:<user>/<repo>.git
git remote -v
git push -u origin main

# branch + pull
git switch -c newbranch
git pull origin main
```

Core idea: **working dir → staging (`git add`) → local repo (`git commit -m`) → remote (`git push`)**.

---

## 1. Random Numbers

### 1.1 Why random numbers? True vs Pseudo RNG

Random numbers are everywhere in physics: thermodynamics, radioactivity, particle collisions, everything in between.

| Type | Idea | Pros | Cons |
|---|---|---|---|
| **True RNG** | Natural phenomenon: coin flip, dice, radioactive decay, thermal noise, atmospheric radio-noise | Truly random | Needs post-processing, **slow**, not useful for regular use |
| **Pseudo RNG (pRNG)** | Deterministic algorithm, generated iteratively | Extremely **fast, portable** | Deterministic, **finite sequence length**, correlated |

> **Basic goal of the course:** write your **own pRNG** and use it for all assignments, exams and DIY.

Sequence length can be made *veryyy* long by proper choice of parameters.

### 1.2 A quick-and-dirty pRNG (logistic map) — why we need tests

$$x_{i+1} = c\,x_i(1-x_i)$$

- $x_0$ is the **seed** — it defines the whole random sequence.
- Exercise in slides: $x_0 = 0.1$ and $c = 3.2,\,3.7,\,3.98$.
- Plot $x_i$ vs $x_{i+5}$ to *see* correlatedness.
- Bad pRNG for many choices of $c$ and seed $x_0$ — settles into regular pattern. For $c=3.98$, $x_0=0.1$ no specific pattern (looks better).
- Lesson: **eyes are good at discerning patterns but can fool us too** → need mathematical tests.

### 1.3 Quantifying good vs bad pRNG

Need mathematical tests: **if pRNG fails a test, don't use it.**

- Basic tests: **correlation, moments**
- Advanced tests: **chi-square, Kolmogorov–Smirnov**
- Ideal randoms: **no correlations**, error purely statistical, i.e. scales as $1/\sqrt{N}$.

**Correlations test:**

$$\epsilon(n,N) = \frac{1}{N}\sum_{i=1}^{N} x_i x_{i+n} - \left(\frac{1}{N}\sum_{i=n}^{N} x_i\right)^2$$

**Connected correlations test:**

$$\epsilon(n,N) = \frac{1}{N}\sum_{i=1}^{N} x_i x_{i+n} - \frac{1}{N}\sum_{i=1}^{N}x_i \cdot \frac{1}{N}\sum_{i=1}^{N}x_{i+n}$$

> If tuplets of random numbers are **not correlated**, $\epsilon(n,N) \to 0$ with statistical error $1/\sqrt{N}$.

### 1.4 Linear Congruential Generator (LCG) — the workhorse

One of the oldest and most common uniform pRNGs:

$$x_{i+1} = (a\,x_i + c)\ \mathrm{mod}\ m \equiv \mathrm{remainder}\left[\frac{a x_i + c}{m}\right]$$

- $m$ determines the **period**: produces integers $0,\dots,m-1$. Then $x_i/m$ gives randoms in $[0,1)$.
- $m$ typically chosen $2^{32}$.
- $a$ = multiplier, $0 < a < m$. Numerical Recipes: $a=1664525$; gcc: $a=1103515245$.
- $c$ = increment, $0 < c < m$. Numerical Recipes: $c=1013904223$; gcc: $c=12345$.

Not all rosy — bad parameter choices collapse immediately:

- $a,c,m,x_0 = 6,7,5,2:\quad 4,1,2,0,2,4,1,2,0,2,\dots$
- $a,c,m,x_0 = 27,11,54,2:\quad 11,38,11,38,\dots$

### 1.5 Hull–Dobell theorem (when does LCG have full period $m$?)

LCG is extremely sensitive to $a,c,x_0,m$. In particular $a$ must be chosen with great care or short/very-short periodicity sets in.

> **Hull–Dobell:** LCG has period $m$ iff $c \ne 0$ **and**:
> 1. $c$ is **coprime** to $m$,
> 2. $a-1$ is a **multiple of $p$** for every prime $p$ dividing $m$,
> 3. $a-1$ is a **multiple of 4** if $m$ is a multiple of 4.

- Works well when $m$ has many repeated prime factors $p$, e.g. a **power of 2**.
- If $m$ is **square-free** (no $n^2$ factor), only $a=1$ is allowed → very bad pRNG.
- LCG is extremely fast, least memory footprint, but period severely limited by $m$: for $m \sim 2^{32} \to 10^9$ randoms — **gets exhausted in seconds!**

**$c=0$ case → Lehmer / Park–Miller pRNG:**

$$x_{i+1} = a\,x_i \ \mathrm{mod}\ m$$

- $m$ can be a prime, or a prime just less than a power of 2 (**Mersenne primes** $2^{31}-1$, $2^{61}-1$, …), or a simple power of 2.

### 1.6 LCG exercises from slides

1. $a=27,\ c=11,\ m=54,\ x_0=10$ (**Ugly**): `0.204, 0.704, 0.204, 0.704, …`
2. $a=57,\ c=1,\ m=256,\ x_0=10$ (**Bad**)
3. $a=1103515245,\ c=12345,\ m=32768,\ x_0=10$ (**Good**)

Try them, histogram them, and run the correlation test on each.

### 1.7 Application 1: Determination of $\pi$

- Unit circle $r=1$ at origin → area $=\pi$.
- Enclose in square of side $2r=2$.
- Confine to first quadrant, randomly generate $(x,y)$, check **inside**: $x^2+y^2 \le 1$ → area $=\pi/4$.
- Over LARGE trials:

$$\frac{\text{circle area}}{\text{square area}} \approx \frac{\text{total inside}}{\text{total trials}} \;\;\Longrightarrow\;\; \pi \approx 4 \times \frac{\text{total inside}}{\text{total trials}}$$

### 1.8 Application 2: Radioactivity

Radioactive decay is intrinsically random at atomic level — can't predict when one atom decays. But a large collection decays at predictable rate $\propto \lambda$.

Decay probability over $\Delta t$:

$$P = \lambda\,\Delta t$$

To simulate $A \to B$ or $A \to B \to C$:

1. Choose random $r \in [0,1]$.
2. If $r < P_A$: $A\to B$, i.e. $N_A \to N_A-1$, $N_B \to N_B+1$.
3. Cycle over all $N_A$ atoms in the time interval.
4. For $A\to B\to C$: repeat over all $N_B$ generated so far, comparing $r < P_B$ and updating $N_B$, $N_C$ accordingly.

### 1.9 Non-uniform distributions from uniform randoms

Standard pRNG gives uniform integers $\ell \in [0,\mathrm{INT_{MAX}}]$ or floats $x = \ell/\mathrm{INT_{MAX}} \in [0,1)$.

**Uniform in $[a,b)$ from $x\in[0,1)$:**

$$u = a + (b-a)\,x$$

**General transformation law:** if $p(x)$ is pdf of uniform $x$ and target pdf is $q(y)$:

$$|q(y)\,dy| = |p(x)\,dx| \;\;\Longrightarrow\;\; q(y) = p(x)\left|\frac{dx}{dy}\right|$$

For uniform $x\in[0,1)$, $p(x)=1$.

**Exponentially distributed RN:**

$$q(y) = a\,e^{-ay},\quad y\ge 0,\ a>0$$

From transformation law:

$$a\,e^{-ay} = \left|\frac{dx}{dy}\right| \;\Longrightarrow\; x = \int_0^y q(y')\,dy' = 1 - e^{-ay}$$

$$\boxed{\,y = -\frac{1}{a}\ln(1-x) \equiv -\frac{1}{a}\ln x\,}$$

because $(1-x)\in[0,1)$ as well whenever $x\in[0,1)$.

### 1.10 How to histogram / call randoms — and what YOU should do

Prof's point: **don't call system RNG in your final work — write and call your own `myrand` seeded only once.**

Histogram in Python (from slides, cleaned):

```python
import matplotlib.pyplot as plt
from mylib import *   # your own library file

if __name__ == "__main__":
    myexprn = []
    myrn = my_exprand(1234)   # seed the pRNG once
    for i in range(2000):     # 2K randoms for histogram
        myrn = my_exprand()
        myexprn.append(myrn)
    plt.hist(myexprn, bins=40, color='skyblue', edgecolor='black')
    plt.show()
```

System RNG calls (for reference only — **don't use for graded work**):

```c
// C — system RNG (reference only)
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char *argv[]) {
    srand(1234); // seed only once
    for (int i = 0; i < 20; i++) {
        float myrand = (float)rand() / (float)RAND_MAX;
        printf("Random # %d : %.4f\n", i, myrand);
    }
    return 0;
}
```

```python
# Python — system RNG (reference only)
import random
random.seed(12345)  # seed only once
for i in range(20):
    print(random.uniform(0, 1))
```

**What you should do — your own LCG in C:**

```c
#include <stdio.h>

float myrand(int seed) {
    static int x0, x1;
    static int a = 1103515245;
    static int c = 12345;   // + m etc.
    static int m = 1 << 31; // example: 2^31
    float rn;
    if (seed) x0 = seed;
    x1 = (a * x0 + c) % m;
    rn = (float)x1 / (float)m;  // convert to [0.0, 1.0]
    x0 = x1;
    return rn;
}

int main(int argc, char *argv[]) {
    float myrn = myrand(1234);  // seed once
    for (int i = 0; i < 20; i++)
        printf("Random # %d : %.4f\n", i, myrand(0));
    return 0;
}
```

**What you should do — your own LCG in Python:**

```python
def myrand(seed=0):
    if not hasattr(myrand, "x0"):
        myrand.x0 = 0
    a = 1103515245
    c = 12345
    m = 2**31
    if seed:
        myrand.x0 = seed
    x1 = (a * myrand.x0 + c) % m
    rn = x1 / m  # convert to [0.0, 1.0]
    myrand.x0 = x1
    return rn

if __name__ == "__main__":
    myrn = myrand(1234)  # seed once
    for i in range(20):
        print(f"Random # {i} : {myrand():.4f}")
```

---

## 2. Linear Equations I — Gauss–Jordan Elimination

### 2.1 The problem $A\cdot x = b$

Commonly encountered:

$$A \cdot x = b$$

- $A$ is $n\times n$ matrix; $x,b$ are $n\times 1$ vectors (linear system) or $n\times n$ matrices (inverse).
- Formal solution:

$$x = A^{-1}b,\qquad A^{-1} = A_c^T/|A|$$

where $|A|$ = determinant, $A_c$ = cofactors. Numerically doable but far more efficient methods exist.

| Family | Methods | Best for |
|---|---|---|
| **Direct elimination** | Gauss–Jordan, LU, Cholesky | Small, dense matrices; memory-intensive, possibly slower, hard to parallelize |
| **Iterative** | Gauss–Jacobi, Gauss–Seidel, Conjugate Gradient, GMRES | Large, esp. sparse matrices; lend themselves to parallelization |

### 2.2 Augmented matrix + RREF idea

For 3 variables:

$$A\cdot X = B \;\Longrightarrow\; \begin{pmatrix} a_{11}&a_{12}&a_{13}\\ a_{21}&a_{22}&a_{23}\\ a_{31}&a_{32}&a_{33} \end{pmatrix}\!\cdot\!\begin{pmatrix}x_1\\x_2\\x_3\end{pmatrix} = \begin{pmatrix}b_1\\b_2\\b_3\end{pmatrix}$$

Write the **augmented matrix**:

$$[A\,|\,B] \equiv \left[\begin{array}{ccc|c} a_{11}&a_{12}&a_{13}&b_1\\ a_{21}&a_{22}&a_{23}&b_2\\ a_{31}&a_{32}&a_{33}&b_3 \end{array}\right]$$

Gauss–Jordan elimination converts it to **reduced row echelon form (RREF)**:

$$\left[\begin{array}{ccc|c} 1&0&0&\tilde b_1\\ 0&1&0&\tilde b_2\\ 0&0&1&\tilde b_3 \end{array}\right] \;\Longrightarrow\; x_1=\tilde b_1,\ x_2=\tilde b_2,\ x_3=\tilde b_3.$$

### 2.3 RREF: definition + elementary row operations

In RREF:

1. All rows with only zero entries are at the bottom.
2. First nonzero entry (**pivot**) of each nonzero row is to the right of the leading entry of the row above.
3. Leading entry (pivot) in any nonzero row is **1**.
4. All other entries in the row or column containing a leading 1 are **zeros**.

Illustration:

$$\begin{pmatrix} 1&-2&3\\ 0&1&-4\\ 0&0&1 \end{pmatrix}$$

Achieve RREF by combining three **elementary row operations**:

- swapping two rows,
- multiplying a row by a nonzero number,
- adding/subtracting a multiple of one row to another.

> Operation involving **rows only** = **partial pivoting**.
> **Full pivoting** (recommended by *Numerical Recipes*): both row and column operations.

### 2.4 Partial pivoting — step by step

Say the augmented matrix is:

$$\left[\begin{array}{ccc|c} 0&a_{12}&a_{13}&b_1\\ a_{21}&a_{22}&a_{23}&b_2\\ a_{31}&a_{32}&a_{33}&b_3 \end{array}\right]$$

$a_{11}$ is $0$ or order(s) of magnitude smaller than $a_{21},a_{31}$ → $a_{31},a_{21} \gg a_{11}$, let $a_{31}>a_{21}$. To get leftmost nonzero entry at top, swap $R_3 \leftrightarrow R_1$:

$$\left[\begin{array}{ccc|c} a_{31}&a_{32}&a_{33}&b_3\\ a_{21}&a_{22}&a_{23}&b_2\\ 0&a_{12}&a_{13}&b_1 \end{array}\right] \equiv \left[\begin{array}{ccc|c} a^0_{11}&a^0_{12}&a^0_{13}&b^0_1\\ a^0_{21}&a^0_{22}&a^0_{23}&b^0_2\\ a^0_{31}&a^0_{32}&a^0_{33}&b^0_3 \end{array}\right]$$

where $a^0_{11}$ = **pivot element**, new $R^0_1$ = **pivot row**, new $C^0_1$ = **pivot column**. Normalize $R^0_1/a^0_{11} \to R^1_1$:

$$\left[\begin{array}{ccc|c} 1&a^1_{12}&a^1_{13}&b^1_1\\ a^0_{21}&a^0_{22}&a^0_{23}&b^0_2\\ a^0_{31}&a^0_{32}&a^0_{33}&b^0_3 \end{array}\right],\quad a^1_{1j}=a^0_{1j}/a^0_{11},\ b^1_1=b^0_1/a^0_{11}.$$

Reduce pivot column $C_1$ by adding/subtracting multiples of pivot row from following rows:

$$\begin{aligned}
a^1_{2j} &= a^0_{2j} - a^0_{21}\,R^1_1, & b^1_2 &= b^0_2 - a^0_{21}\,b^1_1\\
a^1_{3j} &= a^0_{3j} - a^0_{31}\,R^1_1, & b^1_3 &= b^0_3 - a^0_{31}\,b^1_1
\end{aligned}$$

$$\Longrightarrow \left[\begin{array}{ccc|c} 1&a^1_{12}&a^1_{13}&b^1_1\\ 0&a^1_{22}&a^1_{23}&b^1_2\\ 0&a^1_{32}&a^1_{33}&b^1_3 \end{array}\right]$$

- $a^1_{22}$ is the **new pivot**, undergoes same smallness test; if needed, partial pivot again.
- Eventually RREF:

$$\left[\begin{array}{ccc|c} 1&0&0&b^3_1\\ 0&1&0&b^3_2\\ 0&0&1&b^3_3 \end{array}\right]$$

> This is **Gauss–Jordan elimination with partial pivoting**. No matter what steps and in which order, the final augmented matrix is **unique**.

### 2.5 Pivoting tryout (from slides — do it!)

$$\begin{aligned}
y + z - 2w &= -3\\
x + 2y - z &= 2\\
2x + 4y + z - 3w &= -2\\
x - 4y - 7z - w &= -19
\end{aligned}$$

Write matrices in a file, construct augmented matrix, **partial-pivot the first row only**, convert first column to $(1,0,0,0)^T$ and write down the augmented matrix.

> N.B. While pivoting, swap rows **below** the new pivot row.

### 2.6 Full worked example

$$\left.\begin{aligned}
2y + 5z &= 1\\
3x - y + 2z &= -2\\
x - y + 3z &= 3
\end{aligned}\right\} \Longrightarrow \left[\begin{array}{ccc|c} 0&2&5&1\\ 3&-1&2&-2\\ 1&-1&3&3 \end{array}\right]$$

Since $a_{11}=0$, $a_{21}>a_{31}$, swap $R_1\leftrightarrow R_2$:

$$\left[\begin{array}{ccc|c} 3&-1&2&-2\\ 0&2&5&1\\ 1&-1&3&3 \end{array}\right] \xrightarrow{R^0_1/3\,\to\,R^1_1} \left[\begin{array}{ccc|c} 1&-1/3&2/3&-2/3\\ 0&2&5&1\\ 1&-1&3&3 \end{array}\right]$$

$$\xrightarrow{R^0_3-R^1_1\,\to\,R^1_3} \left[\begin{array}{ccc|c} 1&-1/3&2/3&-2/3\\ 0&2&5&1\\ 0&-2/3&7/3&11/3 \end{array}\right] \xrightarrow{R^1_2/2\,\to\,R^2_2} \left[\begin{array}{ccc|c} 1&-1/3&2/3&-2/3\\ 0&1&5/2&1/2\\ 0&-2/3&7/3&11/3 \end{array}\right]$$

$$\xrightarrow{R^1_1+R^2_2/3\,\to\,R^2_1} \left[\begin{array}{ccc|c} 1&0&9/6&-3/6\\ 0&1&5/2&1/2\\ 0&-2/3&7/3&11/3 \end{array}\right] \xrightarrow{R^1_3+2R^2_2/3\,\to\,R^2_3} \left[\begin{array}{ccc|c} 1&0&9/6&-3/6\\ 0&1&5/2&1/2\\ 0&0&12/3&24/6 \end{array}\right]$$

$$\xrightarrow{3R^2_3/12\,\to\,R^3_3} \left[\begin{array}{ccc|c} 1&0&9/6&-3/6\\ 0&1&5/2&1/2\\ 0&0&1&1 \end{array}\right] \xrightarrow{R^2_2-5R^3_3/2\,\to\,R^3_2} \left[\begin{array}{ccc|c} 1&0&9/6&-3/6\\ 0&1&0&-2\\ 0&0&1&1 \end{array}\right]$$

$$\xrightarrow{R^2_1-9R^3_3/6\,\to\,R^3_1} \boxed{\left[\begin{array}{ccc|c} 1&0&0&-2\\ 0&1&0&-2\\ 0&0&1&1 \end{array}\right]} \;\Longrightarrow\; x=-2,\ y=-2,\ z=1.$$

### 2.7 Determinant from Gauss–Jordan (invertibility test)

$$\det\!\begin{pmatrix}a_{11}&a_{12}&a_{13}\\a_{21}&a_{22}&a_{23}\\a_{31}&a_{32}&a_{33}\end{pmatrix} = \det\!\begin{pmatrix}a'_{11}&a'_{12}&a'_{13}\\0&a'_{22}&a'_{23}\\0&0&a'_{33}\end{pmatrix} = (-1)^n\,a'_{11}a'_{22}a'_{33}$$

where $n$ = number of row interchanges, $a'_{ii}\ne 0$, and none of $a'_{ii}$ has to be 1.

### 2.8 Matrix inversion with Gauss–Jordan + cost/accuracy

$$A\cdot B = I \;\to\; \text{augmented }[A\,|\,I] \xrightarrow{\ G\text{-}J\ } [I\,|\,A^{-1}]$$

- **Memory:** $N\times(N+1)$ for linear solve; $N\times(N+N)$ for inverse (ignoring swap/sort).
- **Ops:** $(N-1)\times 2$ (mult/div + add/sub) for each row reduction.
- **Accuracy:** rounding errors $\propto N$ build up quickly.
- For large $N$ not suitable — memory $\propto N^2$, slow. Typically restricted to $\sim N=10$.

---

## 3. Linear Equations II — LU Decomposition & Cholesky

### 3.1 The idea: $A = L\cdot U$

$$A\cdot x = b \;\;\to\;\; A = L\cdot U$$

Factorize $A$ into **lower-triangular $L$** and **upper-triangular $U$**:

$$L=\begin{pmatrix}\ell_{00}&0&0\\\ell_{10}&\ell_{11}&0\\\ell_{20}&\ell_{21}&\ell_{22}\end{pmatrix},\qquad U=\begin{pmatrix}u_{00}&u_{01}&u_{02}\\0&u_{11}&u_{12}\\0&0&u_{22}\end{pmatrix}$$

Looks formidable but very useful and not nearly as hard. Multiplying:

$$\begin{pmatrix}a_{00}&a_{01}&a_{02}\\a_{10}&a_{11}&a_{12}\\a_{20}&a_{21}&a_{22}\end{pmatrix} = \begin{pmatrix} \ell_{00}u_{00}&\ell_{00}u_{01}&\ell_{00}u_{02}\\ \ell_{10}u_{00}&\ell_{10}u_{01}+\ell_{11}u_{11}&\ell_{10}u_{02}+\ell_{11}u_{12}\\ \ell_{20}u_{00}&\ell_{20}u_{01}+\ell_{21}u_{11}&\ell_{20}u_{02}+\ell_{21}u_{12}+\ell_{22}u_{22} \end{pmatrix}$$

> But **without proper ordering, factorization may fail!**

### 3.2 The $9$-equations vs $12$-unknowns catch: Doolittle vs Crout

$L\cdot U$ gives:

$$\begin{aligned}
\ell_{00}u_{00}&=a_{00}, &\ell_{00}u_{01}&=a_{01}, &\ell_{00}u_{02}&=a_{02}\\
\ell_{10}u_{00}&=a_{10}, &\ell_{10}u_{01}+\ell_{11}u_{11}&=a_{11}, &\ell_{10}u_{02}+\ell_{11}u_{12}&=a_{12}\\
\ell_{20}u_{00}&=a_{20}, &\ell_{20}u_{01}+\ell_{21}u_{11}&=a_{21}, &\ell_{20}u_{02}+\ell_{21}u_{12}+\ell_{22}u_{22}&=a_{22}
\end{aligned}$$

Catch: $3\times 3 = 9$ equations but $3\times(3+1)=12$ variables!

Trick: fix the diagonal of one factor to 1:

- **Doolittle:** $\ell_{00}=\ell_{11}=\ell_{22}=1$
- **Crout:** $u_{00}=u_{11}=u_{22}=1$

Straight away:

- Doolittle: $u_{00}=a_{00},\ u_{01}=a_{01},\ u_{02}=a_{02}$
- Crout: $\ell_{00}=a_{00},\ \ell_{10}=a_{10},\ \ell_{20}=a_{20}$

Rest of $\ell_{ij}/u_{ij}$ solved from remaining equations.

### 3.3 Doolittle LU — explicit order + generic form

Remember $u_{0j}=a_{0j}$. Do calculations **in the order they appear**:

$$\begin{array}{ll}
u_{00}=a_{00} & \ell_{00}=1\\
u_{01}=a_{01} & \ell_{01}=0\\
u_{02}=a_{02} & \ell_{02}=0\\
u_{10}=0, & \ell_{10}=a_{10}/u_{00}\\
u_{20}=0, & \ell_{20}=a_{20}/u_{00}\\
u_{11}=a_{11}-\ell_{10}u_{01}, & \ell_{11}=1\\
u_{21}=0, & \ell_{21}=(a_{21}-\ell_{20}u_{01})/u_{11}\\
u_{12}=a_{12}-\ell_{10}u_{02}, & \ell_{12}=0\\
u_{22}=a_{22}-\ell_{20}u_{02}-\ell_{21}u_{12}, & \ell_{22}=1
\end{array}$$

**Generic form** for each $j=0,1,\dots,N-1$, in order, for $i=1,2,\dots,N-1$:

$$\boxed{u_{ij} = a_{ij} - \sum_{k=0}^{i-1}\ell_{ik}u_{kj}\quad \text{for } i\le j}$$

$$\boxed{\ell_{ij} = \left(a_{ij} - \sum_{k=0}^{j-1}\ell_{ik}u_{kj}\right)\Big/u_{jj}\quad \text{for } i>j}$$

### 3.4 LU storage trick — half the memory

Code far more efficiently using $a_{ij}$ wherever $u_{ij}$ or $\ell_{ij}$ appears. Why?

> **Every $a_{ij}$ is used only once and never again** → $u_{ij},\ell_{ij}$ can be stored in the same memory as $a_{ij}$. Memory requirement is **half**.

Stored Doolittle LU matrix (overwrites $A$):

$$LU = \begin{pmatrix} u_{00}&u_{01}&u_{02}\\ \ell_{10}&u_{11}&u_{12}\\ \ell_{20}&\ell_{21}&u_{22} \end{pmatrix} \to A$$

- No need to store $\ell_{ii}=1$; modify loops accordingly.
- Numerical cost same as Gauss–Jordan: $O(N^3)$.
- But **can LU always be done?** (next section)

### 3.5 Can we always LU-decompose?

- If $a_{00}=0$, then either $L$ or $U$ is singular → impossible if $A$ is not. **Solution: row pivot.**
- **Guaranteed if all leading submatrices have nonzero determinant.** Example:

$$A=\begin{pmatrix}1&0&2\\-2&1&0\\1&0&1\end{pmatrix}\Rightarrow A_1=1,\ A_2=\begin{pmatrix}1&0\\-2&1\end{pmatrix},\ A_3=A$$

$$\det A_1=1,\ \det A_2=1,\ \det A_3=-3 \;\checkmark$$

- **Additional pivoting** if any leading submatrix has zero determinant but matrix itself is invertible:

$$A=\begin{pmatrix}1&0&2\\3&0&2\\-1&1&2\end{pmatrix}\Rightarrow \det A_1=1,\ \det A_2=0,\ \det A_3=4$$

- Otherwise no solution exists and you are doomed!

### 3.6 $U$ + backward substitution (recall from Gauss–Jordan)

$$\begin{pmatrix}a_{00}&a_{01}&a_{02}\\a_{10}&a_{11}&a_{12}\\a_{20}&a_{21}&a_{22}\end{pmatrix}\!\begin{pmatrix}x_0\\x_1\\x_2\end{pmatrix}=\begin{pmatrix}b_0\\b_1\\b_2\end{pmatrix} \Rightarrow \begin{pmatrix}u_{00}&u_{01}&u_{02}\\0&u_{11}&u_{12}\\0&0&u_{22}\end{pmatrix}\!\begin{pmatrix}x_0\\x_1\\x_2\end{pmatrix}=\begin{pmatrix}\bar b_0\\\bar b_1\\\bar b_2\end{pmatrix}$$

Solutions from **last row upward (backward substitution)**:

$$\begin{aligned}
u_{22}x_2&=\bar b_2 &\Rightarrow x_2&=\bar b_2/u_{22}\\
u_{11}x_1+u_{12}x_2&=\bar b_1 &\Rightarrow x_1&=(\bar b_1-u_{12}x_2)/u_{11}\\
u_{00}x_0+u_{01}x_1+u_{02}x_2&=\bar b_0 &\Rightarrow x_0&=(\bar b_0-u_{01}x_1-u_{02}x_2)/u_{00}
\end{aligned}$$

Generic $N\times N$ backward substitution:

$$\boxed{x_i = \frac{1}{u_{ii}}\left(\bar b_i - \sum_{j=i+1}^{N}u_{ij}x_j\right),\quad x_N=\bar b_N/u_{NN},\ i=N-2,\dots,0}$$

### 3.7 Solving with LU: forward–backward substitution

Advisable to **begin with partial pivoting**. Then split:

$$A\cdot x = b \;\Longrightarrow\; L\cdot\underbrace{(U\cdot x)}_{y} = b,\qquad U\cdot x = y \;\Longrightarrow\; L\cdot y = b$$

1. **Forward substitution:** solve $L\cdot y = b$ for $y$.
2. **Backward substitution:** solve $U\cdot x = y$ for $x$.

**Forward:**

$$\begin{pmatrix}1&0&0\\\ell_{10}&1&0\\\ell_{20}&\ell_{21}&1\end{pmatrix}\!\begin{pmatrix}y_0\\y_1\\y_2\end{pmatrix}=\begin{pmatrix}b_0\\b_1\\b_2\end{pmatrix}$$

$$\begin{aligned}
y_0 &= b_0\\
\ell_{10}y_0+y_1 &= b_1 \;\Rightarrow\; y_1=b_1-\ell_{10}y_0\\
\ell_{20}y_0+\ell_{21}y_1+y_2 &= b_2 \;\Rightarrow\; y_2=b_2-\ell_{20}y_0-\ell_{21}y_1
\end{aligned}$$

**Backward:**

$$\begin{pmatrix}u_{00}&u_{01}&u_{02}\\0&u_{11}&u_{12}\\0&0&u_{22}\end{pmatrix}\!\begin{pmatrix}x_0\\x_1\\x_2\end{pmatrix}=\begin{pmatrix}y_0\\y_1\\y_2\end{pmatrix}$$

$$\begin{aligned}
x_2 &= y_2/u_{22}\\
x_1 &= (y_1-u_{12}x_2)/u_{11}\\
x_0 &= (y_0-u_{01}x_1-u_{02}x_2)/u_{00}
\end{aligned}$$

**Generic:**

$$\boxed{y_i = b_i - \sum_{j=0}^{i-1}\ell_{ij}y_j,\quad y_0=b_0,\ i=1,\dots,N-1}$$

$$\boxed{x_i = \frac{1}{u_{ii}}\left(y_i-\sum_{j=i+1}^{N}u_{ij}x_j\right),\quad x_N=y_N/u_{NN},\ i=N-2,\dots,0}$$

- **Determinant for free:** $\det A = \det(LU)=\det L\times\det U = (-1)^n\prod_i u_{ii}$ ($n$ = row swaps).
- **Inverse:** iterate through each column of the identity matrix.

### 3.8 Cholesky decomposition ($A = LL^T$)

Symmetries of $A$ often simplify solving. **Cholesky:** factorization of **Hermitian, positive-definite** matrix (common in physics, e.g. covariance matrix) into $L$ and $L^T$.

For $3\times3$ (real):

$$A = LL^\dagger \xrightarrow{\ \mathrm{real}\ } LL^T \;\Longrightarrow\; \begin{pmatrix}a_{00}&a_{01}&a_{02}\\a_{10}&a_{11}&a_{12}\\a_{20}&a_{21}&a_{22}\end{pmatrix}=\begin{pmatrix}\ell_{00}&0&0\\\ell_{10}&\ell_{11}&0\\\ell_{20}&\ell_{21}&\ell_{22}\end{pmatrix}\!\begin{pmatrix}\ell_{00}&\ell_{10}&\ell_{20}\\0&\ell_{11}&\ell_{21}\\0&0&\ell_{22}\end{pmatrix}$$

When $A$ is real and positive-definite:

$$\boxed{\ell_{ii} = \pm\sqrt{a_{ii}-\sum_{j=0}^{i-1}\ell_{ij}^2}}$$

$$\boxed{\ell_{ij} = \frac{1}{\ell_{ii}}\left(a_{ij}-\sum_{k=0}^{i-1}\ell_{ik}\ell_{jk}\right)\quad \text{for } i>j}$$

- Signs before square roots are inconsequential. DIY the complex-matrix version.
- About **TWICE as efficient** as LU for solving linear systems.
- Wikipedia example:

$$\begin{pmatrix}4&12&-16\\12&37&-43\\-16&-43&98\end{pmatrix}=\begin{pmatrix}2&0&0\\6&1&0\\-8&5&3\end{pmatrix}\!\begin{pmatrix}2&6&-8\\0&1&5\\0&0&3\end{pmatrix}$$

- Also used in: nonlinear optimization (multivariable), Monte Carlo (covariance matrix), inversion of Hermitian matrices.

---

## 4. Linear Equations III — Iterative Methods: Jacobi, Gauss–Seidel, SOR

### 4.1 Why iterative?

Gauss–Jordan and LU quickly become unmanageable → **pivoting costly** with increasing size + increasing zeros (**sparse matrix**). Best strategy: **iterative methods: Jacobi & Gauss–Seidel**.

### 4.2 Jacobi method

Works when $A$ has **no zeros on the main diagonal**. Decomposes (not factorizes) $A$ into diagonal $D$ plus $L,U$:

$$A = D + (L+U) = \begin{pmatrix}a_{00}&0&\cdots&0\\0&a_{11}&\cdots&0\\\vdots&&\ddots&\vdots\\0&0&\cdots&a_{n-1,n-1}\end{pmatrix} + \begin{pmatrix}0&a_{01}&\cdots&a_{0,n-1}\\a_{10}&0&\cdots&a_{1,n-1}\\\vdots&&\ddots&\vdots\\a_{n-1,0}&a_{n-1,1}&\cdots&0\end{pmatrix}$$

Then:

$$[D+(L+U)]\cdot x = b \;\Longrightarrow\; D\cdot x = b-(L+U)\cdot x$$

$$\boxed{x^{(k+1)} = D^{-1}\big(b-(L+U)\,x^{(k)}\big)}$$

Component-wise, iterate until $||x^{k+1}-x^{k}||<\epsilon$:

$$\boxed{x^{(k+1)}_i = \frac{1}{a_{ii}}\left(b_i - \sum_{j\ne i}a_{ij}\,x^{(k)}_j\right)}$$

- Fairly fast, but computing $x^{(k+1)}_i$ needs **all** other elements → **cannot overwrite** $x^{(k)}$ with $x^{(k+1)}$ (needs two vectors).
- **Sufficient condition for convergence:** $A$ **strictly diagonally dominant**:

$$\boxed{|a_{ii}| \gtrsim \sum_{j\ne i}|a_{ij}|}$$

- Relatively slower, but **easily parallelized** → handles large matrices, reduces memory + time.

**Worked example:**

$$\begin{aligned}
4x_1-x_2+x_3&=7\\-2x_1+x_2+5x_3&=15\\4x_1-8x_2+x_3&=-21
\end{aligned} \;\Longrightarrow\; A=\begin{pmatrix}4&-1&1\\-2&1&5\\4&-8&1\end{pmatrix} \to \begin{pmatrix}4&-1&1\\4&-8&1\\-2&1&5\end{pmatrix}$$

$A$ is not diagonally dominant, but swapping eqn 2 ↔ 3 (rows 2↔3) makes it so. Then:

$$x_1=\frac{7+x_2-x_3}{4},\quad x_2=\frac{21+4x_1+x_3}{8},\quad x_3=\frac{15+2x_1-x_2}{5}$$

$$\begin{aligned}
x^{(k+1)}_1&=(7+x^{(k)}_2-x^{(k)}_3)/4\\
x^{(k+1)}_2&=(21+4x^{(k)}_1+x^{(k)}_3)/8\\
x^{(k+1)}_3&=(15+2x^{(k)}_1-x^{(k)}_2)/5
\end{aligned}$$

Start $(0,0,0)$:

| $k$ | $x^{(k)}_1$ | $x^{(k)}_2$ | $x^{(k)}_3$ |
|---|---|---|---|
| 0 | 0 | 0 | 0 |
| 1 | 1.75 | 2.625 | 3.0 |
| 2 | 1.656 | 3.875 | 3.175 |
| … | … | … | … |
| 6 | 1.9995 | 3.996 | 2.9987 |

> Try **without** making diagonally dominant: $2x+3y+z=13,\ x+y+z=6,\ 3x+2y+2z=15$ → $(x,y,z)=(3,2,1)$.

### 4.3 Gauss–Seidel method

$$A = L^\star + U \to \begin{pmatrix}a_{00}&0&\cdots&0\\a_{10}&a_{11}&\cdots&0\\\vdots&&\ddots&\vdots\\a_{n-1,0}&\cdots&\cdots&a_{n-1,n-1}\end{pmatrix} + \begin{pmatrix}0&a_{01}&\cdots&a_{0,n-1}\\0&0&\cdots&a_{1,n-1}\\\vdots&&\ddots&\vdots\\0&0&\cdots&0\end{pmatrix}$$

where $L^\star$ = lower-triangular **including diagonal**, $U$ = strict upper. Iterate $A\cdot x=b$ as:

$$\boxed{L^\star x^{(k+1)} = b - U\,x^{(k)}}$$

until $||x^{(k+1)}-x^{(k)}||<\epsilon$:

$$\boxed{x^{(k+1)}_i = \frac{1}{a_{ii}}\left(b_i - \sum_{j=0}^{i-1}a_{ij}x^{(k+1)}_j - \sum_{j=i+1}^{n-1}a_{ij}x^{(k)}_j\right),\quad i=0,\dots,n-1}$$

- Needs **only one storage vector** (vs two in Jacobi): $x^{(k+1)}$ uses already-computed elements + not-yet-updated $x^{(k)}$.
- Convergence guaranteed only if matrix is **strictly diagonally dominant OR symmetric and positive-definite**.

**Explicit $4\times4$ walkthrough (from slides):**

$$RHS = \begin{pmatrix}b_0\\b_1\\b_2\\b_3\end{pmatrix}-\begin{pmatrix}0&a_{01}&a_{02}&a_{03}\\0&0&a_{12}&a_{13}\\0&0&0&a_{23}\\0&0&0&0\end{pmatrix}\!\begin{pmatrix}x^{(k)}_0\\x^{(k)}_1\\x^{(k)}_2\\x^{(k)}_3\end{pmatrix} = \begin{pmatrix}b_0-a_{01}x^{(k)}_1-a_{02}x^{(k)}_2-a_{03}x^{(k)}_3\\b_1-a_{12}x^{(k)}_2-a_{13}x^{(k)}_3\\b_2-a_{23}x^{(k)}_3\\b_3\end{pmatrix} = b_i-\sum_{j=i+1}^{n-1}a_{ij}x^{(k)}_j$$

$$LHS = \begin{pmatrix}a_{00}&0&0&0\\a_{10}&a_{11}&0&0\\a_{20}&a_{21}&a_{22}&0\\a_{30}&a_{31}&a_{32}&a_{33}\end{pmatrix}\!\begin{pmatrix}x^{(k+1)}_0\\x^{(k+1)}_1\\x^{(k+1)}_2\\x^{(k+1)}_3\end{pmatrix} = a_{ii}x^{(k+1)}_i+\sum_{j=0}^{i-1}a_{ij}x^{(k+1)}_j$$

Equating $LHS=RHS$:

$$\begin{aligned}
a_{00}x^{(k+1)}_0 &= b_0-a_{01}x^{(k)}_1-a_{02}x^{(k)}_2-a_{03}x^{(k)}_3\\
a_{10}x^{(k+1)}_0+a_{11}x^{(k+1)}_1 &= b_1-a_{12}x^{(k)}_2-a_{13}x^{(k)}_3\\
a_{20}x^{(k+1)}_0+a_{21}x^{(k+1)}_1+a_{22}x^{(k+1)}_2 &= b_2-a_{23}x^{(k)}_3\\
a_{30}x^{(k+1)}_0+a_{31}x^{(k+1)}_1+a_{32}x^{(k+1)}_2+a_{33}x^{(k+1)}_3 &= b_3
\end{aligned}$$

i.e.

$$\begin{aligned}
x^{(k+1)}_0&=\tfrac{1}{a_{00}}(b_0-a_{01}x^{(k)}_1-a_{02}x^{(k)}_2-a_{03}x^{(k)}_3)\\
x^{(k+1)}_1&=\tfrac{1}{a_{11}}(b_1-a_{10}x^{(k+1)}_0-a_{12}x^{(k)}_2-a_{13}x^{(k)}_3)\\
x^{(k+1)}_2&=\tfrac{1}{a_{22}}(b_2-a_{20}x^{(k+1)}_0-a_{21}x^{(k+1)}_1-a_{23}x^{(k)}_3)\\
x^{(k+1)}_3&=\tfrac{1}{a_{33}}(b_3-a_{30}x^{(k+1)}_0-a_{31}x^{(k+1)}_1-a_{32}x^{(k+1)}_2)
\end{aligned}$$

> Note $x^{(k)}_i$ is used **only once** → can be replaced in-place with $x^{(k+1)}_i$.

### 4.4 Successive Over-Relaxation (SOR) — accelerated Gauss–Seidel

Introduces weight $\omega$ to extrapolate: weighted average of $x^{(k)}_i$ and newly computed Gauss–Seidel $x^{(k+1)}_i$:

$$\boxed{x^{(k+1)}_i = (1-\omega)\,x^{(k)}_i + \frac{\omega}{a_{ii}}\left(b_i-\sum_{j=0}^{i-1}a_{ij}x^{(k+1)}_j-\sum_{j=i+1}^{n-1}a_{ij}x^{(k)}_j\right)}$$

- $\omega=1$: normal Gauss–Seidel
- $1<\omega<2$: **over-relaxation** for slow convergence
- $0<\omega<1$: **under-relaxation** for oscillating/diverging solution

Example (tridiagonal $2$ on diagonal, $-1$ off-diagonal, RHS $=1$):

$$\begin{pmatrix}2&-1&0&\cdots&0\\-1&2&-1&\cdots&0\\0&-1&2&\cdots&0\\\vdots&&&&\vdots\\0&0&0&\cdots&2\end{pmatrix}\!\begin{pmatrix}x_1\\\vdots\\x_n\end{pmatrix}=\begin{pmatrix}1\\\vdots\\1\end{pmatrix}$$

Gauss–Seidel ≈ 150 iterations, SOR ≈ 25 iterations with $\omega=1.57$ for $n\ge 3$.

**Exercise from slides — solve:**

$$\begin{aligned}
4x_1+x_2-x_3+x_4&=-2\\x_1+4x_2-x_3-x_4&=-1\\-x_1-x_2+5x_3+x_4&=0\\x_1-x_2+x_3+3x_4&=1
\end{aligned}$$

---

## 5. Roots of Nonlinear Equations

### 5.1 Setup, bracketing, the five algorithms

Solve in one variable:

$$\boxed{f(x)=0}$$

- Assume $f(x)$ smoothly varying, continuous over $[a,b]$.
- Typical examples: $f(x)=\cos x-x^3,\ 3x+\sin x-e^x,\ xe^x-2,\ x^3+3x-5$, etc.
- If $x_0\in[a_0,b_0]$ satisfies $f(x_0)=0$, $x_0$ is a **root/zero**.
- If $f(a_0)$, $f(b_0)$ have **opposite signs**, by the **intermediate value theorem** $f$ has at least one root in $[a_0,b_0]$.

Finding roots starts with a guess $[a_0,b_0]$ (informed or trial-and-error) where $f$ has opposite signs → **$[a_0,b_0]$ brackets the root**. Iterations produce shrinking intervals $[a_0,b_0]\to[a_i,b_i]$ always containing one root. Good initial guess via (i) plotting $f(x)$ vs $x$ or (ii) informed expectation.

| # | Method | Needs bracket? | Notes |
|---|---|---|---|
| 1 | **Bisection** | Yes | Simplest, slowest, guaranteed |
| 2 | **False position (Regula falsi)** | Yes | Interpolation, faster than bisection |
| 3 | **Newton–Raphson** | No | Needs $f'$; quadratic convergence; multivariate OK |
| 4 | **Fixed-point (Picard)** | No | Needs $x=g(x)$ form; $|g'|<1$; multivariate OK |
| 5 | **Laguerre** | No | For **polynomial** roots + deflation |

> First four are iterative → user precision $\epsilon$, typically $10^{-4}$.

### 5.2 Bisection method

Simplest, relatively slower, but **guaranteed to converge** if bracketing is careful: $f(a)$, $f(b)$ opposite signs, $f$ continuous.

Things can go wrong: (i) $f(x)=0$ at an **extremum** (both $f(a),f(b)$ same sign), or (ii) **multiple roots** in interval.

**Achieving bracketing:**

1. Choose $[a,b]$, compute $f(a),f(b)$. If $f(a)\cdot f(b)<0$, bracketing done.
2. If $f(a)\cdot f(b)>0$ (same sign), check $|f(a)|\lessgtr|f(b)|$:
   - If $|f(a)|<|f(b)|$, shift $a$ left: $a' = a-\beta(b-a)$.
   - If $|f(a)|>|f(b)|$, shift $b$ right: $b' = b+\beta(b-a)$.

**Bisection steps:**

1. Choose bracket $[a,b]$, $a<b$, $f(a)\cdot f(b)<0$. If $|b-a|<\epsilon$, done! Also check $f(a)$ and/or $f(b)<\delta$.
2. Else bisect $c=(a+b)/2$. New interval:
   - $[a,c]\equiv[a_1,b_1]$ if $f(c)\cdot f(a)<0$,
   - $[c,b]\equiv[a_1,b_1]$ if $f(c)\cdot f(b)<0$.
   
   Check $|b_1-a_1|<\epsilon$ along with $f(a_1\text{ or }b_1)<\delta$.
3. Keep bisecting/re-adjusting till $|b_n-a_n|<\epsilon$ **and** $f(a_n\text{ or }b_n)<\delta$.
4. After $n$ bisections:

$$\boxed{|b_n-a_n| = \frac{|b-a|}{2^n} \le \epsilon}$$

### 5.3 Regula falsi (false position) method

Interpolation to converge **faster than bisection**: slope of line joining $[a_0,b_0]$ that brackets the root; point $c_0$ where line crosses $x$-axis becomes new $a$ or $b$ depending on sign of $f(c_0)$.

- Unlike bisection, new boundary $a_n$ or $b_n$ is **directly given by $f(x)$**.
- Always converges, improved speed. As $(a,b)$ get close, can lose significant digits.

**Steps:**

1. Choose bracket $[a_0,b_0]$, $a_0<b_0$, $f(a)\cdot f(b)<0$. If $|b_0-a_0|<\epsilon$, done! Also check $f(a_0)$ and/or $f(b_0)<\delta$.
2. Slope of line joining $f(a_0)$, $f(b_0)$; $c_0$ where $y(c_0)=0$:

$$m = \frac{f(b_0)-f(a_0)}{b_0-a_0} = \frac{f(b_0)-f(c_0)}{b_0-c_0} \;\Longrightarrow\; \boxed{c_0 = b_0-(b_0-a_0)\,\frac{f(b_0)}{f(b_0)-f(a_0)}}$$

   Reference point can as well be $f(a_0)$. If $f(x)$ is convex/concave, one of $a,b$ stays fixed, the other varies. After $n$-th step:

$$\boxed{c_n = b_n-(b_n-a_n)\,\frac{f(b_n)}{f(b_n)-f(a_n)}}$$

3. If $f(a_n)\cdot f(c_n)<0$: root left of $c_n$ → $b_{n+1}=c_n$, $a_{n+1}=a_n$. If $|c_{n-1}-c_n|<\epsilon$ then $c_n$ is the root, $f(c_n)\approx 0$.
4. If $f(b_n)\cdot f(c_n)<0$: root right of $c_n$ → $a_{n+1}=c_n$, $b_{n+1}=b_n$. If $|c_{n-1}-c_n|<\epsilon$ then $c_n$ is the root.

### 5.4 Fixed-point / Picard iteration

Finding root $f(x_0)=0$.

- **Fixed point** of $g(x)$: $x=g(x)$.
- Replace $f(x)=0$ by $x=g(x)$. Instead of where $f$ crosses $x$-axis, find where $y=g(x)$ intersects $y=x$. Intersection $x_0$ = fixed point = root of $f$.
- For convergence (stable/attractive fixed point):

$$\boxed{|g'(x_0)|<1}$$

else unstable/repulsive.

Iteration from initial guess $x_1,x_2,\dots$:

$$\boxed{x_{n+1} = g(x_n)} \;\Longrightarrow\; x_{n+1}\approx g(x_n)\ \text{i.e.}\ |x_{n+1}-x_n|<\epsilon$$

For a given $f(x)=0$ there can be **many equivalent fixed-point forms** with different $g(x)$.

Example:

$$f(x)=x^2-x-2=0 \;\Longrightarrow\; x=g(x)=x^2-2,\ \sqrt{x+2},\ 1+\frac{2}{x},\ \frac{x^2+2}{2x-1},\ \dots$$

> Not all choices converge equally; some introduce unwanted singularities!

Famous fixed-point problem (from slides): devotee starts with some flowers, magic pond **doubles** them each temple visit, offers equal number $X$ each time; at last temple after offering $X$ she has none left. Starting number? How many temples?

$$f(x) = 2x - X \qquad\text{(solve with various }X\text{)}$$

### 5.5 Newton–Raphson method

- Uses **both $f(x)$ and $f'(x)$**, **no bracketing** needed. Based on **Taylor expansion**.
- Works for **multivariate** functions. **Quadratic convergence** → near root, significant digits roughly **double each step**.
- Taylor-expand around initial guess $x_0$:

$$f(x) = f(x_0)+(x-x_0)f'(x_0)+\frac{1}{2!}(x-x_0)^2f''(x_0)+\cdots$$

If close to root, $(x-x_0)^2\approx 0$, stop at $f'$ term:

$$f(x_0)+(x-x_0)f'(x_0)=0 \;\Longrightarrow\; \boxed{x_1 \equiv x = x_0-\frac{f(x_0)}{f'(x_0)}}$$

$x_1$ is better than $x_0$ but needs a derivative — useful if derivative is cheap to evaluate/code. Iterate:

$$\boxed{x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)},\quad n=0,1,\dots}$$

**Secant variant** — finite-difference derivative:

$$f'(x) \approx \frac{1}{2h}\big[f(x+h)-f(x-h)\big]+O(h^3)$$

Needs **two** initial guesses $x_0,x_1$ (i.e. $x\pm h$).

**Newton–Raphson steps:**

1. Good guess $x_0$.
2. Evaluate $f(x)$, $f'(x)$ at $x_0$.
3. Update $x_{n+1}=x_n-f(x_n)/f'(x_n)$.
4. Repeat until $|x_{n+1}-x_n|<\epsilon$ and/or $f(x_n)<\delta$.

### 5.6 Multivariable systems (Fixed-point + Newton–Raphson)

Both fixed-point and Newton–Raphson extend to nonlinear systems:

$${\bf f}({\bf x}) = \begin{pmatrix}f_1(x_1,\dots,x_n)\\f_2(x_1,\dots,x_n)\\\vdots\\f_n(x_1,\dots,x_n)\end{pmatrix}=0,\qquad {\bf x}^T=[x_1,\dots,x_n]$$

> **No bracketing in higher dimensions** → guaranteed methods like bisection **don't exist**!

**Multivariable fixed-point:**

$$\left.\begin{aligned}
x_1&=g_1(x_1,\dots,x_n)\\x_2&=g_2(x_1,\dots,x_n)\\\vdots\\x_n&=g_n(x_1,\dots,x_n)
\end{aligned}\right\} \Longrightarrow {\bf x}={\bf g}({\bf x}) \Longrightarrow \boxed{{\bf x}^{k+1}={\bf g}({\bf x}^k)}$$

Stopping criterion:

$$\boxed{\frac{||{\bf x}^{k+1}-{\bf x}^k||}{||{\bf x}^{k+1}||}<\epsilon},\qquad ||{\bf x}^k||^2=(x_1)_k^2+\cdots+(x_n)_k^2$$

**Fixed-point example from slides:**

$$\begin{aligned}
x_1^2+x_1x_2 &= 10\\
x_2+3x_1x_2^2 &= 57
\end{aligned}\qquad\Rightarrow\qquad {\bf x}=\begin{pmatrix}x_1\\x_2\end{pmatrix}={\bf g}({\bf x})=\begin{pmatrix}g_1\\g_2\end{pmatrix}=\begin{pmatrix}(10-x_1^2)/x_2\\57-3x_1x_2^2\end{pmatrix}$$

From $(1.5,3.5)$ this choice **diverges**:

$$x_1=\tfrac{10-(1.5)^2}{3.5}=2.21429,\quad x_2=57-3(2.21429)(3.5)^2=-24.37516$$
$$x_1=\tfrac{10-(2.21429)^2}{-24.37516}=-0.20910,\quad x_2=57-3(-0.20910)(-24.37516)^2=429.709,\ \dots$$

A different ${\bf g}$ works well:

$$\begin{pmatrix}x_1\\x_2\end{pmatrix}=\begin{pmatrix}\sqrt{10-x_1x_2}\\\sqrt{(57-x_2)/3x_1}\end{pmatrix}$$

| $x_1$ | $x_2$ | $x'_1$ | $x'_2$ | $\epsilon$ |
|---|---|---|---|---|
| 1.5000 | 3.5000 | 2.1794 | 2.8605 | 0.259459 |
| 2.1794 | 2.8605 | 1.9405 | 3.0496 | 0.084286 |
| 1.9405 | 3.0496 | 2.0205 | 2.9834 | 0.028792 |
| … | … | … | … | … |
| 1.9999 | 3.0001 | 2.0000 | 3.0000 | 0.000046 |

**Multivariable Newton–Raphson + Jacobian:**

Taylor-expand $f_1$ about $i$-th iterate:

$$f_1({\bf x}^{(i+1)})\approx f_1({\bf x}^{(i)})+\left.\frac{\partial f_1}{\partial x_1}\right|_{{\bf x}^{(i)}}\!(x^{(i+1)}_1-x^{(i)}_1)+\cdots+\left.\frac{\partial f_1}{\partial x_n}\right|_{{\bf x}^{(i)}}\!(x^{(i+1)}_n-x^{(i)}_n)$$

Similarly for $f_2,\dots,f_n$:

$$\boxed{{\bf f}({\bf x}^{(i+1)})={\bf f}({\bf x}^{(i)})+J({\bf x}^{(i)})({\bf x}^{(i+1)}-{\bf x}^{(i)})=0}$$

$$\boxed{{\bf x}^{(i+1)}={\bf x}^{(i)}-J^{-1}({\bf x}^{(i)})\,{\bf f}({\bf x}^{(i)})}$$

**Jacobian:**

$$\boxed{[J({\bf x}^{(i)})]_{pq}=\left.\frac{\partial f_p}{\partial x_q}\right|_{{\bf x}^{(i)}}=\begin{pmatrix}\partial f_1/\partial x_1&\cdots&\partial f_1/\partial x_n\\\vdots&\ddots&\vdots\\\partial f_n/\partial x_1&\cdots&\partial f_n/\partial x_n\end{pmatrix}}$$

Compare: $x^{(i+1)}=x^{(i)}-(f'(x^{(i)}))^{-1}f(x^{(i)})$ vs ${\bf x}^{(i+1)}={\bf x}^{(i)}-J^{-1}({\bf x}^{(i)}){\bf f}({\bf x}^{(i)})$.

**Newton multivariable worked example (same system):**

$$\begin{aligned}
x_1^2+x_1x_2-10&=0\\x_2+3x_1x_2^2-57&=0
\end{aligned}\qquad\Rightarrow\qquad J=\begin{pmatrix}\partial f_1/\partial x_1&\partial f_1/\partial x_2\\\partial f_2/\partial x_1&\partial f_2/\partial x_2\end{pmatrix}=\begin{pmatrix}2x_1+x_2&x_1\\3x_2^2&1+6x_1x_2\end{pmatrix}$$

Initial guess $x^{(0)}_1=1,\ x^{(0)}_2=1$:

$${\bf f}({\bf x}^{(0)})=\begin{pmatrix}-8\\-53\end{pmatrix}\Rightarrow J=\begin{pmatrix}3&1\\3&7\end{pmatrix}\to J^{-1}=\frac{1}{18}\begin{pmatrix}7&-1\\-3&3\end{pmatrix}$$

$$\begin{pmatrix}x^{(1)}_1\\x^{(1)}_2\end{pmatrix}=\begin{pmatrix}1\\1\end{pmatrix}-\frac{1}{18}\begin{pmatrix}7&-1\\-3&3\end{pmatrix}\!\begin{pmatrix}-8\\-53\end{pmatrix}=\begin{pmatrix}1.1667\\8.5000\end{pmatrix}$$

| $n$ | $x_1$ | $x_2$ |
|---|---|---|
| 1 | 1.1667 | 8.5000 |
| 2 | 1.5670 | 3.6878 |
| 3 | 2.0108 | 2.8824 |
| 4 | 1.9992 | 3.0023 |
| 5 | 2.0000 | 3.0000 |

Converges quadratically, but cost of Jacobian is $N^2$ function evaluations + $N^3$ for solving linear equations.

### 5.7 Laguerre's method (polynomial roots) + deflation

Degree-$n$ polynomial has exactly $n$ roots $\alpha_i$:

$$P(x)=a_0+a_1x+\cdots+a_nx^n=(x-\alpha_1)\cdots(x-\alpha_n)$$

- Real $a_i$ → roots real or complex-conjugate pairs. Course restricts to **real roots only**.
- **Two-step:** (i) Laguerre to find $\alpha_1$, (ii) **deflate** to degree $n-1$:

$$P(x)=(x-\alpha_1)Q(x),\quad Q(x)=(x-\alpha_2)R(x)\ \Rightarrow\ P(x)=(x-\alpha_1)(x-\alpha_2)R(x),\ \text{etc.}$$

Roots of $Q(x)$ = remaining roots of $P(x)$. Each of $n$ steps: Laguerre finds $\alpha_i$, deflation finds remainder polynomial. Deflation uses **synthetic division**.

**Laguerre algorithm:**

1. Initial guess $\beta_0$.
2. If $P(\beta_0)\approx 0$, go to deflation.
3. Else compute:

$$\boxed{G=\frac{P'(\beta_k)}{P(\beta_k)},\quad H=G^2-\frac{P''(\beta_k)}{P(\beta_k)}\quad\Rightarrow\quad a=\frac{n}{G\pm\sqrt{(n-1)(nH-G^2)}}}$$

   Choose sign in denominator to give **larger absolute value**.
4. Set $\beta_{k+1}=\beta_k-a$.
5. Iterate till $|\beta_{k+1}-\beta_k|<\epsilon$, set $\alpha_1=\beta_k$. Check $P(\alpha_1)\approx 0$.
6. Deflate and repeat for $\alpha_2$, etc.:

$$\frac{P(x)}{x-\alpha_1}=Q(x),\quad \frac{Q(x)}{x-\alpha_2}=R(x),\ \text{etc.}$$

**Deflation by synthetic division — worked example:**

$$\frac{P(x)}{x-\alpha_1}=\frac{-x^3+3x^2-4}{x-2}\Rightarrow \text{divisor}=2,\ \text{coeffs}=[-1,3,0,-4]$$

- Arrange descending powers, **0 for missing power**.
- Write divisor + coefficients; bring down leading coefficient $-1$ below line.
- Multiply $-1\times 2=-2$, add to next coefficient $3+(-2)=1$, bring down; repeat till end.

```
-1   3   0   -4
2  + -2   2    4
----------------
-1   1   2    0   <- remainder must be 0 if alpha1=2 is a root
```

$$\frac{-x^3+3x^2-4}{x-2}=-x^2+x+2=Q(x)$$

Repeat with $Q(x)$ till final monomial $(x-\alpha_n)$:

```
-1   1   2
2  + -2  -2
------------
-1  -1   0   => R(x) = -x-1
```

$$P(x)=-x^3+3x^2-4=(x-2)(x-2)(-x-1)$$

**Find all roots (answers in slides):**

- $P(x)=6x^3-11x^2-26x+15$ → $x=3,\,-5/3,\,1/2$
- $P(x)=x^4-x^3-7x^2+x+6$ → $x=1,\,-2,\,3,\,-1$

---

## 6. Numerical Integration

### 6.1 Why integrate numerically?

Definite integrals overwhelm indefinite ones in physics/engineering — fields of charge/mass distributions, shuttle center-of-gravity, aircraft lift/drag, … — and can rarely be done by hand/analytically.

Numerically approximate to desired precision:

$$\boxed{I=\int_a^b f(x)\,dx}$$

$f(x)$ = (piecewise) continuous, well-behaved (smooth) over $[a,b]$. Definite integral = **area under the curve**.

**Reasons:**

1. $f(x)$ known only at certain points.
2. $f(x)$ known but analytic integration too difficult/impossible.
3. Analytic possible but numeric far easier to given accuracy.

> Numerical integration = **quadrature** (esp. 1-D). Evaluate integrand at finite points (equally spaced or not) in $[a,b]$, weighted sum:

$$\boxed{I=\int_a^b f(x)\,dx \approx \sum_{n=1}^{N}w(x_n)\,f(x_n)\equiv I_N}$$

$w(x_n)$ = weight function, $N$ = integration points. Choice of $N$ from desired $|I-I_N|<\epsilon$.

### 6.2 Origin of weights: Lagrange interpolation + strategy

Lagrange interpolating polynomial:

$$p(x)=f(x_0)\ell_0(x)+\cdots+f(x_N)\ell_N(x)$$

Integrating:

$$I_N=\int_a^b p(x)\,dx=\int_a^b\sum_{n=0}^{N}f(x_n)\ell_n(x)\,dx=\sum_{n=0}^{N}\left(f(x_n)\int_a^b\ell_n(x)\,dx\right)\equiv\sum_{n=0}^{N}f(x_n)\,w(x_n)$$

Techniques: **Midpoint, Trapezoidal, Simpson, Monte Carlo, Gaussian Quadrature**.

Strategy: approximate $f(x)$ by $p(x)=a_0+a_1x+a_2x^2+\cdots$:

- Midpoint: only constant $a_0$
- Trapezoidal: up to linear $a_0+a_1x$
- Simpson: up to quadratic $a_0+a_1x+a_2x^2$

### 6.3 Midpoint method

1. Divide $[a,b]$ into $N$ equal parts of width:

$$\boxed{h=(b-a)/N}$$

2. Midpoints (integration points):

$$x_1=\tfrac{(a)+(a+h)}{2},\ x_2=\tfrac{(a+h)+(a+2h)}{2},\ x_3=\tfrac{(a+2h)+(a+3h)}{2},\ \dots$$

   $f$ evaluated only at $x_n$; rectangle $h\times f(x_n)$ ≈ area over sub-interval.
3. Sum:

$$\boxed{M_N=\sum_{n=1}^{N}h\cdot f(x_n)\quad\Rightarrow\quad \lim_{N\to\infty}M_N=\int_a^b f(x)\,dx}$$

where $w(x_n)=w=1$ or $h$, constant for all $x_n$.

### 6.4 Trapezoidal method

Straight lines joining consecutive $f(x_n)$ → trapezoids.

1. $N$ intervals width $h$; endpoints $x_0,\dots,x_N$: $x_0=a,\ x_1=x_0+h,\ \dots,\ x_n=x_0+nh,\ \dots,\ x_N=b$.
2. Area of each trapezoid:

$$\boxed{T_n=\frac{h}{2}\big[f(x_{n-1})+f(x_n)\big]}$$

3. Sum:

$$T_N=\sum_{n=1}^{N}T_n=\frac{h}{2}\Big\{[f(x_0)+f(x_1)]+[f(x_1)+f(x_2)]+\cdots+[f(x_{N-1})+f(x_N)]\Big\}$$

Hence:

$$\boxed{T_N=\frac{h}{2}\big[f(x_0)+2f(x_1)+\cdots+2f(x_{N-1})+f(x_N)\big]=\sum_{n}w(x_n)f(x_n)}$$

$$\lim_{N\to\infty}T_N=\int_a^b f(x)\,dx$$

Weights: $w(x_0)=w(x_N)=1$ or $h/2$; $w(x_1)=\cdots=w(x_{N-1})=2$ or $h$.

> Contrary to expectation, trapezoidal tends to be **less accurate than midpoint**, particularly when curve is strictly concave/convex over limits.

**Exercise:** for $N=4,8,12,20$ use both Midpoint and Trapezoidal, compare with analytic:

$$\int_1^2\frac{1}{x}\,dx=0.69314718,\qquad \int_0^{\pi/2}x\cos x\,dx=\frac{\pi}{2}-1$$

### 6.5 Simpson's method (1/3 rule)

- Midpoint: rectangle (piecewise **constant**).
- Trapezoidal: trapezium (piecewise **linear**).
- **Simpson: piecewise quadratic** → need 3 points: two (sub)interval boundaries + their midpoint:

$$\int_{x_0}^{x_2}f(x)\,dx \to (x_0,f(x_0)),\,(x_2,f(x_2)),\,(x_1,f(x_1)),\quad x_1=\frac{x_0+x_2}{2}$$

With $h=(x_2-x_0)/2$, integrating $a_2x^2+a_1x+a_0$:

$$\begin{aligned}
\int_{x_0}^{x_2}f(x)\,dx &\approx \int_{x_0}^{x_2}(a_2x^2+a_1x+a_0)\,dx = \left[\tfrac{a_2}{3}x^3+\tfrac{a_1}{2}x^2+a_0x\right]_{x_0}^{x_2}\\
&=\tfrac{a_2}{3}(x_2^3-x_0^3)+\tfrac{a_1}{2}(x_2^2-x_0^2)+a_0(x_2-x_0)\\
&=\tfrac{x_2-x_0}{6}\big[2a_2(x_2^2+x_2x_0+x_0^2)+3a_1(x_2+x_0)+6a_0\big]\\
&=\tfrac{h}{3}\big[(a_2x_2^2+a_1x_2+a_0)+(a_2x_0^2+a_1x_0+a_0)\\
&\qquad\qquad +a_2(x_2^2+2x_2x_0+x_0^2)+2a_1(x_2+x_0)+4a_0\big]\\
&=\tfrac{h}{3}\big[f(x_2)+f(x_0)+a_2(2x_1)^2+2a_1(2x_1)+4a_0\big]\\
&=\boxed{\tfrac{h}{3}\big[f(x_0)+4f(x_1)+f(x_2)\big]}
\end{aligned}$$

Similarly $\int_{x_2}^{x_4}f(x)\,dx=\frac{h}{3}[f(x_2)+4f(x_3)+f(x_4)]$.

Summing $[x_0,x_2],[x_2,x_4],\dots,[x_{N-2},x_N]$ → **Simpson's 1/3 rule**:

$$\boxed{S_N=\frac{h}{3}\big[f(x_0)+4f(x_1)+2f(x_2)+4f(x_3)+\cdots+2f(x_{N-2})+4f(x_{N-1})+f(x_N)\big]}$$

$$\boxed{S_N=\frac{h}{3}\left[f(x_0)+4\!\!\sum_{n\ \mathrm{odd}}^{N-1}f(x_n)+2\!\!\sum_{n\ \mathrm{even}}^{N-2}f(x_n)+f(x_N)\right]=\sum_n w(x_n)f(x_n)}$$

$$\lim_{N\to\infty}S_N=\int_a^b f(x)\,dx$$

Weights: $w(x_0)=w(x_N)=1$, $w(x_n)=4$ for odd $n$, $2$ for even $n$.

> N.B. Simpson's 1/3 needs **even number of intervals, odd number of grid points**. There is also Simpson's **3/8 rule** using cubic $a_3x^3$.

Interesting relation:

$$\boxed{S_{2N}=\tfrac{2}{3}M_N+\tfrac{1}{3}T_N}$$

### 6.6 Error in numerical integration

Can't define absolute $|I_N-I|$ / relative $|I_N-I|/|I|\times 100\%$ since $I$ unknown (that's why we integrate numerically). Instead define **upper bound** of error. $f$ continuous over $[a,b]$:

- finite $f''(x)$ for Midpoint/Trapezoidal,
- finite $f''''(x)$ for Simpson.

$$\boxed{\text{Midpoint: Error in }M_N \le \frac{(b-a)^3}{24N^2}\,|f''(x)|_{\max}}$$

$$\boxed{\text{Trapezoidal: Error in }T_N \le \frac{(b-a)^3}{12N^2}\,|f''(x)|_{\max}}$$

$$\boxed{\text{Simpson: Error in }S_N \le \frac{(b-a)^5}{180N^4}\,|f''''(x)|_{\max}}$$

**What if $f''\cdots(x)=0$?** Doesn't mean error is zero — just means error-bound estimate can't be made.

**Worked example:** $\int_0^1 x^2\,dx=0.333$. Ignore bound for a moment. $N=4$: $[0,\tfrac14],[\tfrac14,\tfrac24],[\tfrac24,\tfrac34],[\tfrac34,1]$, $h=\tfrac14$. Midpoints $\tfrac18,\tfrac38,\tfrac58,\tfrac78$:

$$M_4=\tfrac14\left[f(\tfrac18)+f(\tfrac38)+f(\tfrac58)+f(\tfrac78)\right]=\tfrac14\left[(\tfrac18)^2+(\tfrac38)^2+(\tfrac58)^2+(\tfrac78)^2\right]=\tfrac{21}{64}=0.328$$

Absolute error $|0.333-0.328|\approx 0.0052$, relative $1.56\%$.

For upper bound $0.001$: $f(x)=x^2\Rightarrow f''=2\Rightarrow|f''|_{\max}=2$:

$$0.001=\frac{(1-0)^3}{24N^2}\times 2 \;\Rightarrow\; N=9$$

**Choosing $N$ from the bound — two things:**

1. $N$ = smallest integer $\ge$ computed $N$ (i.e. `ceil`).
2. Actual estimate may be **much better** than the bound indicates.

> Strange but often true: taking $N$ **larger** than from error bounds can actually **deteriorate** accuracy, contrary to $N\to\infty$ limit.

**Comparison example:** $f(x)=1/x$ on $[1,2]$: $f''(x)=2/x^3$, $f''''(x)=24/x^5$ → $|f''|_{\max}=2$, $|f''''|_{\max}=24$. For max error $0.0001$:

- Midpoint: $0.0001=\frac{(2-1)^3}{24N^2}\times 2 \to N=29$
- Trapezoidal: $0.0001=\frac{(2-1)^3}{12N^2}\times 2 \to N=41$
- Simpson: $0.0001=\frac{(2-1)^5}{180N^4}\times 24 \to N=7$

### 6.7 Monte Carlo integration

A class using **repeated random sampling**. Widely used (perhaps the **only**) method for **higher-dimensional / multivariable** integrals.

To estimate $\int_a^b f(x)\,dx$, define an **estimator**. For random $X\sim p(X)$:

$$\boxed{F_N \equiv \frac{1}{N}\sum_{i=1}^{N}\frac{f(X_i)}{p(X_i)}\quad\Rightarrow\quad \langle F_N\rangle=\int_a^b F_N\,p(X)\,dX\equiv\int_a^b f(x)\,dx}$$

For **uniform** PDF on $[a,b]$: $\int_a^b p(X)\,dX=\int_a^b C\,dX=1\Rightarrow p(X)=C=1/(b-a)$:

$$\boxed{F_N=\frac{b-a}{N}\sum_{i=1}^{N}f(X_i)}$$

Looks very similar to Midpoint expression.

**Steps:**

1. Choose $N$ (10, 20, 50, whatever).
2. Draw $N$ random $X_i$ in $[a,b]$. If pRNG gives $\xi\in[0,1]$, convert:

$$\boxed{X=a+(b-a)\,\xi}$$

3. For each $X_i$ compute $f(X_i)$, then $F_N$ and $\sigma_f$:

$$\boxed{\sigma_f^2=\frac{1}{N}\sum_{i=1}^{N}f(X_i)^2-\left(\frac{1}{N}\sum_{i=1}^{N}f(X_i)\right)^2}$$

4. Tabulate/plot $F_N$ vs $N$, track $\sigma_f$. As $N$ grows, $F_N$ converges but $\sigma_f$ falls slowly $\sim 1/\sqrt{N}$.
5. Increase $N$ by $10\times$ or whatever, repeat.

### 6.8 Gaussian Quadrature

Midpoint/Trapezoidal/Simpson = collectively **closed Newton–Cotes**: integrate Lagrange polynomial interpolating $f$ at $N$ **equally spaced** points:

$$I=\int_a^b f(x)\,dx\approx\sum_{n=1}^{N}w(x_n)f(x_n),\quad w(x_n)=\int_a^b L_n(x)\,dx$$

$$L_n(x)=\prod_{k=1,k\ne n}^{N}\frac{x-x_n}{x_k-x_n},\qquad x_n=\text{nodes},\ w(x_n)=\text{weights}$$

If points are **not equally spaced** (only must lie in $[a,b]$), we have $2n$ parameters: $x_n$ **and** $w(x_n)$.

**Gaussian quadrature rule:** express integrand in **orthogonal basis polynomials** over the interval. E.g.:

$$\int_{-1}^{1}f(x)\,dx\approx\sum_{n=1}^{N}w(x_n)f(x_n)$$

$x_n$ = **roots of an orthogonal polynomial** for that interval. Optimal — fits all polynomials up to degree **$2n-1$** exactly (defined by $2n$ parameters).

On $[-1,1]$, $f$ well-approximated by **Legendre** $P_n(x)$; $x_n$ = roots, weights:

| $n$ | $x_n$ (roots of $P_n$) | weights $w(x_n)$ |
|---|---|---|
| 1 | $0$ | $2$ |
| 2 | $\pm\sqrt{1/3}$ | $1$ |
| 3 | $0,\ \pm\sqrt{3/5}$ | $8/9,\ 5/9$ |
| 4 | $\pm\sqrt{\tfrac37-\tfrac27\sqrt{\tfrac65}},\ \pm\sqrt{\tfrac37+\tfrac27\sqrt{\tfrac65}}$ | $(18+\sqrt{30})/36,\ (18-\sqrt{30})/36$ |
| 5 | $\cdots$ | $\cdots$ |

For general $\int_a^b$ (map $[a,b]\to[-1,1]$):

$$\int_a^b f(t)\,dt=\frac{b-a}{2}\int_{-1}^{1}f\!\left(\frac{b-a}{2}x+\frac{b+a}{2}\right)dx$$

$$\boxed{\int_a^b f(t)\,dt\approx\frac{b-a}{2}\sum_{n=1}^{N}w(x_n)\,f\!\left(\frac{b-a}{2}x_n+\frac{b+a}{2}\right)}$$

For $\int_0^\infty$ or $\int_{-\infty}^{\infty}$ instead of $\int_{-1}^{1}$, use **Laguerre** $L_n(x)$ and **Hermite** $H_n(x)$ polynomials.

**Steps:** instead of midpoints/boundaries:

1. Evaluate $f$ at roots $x_0,\dots,x_N$ of the orthogonal polynomial.
2. Weighted sum $\sum_n w(x_n)f(x_n)\to\int_a^b f\,dx$.

**Efficiency demo:** $\int_{-1}^{1}x\,e^x\,dx=2/e=0.73575888$:

| Simpson $N$ | $I_S$ | Gauss $N$ | $I_G$ |
|---|---|---|---|
| 20 | 0.73576450 | 4 | 0.73575650 |
| 40 | 0.73575923 | 5 | 0.73575887 |
| 60 | 0.73575895 | 6 | 0.73575888 |
| 80 | 0.73575890 | | |
| 100 | 0.73575888 | | |
| 120 | 0.73575888 | | |

> Gaussian quadrature hits full accuracy with $N=6$ where Simpson needs $N\approx 100$.

---


---

## 7. Ordinary Differential Equations (ODEs)

> Source: `diffeq.pdf` (Part I) & Course Assignments 7A–7C.

### 7.1 Initial Value Problems (IVPs) & Problem Formulation
A first-order ordinary differential equation with initial value is given by:
$$\frac{dy}{dx} = f(x, y), \quad y(x_0) = y_0$$
The derivative $dy/dx$ defines the slope of the solution curve $y(x)$ at $(x_0, y_0)$.

Higher-order ODEs are converted into coupled systems of first-order ODEs. For example, a 2nd-order equation:
$$\frac{d^2x}{dt^2} = g(t, x, v) \implies \begin{cases} \frac{dx}{dt} = v \\ \frac{dv}{dt} = g(t, x, v) \end{cases}$$
State vector $\mathbf{y}(t) = [x(t), v(t)]^T$ is evolved simultaneously.

---

### 7.2 Forward (Explicit) Euler Method
Based on the forward Taylor expansion about $x_n$:
$$y(x_n + h) = y(x_n) + h y'(x_n) + \frac{h^2}{2} y''(x_n) + \dots = y_n + h f(x_n, y_n) + O(h^2)$$

- **Algorithm:** $y_{n+1} = y_n + h f(x_n, y_n)$
- **Order:** Local truncation error $O(h^2)$; global accumulated error $O(h)$ (first-order).
- **Stability:** Conditionally stable. For test equation $y' = -\lambda y$, numerical stability requires $|1 - h\lambda| \le 1 \implies h \le \frac{2}{\lambda}$. If $h$ exceeds this limit, numerical oscillations explode.

---

### 7.3 Backward (Implicit) Euler Method
Re-defines $dy/dx$ using the backward difference quotient:
$$y(x_n + h) = y(x_n) + h f(x_{n+1}, y(x_n + h)) \implies y_{n+1} = y_n + h f(x_{n+1}, y_{n+1})$$

Because $y_{n+1}$ appears on both sides:
1. **Linear ODEs:** Rearrange analytically for $y_{n+1}$. E.g., for $y' = -\lambda y$:
   $$y_{n+1} = \frac{y_n}{1 + h\lambda}$$
   Unconditionally stable for all $h > 0$ since $|1 + h\lambda| > 1$.
2. **Non-Linear ODEs:** Solved at each time step using 1D Newton-Raphson:
   $$g(y) = y - y_n - h f(x_{n+1}, y) = 0, \quad g'(y) = 1 - h \frac{\partial f}{\partial y}(x_{n+1}, y)$$

---

### 7.4 Heun's Predictor-Corrector Method
Combines forward Euler prediction with trapezoidal slope averaging:
1. **Predictor (Forward Euler):**
   $$y_{n+1}^p = y_n + h f(x_n, y_n)$$
2. **Corrector (Trapezoidal Average):**
   $$y_{n+1} = y_n + \frac{h}{2} \left[ f(x_n, y_n) + f(x_{n+1}, y_{n+1}^p) \right]$$
- **Global Error:** $O(h^2)$ (second-order accuracy). Significantly more accurate than simple Forward Euler with minimal overhead (two function evaluations per step).

---

### 7.5 Runge-Kutta 4th-Order (RK4) Method
The gold standard for initial-value problems. Evaluates four trial slopes across the step $[x_n, x_n + h]$ and combines them with Simpson's 1/3-rule weights $(1, 2, 2, 1)$:

$$\begin{aligned}
k_1 &= f(x_n, y_n) && \text{(slope at beginning)} \\
k_2 &= f\left(x_n + \frac{h}{2}, y_n + \frac{h}{2} k_1\right) && \text{(midpoint slope via } k_1 \text{)} \\
k_3 &= f\left(x_n + \frac{h}{2}, y_n + \frac{h}{2} k_2\right) && \text{(midpoint slope via } k_2 \text{)} \\
k_4 &= f(x_n + h, y_n + h k_3) && \text{(endpoint slope via } k_3 \text{)}
\end{aligned}$$

**Weighted Update:**
$$y_{n+1} = y_n + \frac{h}{6} (k_1 + 2k_2 + 2k_3 + k_4)$$
- **Accuracy:** Fourth order $O(h^4)$ global error, $O(h^5)$ local truncation error.
- **Vector Formulation:** Evaluates vectors $\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3, \mathbf{k}_4 \in \mathbb{R}^m$ for coupled dynamical systems.

---

### 7.6 Coupled ODE Applications
1. **Damped Simple Harmonic Oscillator:**
   $$\ddot{x} + \mu \dot{x} + \omega^2 x = 0, \quad x(0) = x_0, \; v(0) = v_0$$
   Total mechanical energy $E(t) = \frac{1}{2} m v^2 + \frac{1}{2} k x^2$. RK4 accurately captures both exponential amplitude decay and energy dissipation $\frac{dE}{dt} = -m \mu v^2 \le 0$.
2. **Chaotic Dynamics (Rössler Attractor):**
   $$\frac{dx}{dt} = -y - z, \quad \frac{dy}{dt} = x + ay, \quad \frac{dz}{dt} = b + z(x - c)$$
   For $a = 0.2, b = 0.2, c = 5.7$, RK4 traces a strange attractor with sensitive dependence on initial conditions.

---

### 7.7 Two-Point Boundary Value Problems (BVPs) & Shooting Method
Given $y'' = f(x, y, y')$ with Dirichlet boundary conditions $y(a) = \alpha, y(b) = \beta$:
1. Convert to IVP by guessing the initial slope $z = y'(a)$.
2. Pick two trial slopes $z_1$ and $z_2$, and integrate using RK4 to obtain boundary values $y(b; z_1)$ and $y(b; z_2)$.
3. Use linear interpolation (Lagrange formula) to find the exact initial slope $z^*$:
   $$z^* = z_1 + (z_2 - z_1) \frac{\beta - y(b; z_1)}{y(b; z_2) - y(b; z_1)}$$
4. One final RK4 integration with $y'(a) = z^*$ satisfies the target boundary condition $y(b) = \beta$ exactly (for linear ODEs; iterative shooting for non-linear ODEs).

---

### 7.8 1D Finite Element Method (FEM)
Solves boundary value problems (e.g. Poisson $y''(x) = f(x)$ on $[a, b]$, $y(a) = y_a, y(b) = y_b$) by minimizing the energy functional / weak variational form:
$$\int_a^b \left( -y'(x) v'(x) - f(x) v(x) \right) dx = 0$$
Using piecewise linear hat functions $\phi_i(x)$ on a mesh of $n$ elements:
$$K \mathbf{y}_{\text{int}} = \mathbf{F}$$
Stiffness matrix $K$ is symmetric positive-definite tridiagonal:
$$K_{ii} = \frac{2}{h}, \quad K_{i, i\pm 1} = -\frac{1}{h}$$
Solved efficiently and stably using Cholesky factorization ($K = L L^T$).

---

## 8. Partial Differential Equations (PDEs)

> Source: `diffeq.pdf` (Part II) & Assignment 7D.

### 8.1 PDE Classification
- **Parabolic:** $\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2}$ (Heat / diffusion equations; dissipation, smoothing).
- **Hyperbolic:** $\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$ (Wave equations; wave propagation, finite speed $c$).
- **Elliptic:** $\nabla^2 u = \rho(x, y)$ (Laplace / Poisson equations; steady-state equilibrium, boundary-driven).

---

### 8.2 1D Heat Equation: Explicit FTCS Scheme
Discretizing $u_t = u_{xx}$ with spatial step $h_x = L / N_x$ and time step $h_t = T / N_t$:
$$\frac{u_i^{j+1} - u_i^j}{h_t} = \frac{u_{i+1}^j - 2 u_i^j + u_{i-1}^j}{h_x^2}$$
$$u_i^{j+1} = \alpha (u_{i+1}^j + u_{i-1}^j) + (1 - 2\alpha) u_i^j, \quad \alpha = \frac{h_t}{h_x^2}$$

- **Von Neumann Stability Analysis:** Substituting Fourier mode $u_i^j = \xi^j e^{I k i h_x}$ yields amplification factor $\xi = 1 - 4\alpha \sin^2(k h_x / 2)$.
- **Stability Criterion:** $|\xi| \le 1 \iff \alpha \le 0.5$. If $\alpha > 0.5$, high-frequency modes explode with alternating sign oscillations.

---

### 8.3 1D Heat Equation: Implicit BTCS Scheme & Thomas Algorithm
Backward-Time Central-Space evaluates the second spatial derivative at the forward time level $j+1$:
$$\frac{u_i^{j+1} - u_i^j}{h_t} = \frac{u_{i+1}^{j+1} - 2 u_i^{j+1} + u_{i-1}^{j+1}}{h_x^2}$$
$$-\alpha u_{i-1}^{j+1} + (1 + 2\alpha) u_i^{j+1} - \alpha u_{i+1}^{j+1} = u_i^j$$

- **Stability:** Unconditionally stable for all $\alpha > 0$.
- **Solution Method:** Forms a strictly diagonally dominant tridiagonal system:
  $$\begin{bmatrix} 1+2\alpha & -\alpha & & \\ -\alpha & 1+2\alpha & -\alpha & \\ & \ddots & \ddots & \ddots \end{bmatrix} \begin{bmatrix} u_1^{j+1} \\ u_2^{j+1} \\ \vdots \end{bmatrix} = \mathbf{d}^j$$
  Solved in $O(N_x)$ time via the **Thomas Algorithm (TDMA)**:
  1. Forward sweep: $c'_0 = c_0 / b_0$, $c'_i = c_i / (b_i - a_i c'_{i-1})$; $d'_0 = d_0 / b_0$, $d'_i = (d_i - a_i d'_{i-1}) / (b_i - a_i c'_{i-1})$.
  2. Back substitution: $x_{n-1} = d'_{n-1}$, $x_i = d'_i - c'_i x_{i+1}$.

---

### 8.4 Crank-Nicolson Scheme
Averages explicit (level $j$) and implicit (level $j+1$) spatial differences:
$$\frac{u_i^{j+1} - u_i^j}{h_t} = \frac{1}{2} \left[ \frac{u_{i+1}^{j+1} - 2u_i^{j+1} + u_{i-1}^{j+1}}{h_x^2} + \frac{u_{i+1}^j - 2u_i^j + u_{i-1}^j}{h_x^2} \right]$$
$$-\frac{\alpha}{2} u_{i-1}^{j+1} + (1 + \alpha) u_i^{j+1} - \frac{\alpha}{2} u_{i+1}^{j+1} = \frac{\alpha}{2} u_{i-1}^j + (1 - \alpha) u_i^j + \frac{\alpha}{2} u_{i+1}^j$$

- **Accuracy:** Second order in time and space: $O(h_t^2 + h_x^2)$.
- **Stability:** Unconditionally stable for all $\alpha > 0$.
- **Matrix Evolution Form:** $A \mathbf{u}^{j+1} = B \mathbf{u}^j + \mathbf{b}_{\text{bc}} \implies \mathbf{u}^{j+1} = (A^{-1} B) \mathbf{u}^j + A^{-1} \mathbf{b}_{\text{bc}}$.

---

### 8.5 1D Hyperbolic Wave Equation
$$u_{tt} = c^2 u_{xx}, \quad u(x, 0) = u_0(x), \; u_t(x, 0) = v_0(x)$$
Central differences in both time and space:
$$\frac{u_i^{j+1} - 2 u_i^j + u_i^{j-1}}{h_t^2} = c^2 \frac{u_{i+1}^j - 2 u_i^j + u_{i-1}^j}{h_x^2}$$
$$u_i^{j+1} = 2(1 - r^2) u_i^j + r^2(u_{i+1}^j + u_{i-1}^j) - u_i^{j-1}, \quad r = \frac{c h_t}{h_x}$$

- **Courant-Friedrichs-Lewy (CFL) Condition:** $r = c h_t / h_x \le 1.0$. The numerical domain of dependence must contain the physical domain of dependence. If $r > 1$, energy grows without bound.
- **First Step ($j = 1$):** Taylor expansion using velocity: $u_i^1 = u_i^0 + h_t v_0(x_i) + \frac{1}{2} r^2 (u_{i+1}^0 - 2 u_i^0 + u_{i-1}^0)$.

---

### 8.6 2D Elliptic Equations (Laplace / Poisson)
$$\nabla^2 u = u_{xx} + u_{yy} = \rho(x, y)$$
5-point finite-difference stencil on square grid ($h_x = h_y = h$):
$$u_{i+1, j} + u_{i-1, j} + u_{i, j+1} + u_{i, j-1} - 4 u_{i, j} = h^2 \rho_{i, j}$$
Solved using Successive Over-Relaxation (SOR):
$$u_{i, j}^{(k+1)} = (1 - \omega) u_{i, j}^{(k)} + \frac{\omega}{4} \left( u_{i+1, j}^{(k)} + u_{i-1, j}^{(k+1)} + u_{i, j+1}^{(k)} + u_{i, j-1}^{(k+1)} - h^2 \rho_{i, j} \right)$$
Optimal $\omega \in (1, 2)$ reduces iteration counts by an order of magnitude.

---

## 9. Eigenvalues & Eigenvectors

> Source: `eigv.pdf` & Assignment 8.

### 9.1 The Algebraic Eigenvalue Problem
For an $n \times n$ matrix $A$:
$$A \mathbf{x} = \lambda \mathbf{x} \iff (A - \lambda I) \mathbf{x} = \mathbf{0}$$
Analytically solved via the characteristic polynomial $p_A(\lambda) = \det(A - \lambda I) = 0$. However:
> **Word of wisdom (from lecture):** In general, eigenvalues cannot be determined using a finite number of rational operations (Abel-Ruffini theorem for $n \ge 5$). Numerical methods must be iterative.

Small changes in general matrices can cause massive changes in eigenvalues. If $A$ is **symmetric**, eigenvalues are well-conditioned and robust against perturbations.

---

### 9.2 Power Iteration Method
An iterative algorithm designed to isolate the **dominant eigenvalue** (largest in absolute magnitude $|\lambda_1| > |\lambda_2| \ge \dots$) and its eigenvector $\mathbf{v}_1$.

#### Algorithm:
1. Start with initial guess $\mathbf{x}_0$ (not orthogonal to $\mathbf{v}_1$). Normalize: $\mathbf{x}_0 \leftarrow \mathbf{x}_0 / \|\mathbf{x}_0\|_2$.
2. Repeated multiplication: $\mathbf{x}_{k+1} = A \mathbf{x}_k$.
3. Expand in eigenvector basis:
   $$\mathbf{x}_k = c_1 \lambda_1^k \mathbf{v}_1 + c_2 \lambda_2^k \mathbf{v}_2 + \dots = \lambda_1^k \left[ c_1 \mathbf{v}_1 + c_2 \left( \frac{\lambda_2}{\lambda_1} \right)^k \mathbf{v}_2 + \dots \right]$$
   As $k \to \infty$, since $|\lambda_2 / \lambda_1| < 1$, all subdominant components vanish.
4. **Rayleigh Quotient:** Evaluates $\lambda_1$ using test vector $\mathbf{y}$ (typically $\mathbf{y} = \mathbf{x}_k$):
   $$\lambda_1^{(k)} = \frac{A \mathbf{x}_k \cdot \mathbf{y}}{\mathbf{x}_k \cdot \mathbf{y}} = \frac{\mathbf{x}_{k+1} \cdot \mathbf{y}}{\mathbf{x}_k \cdot \mathbf{y}}$$
5. Normalize $\mathbf{x}_{k+1} \leftarrow \mathbf{x}_{k+1} / \|\mathbf{x}_{k+1}\|_2$ at **every step** to avoid floating-point overflow.
6. Stop when $|\lambda_1^{(k)} - \lambda_1^{(k-1)}| < \epsilon$.

---

### 9.3 Matrix Deflation for Subdominant Eigenvalues
Once $\lambda_1$ and normalized eigenvector $\mathbf{u}_1 = \mathbf{v}_1 / \|\mathbf{v}_1\|_2$ are known, deflate the matrix by subtracting the dominant projection:
$$A' = A - \lambda_1 \mathbf{u}_1 \mathbf{u}_1^T$$
- **Property:** $A'$ has eigenvalues $0, \lambda_2, \lambda_3, \dots, \lambda_n$ with the same eigenvectors as $A$.
- Power iteration on $A'$ converges to the second dominant eigenvalue $\lambda_2$.
- *Caution:* Deflation beyond $\lambda_2$ accumulates roundoff error; QR iteration or Jacobi eigenvalue methods are preferred for full spectra.

---

### 9.4 Tridiagonal Toeplitz Spectrum (Analytical Benchmark)
For an $n \times n$ symmetric tridiagonal matrix with subdiagonal $a$, main diagonal $b$, and superdiagonal $c = a$:
$$\lambda_k = b + 2\sqrt{ac} \cos\left( \frac{k\pi}{n + 1} \right), \quad v_{k, i} = \sin\left( \frac{ik\pi}{n + 1} \right), \quad k = 1, \dots, n$$
For $a = c = -1, b = 2, n = 5$:
- $k=1$: $\lambda_1 = 2 + 2\cos(\pi/6) = 2 + \sqrt{3} \approx 3.732051$
- $k=2$: $\lambda_2 = 2 + 2\cos(2\pi/6) = 2 + 1 = 3.000000$
- $k=3$: $\lambda_3 = 2 + 2\cos(3\pi/6) = 2.000000$

---

## 10. Interpolation & Data Fitting

> Source: `datafit.pdf` & Assignments 9A–9B.

### 10.1 Lagrange Interpolation
Given $N+1$ points $(x_0, y_0), (x_1, y_1), \dots, (x_N, y_N)$ with distinct $x_i$, the unique polynomial $P_N(x)$ of degree $\le N$ passing through all points is:
$$P_N(x) = \sum_{i=0}^N y_i \ell_i(x), \quad \ell_i(x) = \prod_{k \ne i} \frac{x - x_k}{x_i - x_k}$$
- **Kronecker Property:** $\ell_i(x_j) = \delta_{ij} \implies P_N(x_i) = y_i$.
- Polynomial coefficients in ascending power order $[c_0, c_1, \dots, c_N]$ where $P_N(x) = \sum c_m x^m$ are extracted via polynomial convolution.

---

### 10.2 Principle of Weighted Least Squares
When fitting experimental data $(x_i, y_i)$ with measurement uncertainties $\sigma_i$, we minimize the chi-square sum:
$$\chi^2 = \sum_{i=1}^N \left( \frac{y_i - f(x_i; \mathbf{a})}{\sigma_i} \right)^2 = \sum_{i=1}^N w_i \left( y_i - f(x_i; \mathbf{a}) \right)^2, \quad w_i = \frac{1}{\sigma_i^2}$$

---

### 10.3 Linear Regression ($y = mx + c$)
Minimizing $\chi^2$ with respect to slope $m$ and intercept $c$:
$$\begin{aligned}
S &= \sum w_i, & S_x &= \sum w_i x_i, & S_y &= \sum w_i y_i \\
S_{xx} &= \sum w_i x_i^2, & S_{xy} &= \sum w_i x_i y_i, & S_{yy} &= \sum w_i y_i^2
\end{aligned}$$
$$\Delta = S S_{xx} - S_x^2, \quad m = \frac{S S_{xy} - S_x S_y}{\Delta}, \quad c = \frac{S_{xx} S_y - S_x S_{xy}}{\Delta}$$

**Pearson Correlation Coefficient ($r^2$):**
$$r^2 = \frac{(S S_{xy} - S_x S_y)^2}{\Delta (S S_{yy} - S_y^2)}$$
Linearized fitting for non-linear physical models:
- **Power law:** $y = a x^b \implies \ln y = \ln a + b \ln x$
- **Exponential:** $y = a e^{-bx} \implies \ln y = \ln a - b x$
Comparing $r^2$ identifies which physical model best captures the data.

---

### 10.4 Polynomial Least Squares & Normal Equations
Fitting polynomial $P_k(x) = \sum_{j=0}^k a_j x^j$ of order $k$:
$$(A^T W A) \mathbf{a} = A^T W \mathbf{y}$$
- $A \in \mathbb{R}^{N \times (k+1)}$ is the Vandermonde design matrix: $A_{ij} = x_i^j$.
- $M = A^T W A$ is the **Hankel normal matrix**:
  $$M_{jl} = \sum_{i=1}^N w_i x_i^{j+l}, \quad b_j = \sum_{i=1}^N w_i y_i x_i^j$$
Solved directly using Doolittle LU decomposition without computing matrix inverse.

---

### 10.5 Matrix Ill-Conditioning & Shifted Chebyshev Basis
As polynomial degree $k$ increases, monomial powers $1, x, x^2, x^3, \dots$ become nearly linearly dependent on $[0, 1]$.
- **Frobenius Condition Number:** $\kappa(M) = \|M\|_F \|M^{-1}\|_F$. For standard monomials with $k=3$, $\kappa(M) \sim 10^4 - 10^6$.
- **Shifted Chebyshev Polynomials:** Defined on $[0, 1]$ via $T_m^*(x) = T_m(2x - 1)$:
  $$\phi_0(x) = 1, \quad \phi_1(x) = 2x - 1, \quad \phi_2(x) = 8x^2 - 8x + 1, \quad \phi_3(x) = 32x^3 - 48x^2 + 18x - 1$$
- **Recurrence Relation:**
  $$\phi_{m+1}(x) = 2(2x - 1) \phi_m(x) - \phi_{m-1}(x)$$
- **Condition Number Benefit:** The normal matrix formed using Chebyshev basis functions has condition number orders of magnitude lower than the monomial Hankel matrix, preventing catastrophic loss of precision.

## Appendix A — Formula Cheat Sheet

### Random

- LCG: $x_{i+1}=(ax_i+c)\bmod m$; $[0,1)$ via $x_i/m$.
- Hull–Dobell: $c\ne0$; $c$ coprime $m$; $a-1$ multiple of all $p|m$; $a-1$ multiple of 4 if $4|m$.
- Park–Miller: $x_{i+1}=ax_i\bmod m$.
- Uniform $[a,b)$: $u=a+(b-a)x$. Exponential: $y=-\frac1a\ln x$.
- Correlation: $\epsilon(n,N)=\frac1N\sum x_ix_{i+n}-\frac1N\sum x_i\cdot\frac1N\sum x_{i+n}\to 0$.
- $\pi\approx 4\times \mathrm{inside}/\mathrm{trials}$; decay $P=\lambda\Delta t$.

### Gauss–Jordan

- $[A|b]\to[I|\tilde b]$ via row ops; pivot normalize + column elimination.
- $\det A=(-1)^n\prod a'_{ii}$; $[A|I]\to[I|A^{-1}]$.
- Memory $N(N+1)$ solve / $N(2N)$ inverse; $\sim N=10$ practical limit.

### LU

- Doolittle ($\ell_{ii}=1$): $u_{ij}=a_{ij}-\sum_{k<i}\ell_{ik}u_{kj}$ ($i\le j$); $\ell_{ij}=(a_{ij}-\sum_{k<j}\ell_{ik}u_{kj})/u_{jj}$ ($i>j$).
- Solve $Ly=b$ forward $y_i=b_i-\sum_{j<i}\ell_{ij}y_j$, then $Ux=y$ backward $x_i=(y_i-\sum_{j>i}u_{ij}x_j)/u_{ii}$.
- $\det A=(-1)^n\prod u_{ii}$. Cost $O(N^3)$, half memory via overwrite.
- Cholesky ($A=LL^T$, pos-def): $\ell_{ii}=\sqrt{a_{ii}-\sum_{j<i}\ell_{ij}^2}$, $\ell_{ij}=(a_{ij}-\sum_k\ell_{ik}\ell_{jk})/\ell_{ii}$; ~2× faster than LU.

### Iterative

- Jacobi: $x^{(k+1)}_i=(b_i-\sum_{j\ne i}a_{ij}x^{(k)}_j)/a_{ii}$. Needs 2 vectors. Converges if strictly DD: $|a_{ii}|>\sum_{j\ne i}|a_{ij}|$.
- Gauss–Seidel: $x^{(k+1)}_i=(b_i-\sum_{j<i}a_{ij}x^{(k+1)}_j-\sum_{j>i}a_{ij}x^{(k)}_j)/a_{ii}$. One vector. Converges if strictly DD or symm. pos-def.
- SOR: $x^{(k+1)}_i=(1-\omega)x^{(k)}_i+\omega\cdot(\text{GS update})$; $\omega=1$ GS; $1<\omega<2$ over; $0<\omega<1$ under.

### Roots

- Bisection: $c=(a+b)/2$; $|b_n-a_n|=|b-a|/2^n\le\epsilon$; bracket $f(a)f(b)<0$.
- Regula falsi: $c_n=b_n-(b_n-a_n)f(b_n)/[f(b_n)-f(a_n)]$.
- Fixed-point: $x_{n+1}=g(x_n)$; converges if $|g'|<1$.
- Newton: $x_{n+1}=x_n-f(x_n)/f'(x_n)$; secant $f'\approx[f(x+h)-f(x-h)]/2h$.
- Multi-NR: ${\bf x}^{(i+1)}={\bf x}^{(i)}-J^{-1}{\bf f}$; $J_{pq}=\partial f_p/\partial x_q$.
- Laguerre: $G=P'/P$, $H=G^2-P''/P$, $a=n/[G\pm\sqrt{(n-1)(nH-G^2)}]$, $\beta_{k+1}=\beta_k-a$; then synthetic division.

### Integration

- Midpoint: $M_N=\sum h\,f(x_n)$, $h=(b-a)/N$, error $\le(b-a)^3|f''|_{\max}/24N^2$.
- Trapezoidal: $T_N=\frac h2[f_0+2f_1+\cdots+2f_{N-1}+f_N]$, error $\le(b-a)^3|f''|_{\max}/12N^2$.
- Simpson: $S_N=\frac h3[f_0+4f_1+2f_2+\cdots+4f_{N-1}+f_N]$, error $\le(b-a)^5|f''''|_{\max}/180N^4$; $N$ even; $S_{2N}=\frac23M_N+\frac13T_N$.
- Monte Carlo: $F_N=\frac{b-a}{N}\sum f(X_i)$; $X=a+(b-a)\xi$; $\sigma\sim1/\sqrt N$.
- Gauss–Legendre on $[a,b]$: $\frac{b-a}{2}\sum w_n f(\frac{b-a}{2}x_n+\frac{b+a}{2})$; exact to degree $2n-1$.

### ODEs

- Explicit Forward Euler: $y_{n+1}=y_n+h f(x_n,y_n)$; $O(h)$ global; conditional stability $|1+h\lambda|\le 1$.
- Implicit Backward Euler: $y_{n+1}=y_n+h f(x_{n+1},y_{n+1})$; solve via Newton $g(y)=y-y_n-hf(x_{n+1},y)=0$; unconditionally stable.
- Heun's PC: $y_{n+1}^p=y_n+hf_n$, $y_{n+1}=y_n+\frac h2[f_n+f(x_{n+1},y_{n+1}^p)]$; $O(h^2)$ global.
- RK4: $k_1=f(x_n,y_n)$, $k_2=f(x_n+h/2,y_n+hk_1/2)$, $k_3=f(x_n+h/2,y_n+hk_2/2)$, $k_4=f(x_n+h,y_n+hk_3)$; $y_{n+1}=y_n+\frac h6(k_1+2k_2+2k_3+k_4)$; $O(h^4)$.
- Shooting method slope interpolation: $z^*=z_1+(z_2-z_1)[\beta-y(b;z_1)]/[y(b;z_2)-y(b;z_1)]$.
- 1D FEM Poisson: $K_{ii}=2/h, K_{i,i\pm 1}=-1/h$; solve $K y = F$ via Cholesky $K=LL^T$.

### PDEs

- 1D Heat FTCS: $u_i^{j+1}=\alpha(u_{i+1}^j+u_{i-1}^j)+(1-2\alpha)u_i^j$; $\alpha=h_t/h_x^2\le 0.5$ stability.
- 1D Heat BTCS: $-\alpha u_{i-1}^{j+1}+(1+2\alpha)u_i^{j+1}-\alpha u_{i+1}^{j+1}=u_i^j$; unconditionally stable; TDMA in $O(N)$.
- Crank-Nicolson: $-\frac\alpha2 u_{i-1}^{j+1}+(1+\alpha)u_i^{j+1}-\frac\alpha2 u_{i+1}^{j+1}=\frac\alpha2 u_{i-1}^j+(1-\alpha)u_i^j+\frac\alpha2 u_{i+1}^j$; $O(h_t^2+h_x^2)$.
- 1D Wave explicit: $u_i^{j+1}=2(1-r^2)u_i^j+r^2(u_{i+1}^j+u_{i-1}^j)-u_i^{j-1}$; CFL $r=ch_t/h_x\le 1$.
- Thomas algorithm (TDMA): $c'_i=c_i/(b_i-a_i c'_{i-1})$, $d'_i=(d_i-a_i d'_{i-1})/(b_i-a_i c'_{i-1})$; $x_i=d'_i-c'_i x_{i+1}$.

### Eigenvalues

- Power iteration: $x_{k+1}=A x_k$; Rayleigh quotient $\lambda_1=(A x_k \cdot y)/(x_k \cdot y)$; normalize $x_{k+1}/\|x_{k+1}\|_2$; error scales as $|\lambda_2/\lambda_1|^k$.
- Deflation: $A'=A-\lambda_1 u_1 u_1^T$ where $u_1=v_1/\|v_1\|_2$.
- Tridiagonal spectrum: $\lambda_k=b+2\sqrt{ac}\cos[k\pi/(n+1)]$, $v_{k,i}=\sin[ik\pi/(n+1)]$.

### Data Fitting

- Lagrange interpolation: $P_N(x)=\sum y_i \prod_{k\ne i}(x-x_k)/(x_i-x_k)$.
- Weighted linear regression: $m=(S S_{xy}-S_x S_y)/\Delta$, $c=(S_{xx} S_y-S_x S_{xy})/\Delta$; $\Delta=S S_{xx}-S_x^2$; $r^2=(S S_{xy}-S_x S_y)^2/[\Delta(S S_{yy}-S_y^2)]$.
- Polynomial normal equations: $(A^TW A)a=A^TWy$; Hankel matrix $M_{jl}=\sum w_i x_i^{j+l}$; condition number $\kappa(M)=\|M\|_F \|M^{-1}\|_F$.
- Shifted Chebyshev basis on $[0,1]$: $\phi_0=1, \phi_1=2x-1, \phi_{m+1}=2(2x-1)\phi_m-\phi_{m-1}$.

---

## Appendix B — All Exercises / Tryouts in One List

1. **pRNG:** logistic map $x_0=0.1$, $c=3.2,3.7,3.98$; plot $x_i$ vs $x_{i+5}$; correlation tests.
2. **LCG:** ugly ($27,11,54,10$), bad ($57,1,256,10$), good ($1103515245,12345,32768,10$) — generate, histogram, compare.
3. **$\pi$ by Monte Carlo** + **radioactivity** $A\to B\to C$ simulation.
4. **Distributions:** uniform $[a,b)$, exponential $y=-\ln(x)/a$; histogram 2000 draws.
5. **Gauss–Jordan pivoting tryout:** 4-equation $y+z-2w$ system — first-column pivot only.
6. **Gauss–Jordan full:** $2y+5z=1$ system → $(-2,-2,1)$; determinant + inverse practice.
7. **LU:** Doolittle by hand for $3\times3$; forward/backward solve; determinant; Cholesky example check.
8. **Jacobi/GS:** diagonally-dominant swap example; non-DD $2x+3y+z$ system; SOR tridiagonal; 4-variable exercise in §4.4.
9. **Roots:** bracket + bisect/regulafalsi any of $\cos x-x^3$, $3x+\sin x-e^x$, $xe^x-2$; fixed-point forms of $x^2-x-2$; devotee-flowers $f(x)=2x-X$; multivariable $(2,3)$ system via both fixed-point and Newton; Laguerre + synthetic division on $-x^3+3x^2-4$, $6x^3-11x^2-26x+15$, $x^4-x^3-7x^2+x+6$.
10. **Integration:** $N=4,8,12,20$ Midpoint vs Trapezoidal on $\int_1^2 dx/x$, $\int_0^{\pi/2}x\cos x\,dx$; $\int_0^1x^2dx$ by hand ($M_4=21/64$); $N$ from error bounds ($29/41/7$ for $1/x$); Monte Carlo convergence plot; Gauss–Legendre $\int_{-1}^1xe^xdx$ with $N=4,5,6$ vs Simpson.

---

---

## Appendix C — Assignments (Repo): Questions & How They're Solved

> Source repo: [Aryans-lab/Computational_Physics_Lab](https://github.com/Aryans-lab/Computational_Physics_Lab)
> This appendix adds **every assignment question** from the repo and **how it was solved** — method, key code idea, inputs, and actual outputs.
> Folder names below match the repo exactly. Library convention: from Week 01 onward every routine lives in `mylib.py`; drivers only do I/O + calls + verification.

### How to use this appendix in exams

- **Q1-type "write library routine" questions** → memorize the 5–10 line core loop of each method (given under "Core logic").
- **Q2+ "apply + verify" questions** → memorize input format, verification (`||Ax−b||`, residual plots, iteration counts), and the final numbers as sanity checks.
- **Plots** → in `assign1_rand` and `assign6b_simmc` the plot *is* the answer (correlation scatter, π convergence, decay curves, MC convergence).

---

### C.0 Warmup — `Assgn0/` (= `assignments/week00_warmup/`)

#### C.0.1 Q1 — Loops & arithmetic

**Question:** Sum of first 20 even numbers and factorial of 8, using a loop.

**How solved:**

```python
sum_even = sum(2*i for i in range(1, 21))   # 2+4+...+40
fact = 1
for i in range(1, 9): fact *= i             # 8!
```

**Output:**

```
Sum of first 20 even numbers: 420
Factorial of 8: 40320
```

**Note:** `week00` version reads `N M` from `data/q1_input.txt` (`20 8`) and writes `output/q1_output.txt` (prof's Rule 2+3).

#### C.0.2 Q2 — GP & HP sums (no closed forms!)

**Question:** Sum $N=15$ terms of GP and HP with $t_0=1.25$, $r=0.5$, $d=1.5$. Analytical formulae **not allowed** — accumulate term by term.

**How solved:**

```python
# GP: t, t*r, t*r^2, ...
term = t0
for i in range(n): sum_gp += term; term *= r
# HP: 1/t0 + 1/(t0+d) + 1/(t0+2d) + ...
term = t0
for i in range(n): sum_hp += 1/term; term += d
```

**Output:**

```
Sum of GP series: 2.4999237060546875   # → t0/(1-r) = 2.5, truncation err ≈ 7.8e-5
Sum of HP series: 2.4139570733659186
```

#### C.0.3 Q3 — Matrix I/O: find $AB$, $BC$, $D\cdot C$

**Question:** Read ASCII files `asgn0_matA`, `asgn0_matB`, `asgn0_vecC`, `asgn0_vecD`, compute $AB$, $BC$, $D\cdot C$.

**Inputs:**

```
A =  2.0  -3.0   1.4      B =  0.0  -1.0   1.0
     2.5   1.0  -2.0           1.5   0.5  -2.0
    -0.8   0.0   3.1           3.0   0.0  -2.0
C = [-2, 0.5, 1.5]^T    D = [1, 0, -1]^T
```

**How solved:** hand-written triple loop `result[i][j] += X[i][k]*Y[k][j]`; dot product `Σ X[i][0]*Y[i][0]`. No numpy.

**Output:**

```
Matrix AB:
[-0.3000000000000007, -3.5, 5.2]
[-4.5, -2.0, 4.5]
[9.3, 0.8, -7.0]
Matrix BC:
[1.0]
[-5.75]
[-9.0]
Dot product of D.C is: -3.5
```

**Exam note:** `AB[0][0]` shows ~1e-16 float error — classic IEEE-754 artefact, mentioned in repo.

#### C.0.4 Q5 — `MyComplex` class

**Question:** Define `MyComplex` class; sum, difference, product, modulus of $(1.3-2.2j)$ and $(-0.8+1.7j)$.

**How solved:** class with `real, image`; `add/subtract/multiply/modulus`:

```python
# (a+bi)(c+di) = (ac-bd) + (ad+bc)i ; |c| = sqrt(a²+b²)
```

**Output:**

```
Sum:        0.5 - 0.5j
Difference: 2.1 - 3.9j
Product:    2.7 + 3.97j
|c1| = 2.5553864678
|c2| = 1.8788294228
```

---

### C.1 Random Numbers — `assign1_rand/`

Library: `lcg.py` / `mylib.py::myrand` — LCG with **gcc params** $a=1103515245$, $c=12345$, $m=32768$:

```python
def myrand(seed=0):
    if not hasattr(myrand, 'x0'): myrand.x0 = 0
    a, c, m = 1103515245, 12345, 32768
    if seed: myrand.x0 = seed
    x1 = (a*myrand.x0 + c) % m
    rn = x1/m            # → [0,1)
    myrand.x0 = x1
    return rn
```

#### C.1.1 Q1 — Logistic-map "pRNG" + correlation plots

**Question:** Use $x_{i+1}=c\,x_i(1-x_i)$ to generate 1,000 numbers. Show correlation by plotting $x_i$ vs $x_{i+k}$. Seed $x_0=0.1$; choose 3 different $c$; try $k=3,5,10$.

**How solved:** repo chose $c\in\{3.85, 3.95, 4.0\}$; $3\times3$ grid (`c` × `k`), `scatter(s=2)`, saved as `corrolation_question1.png`.

**What to conclude (viva point):** structured bands/patterns for smaller $c$ = correlated/bad; $c=4.0$ fills space more uniformly (chaotic regime). This motivates LCG + formal tests (§1.3).

#### C.1.2 Q2 — Own LCG + correlation check ($k=5$)

**Question:** Write own LCG with $a=1103515245$, $c=12345$, $m=32768$; check correlation by plotting for $k=5$.

**How solved:** `n=2000`, `x0=1.2`, scatter $x_i$ vs $x_{i+5}$ → `correlation_Question2.png`. Uniform cloud with no bands = good (vs Q1 patterns).

#### C.1.3 Q3 — $\pi$ by throwing method

**Question:** Determine $\pi$ by throwing method with quarter circle of unit radius in first quadrant.

**How solved:**

```python
# N throws: x,y ~ U[0,1); inside if x²+y² ≤ 1
pi_est = 4 * (inside / N)
```

$N=20\ldots5000$; two-panel plot (estimate vs $N$ + residual vs $N$) → `pi_estimate_residual.png`.

**Output:** `output_Q3.txt`:

```
Final estimated value of π at N=5000: 3.1216
```

**Viva point:** error shrinks $\sim 1/\sqrt{N}$ (slow!); residual plot oscillates around 0.

#### C.1.4 Q4 — Radioactive decay $A\to B\to C$

**Question:** Simulate $A\to B\to C$ with $\lambda_A=0.0347$, $\lambda_B=0.0198$, $N_A=500$, $N_B=N_C=0$; overlay plot of $N_A,N_B,N_C$ vs time. **No libraries** (own pRNG).

**How solved:**

```python
# per step dt: each A-nucleus decays if r < λA*dt; each B if r < λB*dt
# NA -= decays_A ; NB += decays_A - decays_B ; NC += decays_B
```

Repo: $t_{\max}=400$, $dt=0.1$ (4000 steps); pre-generated pool `NA0*n_steps` randoms; scatter `s=2` → `radioactive_simulation.png`.

**Expected curves:** $N_A$ exponential fall; $N_B$ rises then falls (transient equilibrium); $N_C$ saturating growth. Governs $P=\lambda\Delta t$ (§1.8).

#### C.1.5 Q5 — Exponential deviates from uniform

**Question:** Generate $\geq 5{,}000$ exponential $e^{-x}$ deviates from $U[0,1)$ pRNG.

**How solved:** inverse transform (§1.9) $y=-\ln(x)$, $x\sim U[0,1)$; histogram 40 bins → `exponential_distribution.png`:

```python
y_rn = [-log(myrand()) for _ in range(5000)]
plt.hist(y_rn, bins=40, color='skyblue', edgecolor='black')
```

---

### C.2 Direct Solvers — `assign2_gj_lu/`

#### C.2.1 Q1 — Library: Gauss–Jordan (row pivoting) + LU (Doolittle)

**Question:** Prepare library routines for Gauss–Jordan with row pivoting and LU decomposition; state Doolittle vs Crout clearly.

**Choice:** **Doolittle** ($L_{ii}=1$).

**Core logic — GJ (`gauss_jordan_elimination_augmented`):**

```python
# for each col i: pick max|.| pivot row below (partial pivot), swap,
# normalize pivot row, eliminate all other rows k≠i
# raise on all-zero row: no-solution vs infinite-solutions
```

**Core logic — LU (`lu_decomposition`):**

```python
# U[0] = A[0]; for i: U[i][j≥i] = A[i][j]-Σ L[i][k]U[k][j]
# L[j>i][i] = (A[j][i]-Σ L[j][k]U[k][i]) / U[i][i]
```

#### C.2.2 Q2 — GJ check on $3\times3$ system

**Question:** Solve via Gauss–Jordan: $2y+5z=1$, $3x-y+2z=-2$, $x-y+3z=3$ (matrix in `asgn2_mat1`).

**Input `asgn2_mat1`** (augmented):

```
0  2   5 |  1
3 -1   2 | -2
1 -1   3 |  3
```

**How solved:** `read_matrix_from_file` → `gauss_jordan_elimination_augmented` → print. Output (`output2.txt`):

```
# The values of x, y, z are: -2.0 -2.0 1.0
```

Matches hand calculation in §2.6.

#### C.2.3 Q3 — LU decomposition + verify

**Question:** LU-decompose $A=\begin{pmatrix}1&2&4\\3&8&14\\2&6&13\end{pmatrix}$ and verify $L\cdot U=A$.

**How solved:** `lu_decomposition` → `matrix_multiply(L,U)` → compare. Output (`output3.txt`):

```
L = [1,0,0]  U = [1,2,4]
    [3,1,0]      [0,2,2]
    [2,1,1]      [0,0,3]
L·U = original A ✓
```

---

### C.3 LU Solve + Cholesky — `assign3_lu_chsk/`

#### C.3.1 Q1 — Library: `lu_forback` + `cholesky_decomposition`

**Question:** Library routines for Cholesky (must check symmetry) + LU forward-backward solve.

**Core logic:**

```python
# lu_forback: y[i] = d[i]-Σ_{j<i} L[i][j]y[j]      (forward)
#             x[i] = (y[i]-Σ_{j>i} U[i][j]x[j])/U[i][i]  (backward)
# cholesky: check square + symmetric (A[i][j]==A[j][i]) else raise
# L[i][j] = (A[i][j]-Σ L[i][k]L[j][k]) ; diag = sqrt(value), require value>0 (pos-def)
```

#### C.3.2 Q2 — Solve $6\times6$ via LU forward-backward

**Question:** Solve the given $6\times6$ system by LU + forward-backward substitution (Doolittle or Crout).

**Inputs** (`asgn3_mat1`, `asgn3_vec1`):

```
A = 1 -1  4  0  2  9      b = [19, 2, 13, -7, -9, 2]^T
    0  5 -2  7  8  4
    1  0  5  7  3 -2
    6 -1  2  3  0  8
   -4  2  0  5 -5  3
    0  7 -1  5  4 -2
```

**How solved:** `lu_forback('asgn3_mat1','asgn3_vec1')`. Output (`output2.txt`):

```
a1 = -1.761817043997862
a2 =  0.896228033874014
a3 =  4.051931404116158
a4 = -1.6171308025395423
a5 =  2.041913538501913
a6 =  0.15183248715593523
```

Verify with $||Ax-b||_2\approx 0$.

#### C.3.3 Q3 — Cholesky factorize + solve

**Question:** Cholesky-factorize $A$ and solve $Ax=b$.

**How solved:** `L = cholesky_decomposition(A)` → `x = cholesky_forback(L, b)` (i.e. $Ly=b$, then $L^Tx=y$). Output (`output3.txt`):

```
L = [2.0, 0, 0, 0]
    [0.5, 1.6583, 0, 0]
    [0.5, -0.7538, 1.0871, 0]
    [0.5, 0.4523, 0.0836, 1.2403]
Solution: a1 ≈ 0 (−1.1e−16), a2 = 1.0, a3 = 1.0, a4 = 1.0
```

(Clean integers → good self-check that Cholesky path is right.)

---

### C.4 Iterative Solvers — `assign4_jac_gs/`

#### C.4.1 Q1 — Library: Jacobi, Gauss–Seidel, SOR

**Question:** Add Jacobi, Gauss–Seidel, SOR to `mylib.py`.

**Core logic (all three):**

```python
# optional row swap toward diagonal dominance if |A[i][i]| < Σ_{j≠i}|A[i][j]|
# Jacobi:      x[i] = (b[i]-Σ_{j≠i} A[i][j]x_prev[j])/A[i][i]   (uses old vector)
# Gauss-Seidel:x[i] = (b[i]-Σ_{j<i}A[i][j]x[j]-Σ_{j>i}A[i][j]x[j])/A[i][i]  (in-place)
# SOR:         x[i] = ω*(GS update) + (1-ω)*x_prev[i]
# stop when all |x[i]-x_prev[i]| < tol
```

#### C.4.2 Q2 — Jacobi vs Gauss–Seidel on $6\times6$

**Question:** Solve with Jacobi and Gauss–Seidel (files `asgn4_mat1`, `asgn4_vec1` at repo root).

**Inputs:**

```
A = 4 -1  0 -1  0  0      b = [2,1,2,2,1,2]^T
   -1  4 -1  0 -1  0
    0 -1  4  0  0 -1
   -1  0  0  4 -1  0
    0 -1  0 -1  4 -1
    0  0 -1  0 -1  4
```

**How solved:** `jacobi_it(A,b,tol=1e-6)`, `gauss_seidel(A,b,tol=1e-6)`. Output:

```
Converged after 16 iterations
Jacobi:      [0.99999975, 0.99999979, 0.99999991, 0.99999985, 0.99999987, 0.99999995]
Converged after 16 iterations
Gauss-Seidel:[same to ~1e-16]
```

True solution is $[1,1,1,1,1,1]^T$. (Matrix is strictly DD → both converge.)

#### C.4.3 Q3 — SOR ($\omega=1.57$, $n=10$) vs Gauss–Seidel

**Question:** Solve tridiagonal $A_{ii}=2$, $A_{i,i\pm1}=-1$, $b_i=1$ ($n=10$) by Gauss–Seidel and SOR with $\omega=1.57$; compare convergence rate.

**How solved:** build $A$ in code → `sor_gauss_seidel(A,b,omega=1.57,tol=1e-16)`. Output:

```
Converged after 71 iterations
SOR solution: [5.0, 9.0, 12.0, 14.0, 15.0, 15.0, 14.0, 12.0, 9.0, 5.0]
```

**Viva point:** symmetric parabola (max 15 at center) — expected for this discrete-Poisson-like system; SOR $\approx$ 2–6× fewer iterations than plain GS (slides quote ~150 → ~25 for looser tol).

---

### C.5 Bisection & Regula Falsi — `assign5_bi_rf/`

#### C.5.1 Q1 — Library: `bisection` + `regula_falsi`

**Question:** Code library functions Bisection and Regula Falsi.

**Core logic (both):**

```python
# auto-bracket: while f(a)*f(b) ≥ 0: expand side with smaller |f| by m*(b-a)
# bisection: c=(a+b)/2 ; regula falsi: c=(a*f(b)-b*f(a))/(f(b)-f(a))
# keep half with sign change; stop when |f(c)|<acc AND |b-a|<acc
```

#### C.5.2 Q2 — Root of $\log(x/2)-\sin(5x/2)$

**Question:** Root in $[1.5,3.0]$ (header) to $10^{-6}$ via Bisection + Regula Falsi. *(Code actually starts from $[0.5,1.0]$ and auto-expands — same root found.)*

$$f(x)=\log(x/2)-\sin(5x/2)$$

**How solved:** `bisection(f,0.5,1.0,1e-6)`, `regula_falsi(f,0.5,1.0,1e-6)`. Output (`output2.txt`):

```
Bisection:   root ≈ 1.4019301055908207,  22 iterations
Regula Falsi: root ≈ 1.4019299316146128,   9 iterations
```

**Takeaway:** RF ~2.4× fewer iterations here (interpolation pays off on smooth $f$).

#### C.5.3 Q3 — Auto-bracket $-x-\cos x$ from $[2,3]$, then solve

**Question:** For $f(x)=-x-\cos x$, find a valid bracket starting from $[2,3]$, then find root.

**How solved:** expand-toward-smaller-$|f|$ loop → bracket → `regula_falsi`. Output (`output3.txt`):

```
Bracket found: [-2.0625, 5.53125]
Root ≈ -0.7390851332151607,  22 iterations
```

($-0.739085\ldots$ is the famous Dottie number, root of $x+\cos x=0$.)

---

### C.6 Fixed-Point & Newton–Raphson — `assign5_fx_nr/`

#### C.6.1 Q1 — Library: `fixed_point`, `newton_raphson` (1D), `newton_raphson_system`

**Question:** Library functions for Fixed-point and Newton–Raphson.

**Core logic:**

```python
# fixed_point: x_next = g(x); stop |x_next-x| < acc
# newton 1D: x_next = x - f(x)/f'(x); f' analytic if given else (f(x+h)-f(x-h))/2h
#            guard |f'| < 1e-12; stop |x_next-x|<acc AND |f(x_next)|<acc
# newton system: x -= J^{-1}·F ; J^{-1} via gauss_jordan_inverse
```

#### C.6.2 Q2 — $3x+\sin x-e^x=0$: Bisection vs RF vs NR

**Question:** Root in $[-1.5,1.5]$ to $10^{-6}$ via Bisection, Regula Falsi, Newton–Raphson; compare convergence. $f(x)=3x+\sin x-e^x$, $f'(x)=3+\cos x-e^x$, NR $x_0=0$.

**How solved:** all three from `mylib`; NR run twice (analytic $f'$ + numerical $f'$). Output (`question2_out.txt`):

```
Bisection:    0.3604220151901245,  23 iterations
Regula Falsi: 0.3604217029603244,  25 iterations
Newton (analytic):  0.36042168047601975,  3 iterations
Newton (numeric):   0.3604216804760228,   3 iterations
```

**Takeaway (the money line for viva):** NR converges in **3 iterations** vs ~24 — quadratic convergence; analytic vs numeric derivative agree to 15 digits.

#### C.6.3 Q3 — Fixed-point for $x^2-2x-3=0$

**Question:** Root via fixed-point. (Analytical roots $-1$, $3$.)

**How solved:** rewrite for root near 3 as $x=g(x)=\sqrt{2x+3}$ ($|g'(3)|=1/3<1$ ✓), $x_0=1$:

```python
g1 = lambda x: (2*x+3)**0.5
fixed_point(g1, x0=1, accuracy=1e-6)
```

**Output:** `root ≈ 2.9999998290170797, 15 iterations`.

**Viva point:** choice of $g$ matters — this $g$ converges to 3; a different $g$ would be needed for $-1$ (cf. §5.4).

#### C.6.4 Q4 — Nonlinear system via multivariate Newton

**Question:** Solve $f_1=3x^2-2xy-3=0$, $f_2=3y^2-4xy=0$ via Newton–Raphson or fixed-point.

**How solved:** multivariate NR with **numerical Jacobian** (central differences) + GJ inverse, $x^{(0)}=(2.5,3.1)$:

```python
# J[i][j] ≈ (fi(x+h·ej)-fi(x-h·ej))/2h ; x -= J^{-1}·F
```

**Output:** `Final root approximation: (3.0000, 4.0000) after 4 iterations.`

Check: $3·9-2·3·4-3=27-24-3=0$ ✓; $3·16-4·3·4=48-48=0$ ✓.

---

### C.7 Polynomials — `assign5_poly/`

#### C.7.1 Q1 — Library: Laguerre + deflation (real roots, coeffs only)

**Question:** Library functions for Laguerre's method + deflation, real roots only; pass only coefficients (descending powers $[a_n,\dots,a_0]$).

**Core logic:**

```python
# Horner: P(x) = (...((a_n)x + a_{n-1})x + ... )
# P', P'' coeff lists by power rule
# G=P'/P, H=G²-P''/P, D=(n-1)(nH-G²)
# a = n/(G ± √D)  [larger-|.| denominator]; b -= a
# stop |b_new-b|<acc1 AND |P(b_new)|<acc2
# synthetic_division: bring-down + multiply-add; remainder must be ~0
# laguerre_roots: loop Laguerre→deflate; finish linear ax+b=0 directly
```

#### C.7.2 Q2 — All real roots of three polynomials

**Question:** Laguerre + deflation for:

$$P_1=x^4-x^3-7x^2+x+6,\quad P_2=x^4-5x^2+4,\quad P_3=2x^5-19.5x^3+0.5x^2+13.5x-4.5$$

Coefficients: `P1=[1,-1,-7,1,6]`, `P2=[1,0,-5,0,4]`, `P3=[2,0,-19.5,0.5,13.5,-4.5]`; $b_0=1$. Output (`output2.txt`):

```
P1 roots:  1.0 (1 it), -1.0 (5 it), -2.0 (2 it), 3.0 (linear finish)
P2 roots:  1.0 (1 it), 2.0 (3 it), -1.0 (2 it), -2.0 (linear finish)
P3 roots:  0.50001747 (9 it), 0.49998253 (2 it), -1.0 (4 it), -3.0 (2 it), 3.0 (finish)
```

**Viva points:** $P_1$ matches slide answer $\{1,-2,3,-1\}$ (§5.7); $P_3$ has a **double root at 0.5** (splits into $0.50002/0.49998$ numerically — typical deflation behavior on multiple roots).

---

### C.8 Midpoint & Trapezoidal — `assign6_midtrap/`

#### C.8.1 Q1 — Library: `midpoint` + `trapezoidal`

**Question:** Library functions for Midpoint and Trapezoidal integration.

**Core logic:**

```python
# midpoint: h=(b-a)/N; Σ f(a+(i+0.5)h); × h
# trapezoidal: h=(b-a)/N; (h/2)[f(a)+f(b)+2Σ_{i=1}^{N-1} f(a+ih)]
```

#### C.8.2 Q2 — Three integrals, $N=4,8,15,20$, tabular

**Question:** Midpoint + Trapezoidal for $N=4,8,15,20$; tabular comparison:

$$I_1=\int_1^2\frac{dx}{x}\;(=0.69314718),\quad I_2=\int_0^{\pi/2}x\cos x\,dx\;(=\pi/2-1\approx0.57079633),\quad I_3=\int_0^1x\tan^{-1}x\,dx\;(=\pi/4-1/2\approx0.28539816)$$

**How solved:** loop over $N$, both rules, `:.8f` table. Output (`output2.txt`):

```
Function 1: 1/x from 1 to 2
N     Midpoint      Trapezoidal
4     0.69121989    0.69702381
8     0.69266055    0.69412185
15    0.69300843    0.69342480
20    0.69306910    0.69330338

Function 2: x·cos(x) from 0 to pi/2
N     Midpoint      Trapezoidal
4     0.58744792    0.53760713
8     0.57493427    0.56252752
15    0.57197166    0.56844624
20    0.57145729    0.56947459

Function 3: x·tan⁻¹(x) from 0 to 1
N     Midpoint      Trapezoidal
4     0.28204605    0.29209835
8     0.28456102    0.28707220
15    0.28516010    0.28587426
20    0.28526426    0.28566596
```

**Takeaway:** Midpoint beats Trapezoidal at every $N$ here (matches slide warning §6.4); both crawl toward analytic values as $N$ grows.

---

### C.9 Simpson & Monte Carlo — `assign6b_simmc/`

#### C.9.1 Q1 — Library: `simpson` + `monte_carlo`

**Question:** Library function for Simpson's 1/3 rule and Monte Carlo (function-passing style).

**Core logic:**

```python
# simpson: require N even; (h/3)[f(a)+f(b)+4Σodd+2Σeven]; return (integral, evals=N+1)
# monte_carlo: seed own myrand; x=a+(b-a)·ξ; mean_f, σf²=mean(f²)-mean(f)²
#   integral=(b-a)·mean_f ; σ_int=(b-a)·σf/√N ; return all three
```

#### C.9.2 Q2 — 6-decimal accuracy via error-bound $N$

**Question:** Midpoint + Simpson accurate to 6 decimals; choose $N$ from error-bound formulae. $I_1=\int_1^2 dx/x$, $I_2=\int_0^{\pi/2}x\cos x\,dx$.

**How $N$ chosen (the actual exam skill):** invert §6.6 bounds for target $10^{-6}$:

- $I_1$: $|f''|_{\max}=2$ on $[1,2]$ → Midpoint $N=289$; $|f''''|_{\max}=24$ → Simpson $N=20$.
- $I_2$: bounds give Midpoint $N=610$, Simpson $N=22$.

**How solved:** `midpoint(f,·,N)`, `simpson(f,·,N)`. Output (`output2.txt`):

```
Integral of 1/x from 1 to 2
Midpoint = 0.6931468064035278
Simpson  = 0.6931473746651161  (21 evaluations)

Integral of x·cos(x) from 0 to pi/2
Midpoint = 0.5707970370864699
Simpson  = 0.570796987316687   (23 evaluations)
```

**Takeaway:** Simpson needs **~15–25× fewer points** than Midpoint for the same 6 decimals ($O(h^4)$ vs $O(h^2)$).

#### C.9.3 Q3 — Monte Carlo $\int_{-1}^{1}\sin^2x\,dx$ to 4 decimals + plots

**Question:** Monte Carlo (own pRNG) to 4 decimals; plot integral-vs-$N$ and deviation-vs-$N$.

$$I=\int_{-1}^{1}\sin^2x\,dx = 1-\frac{\sin 2}{2}\approx 0.54535129$$

**How solved:** $N=1100,1200,\dots,50000$ (step 100), **fresh seed per run**; record integral, $|\text{est}-\text{exact}|$, $\sigma_{\text{int}}$; find first $N$ with deviation $<10^{-4}$; save `integral_vs_N.png`, `deviation_vs_N.png`.

**Viva points:** integral oscillates around exact value with shrinking envelope $\sim 1/\sqrt{N}$; deviation plot crosses the $10^{-4}$ line back and forth (MC error is statistical, not monotonic — unlike deterministic rules).

---

### C.10 Gaussian Quadrature — `assign6c_gq/`

#### C.10.1 Q1 — Library: unified Gauss–Legendre + Gauss–Laguerre ($n=1\ldots6$)

**Question:** One Gauss quadrature routine with hardwired zeros/weights tables (at least $n=5$); dispatch on Legendre vs Laguerre.

**How solved:** dicts `GAUSS_LEGENDRE[N] = (nodes, weights)`, `GAUSS_LAGUERRE[N] = (nodes, weights)` for $N=1\ldots6$; single `gaussian_quadrature(f,a,b,N,method)`:

```python
# legendre: x = ((b-a)/2)t + (b+a)/2 ; integral = ((b-a)/2)·Σ wᵢf(xᵢ)
# laguerre: integral = Σ wᵢf(tᵢ)   [f EXCLUDES e^{-t}; weight built into wᵢ]
```

#### C.10.2 Q2 — 4-pt Legendre: $\int_{-1}^{1}x^2/(1+x^4)\,dx$ vs Simpson

**Question:** 4-point GQ for $\int_{-1}^{1}\frac{x^2}{1+x^4}dx\approx 0.487495494$; compare with Simpson 1/3 to same precision.

**How solved:** `gaussian_quadrature(f,-1,1,4,'legendre')` vs `simpson(f,-1,1,8)`. Output (`output2.txt`):

```
Gaussian Quadrature (N=4) = 0.4816354816354816
Simpson's Rule (N=8)      = 0.4881357142840971  (9 evaluations)
```

**Takeaway:** Simpson($N=8$, 9 evals) lands closer to 0.48749… than 4-pt Gauss here — rational (non-polynomial) integrand doesn't enjoy Gauss's $2n-1$ exactness; still, Gauss uses only 4 evals.

#### C.10.3 Q3 — 5-pt Laguerre: $\int_0^{\infty}e^{-x}/(1+x)\,dx$

**Question:** 5-pt Gaussian quadrature for $\int_0^{\infty}\frac{e^{-x}}{1+x}dx\approx 0.59635$.

**How solved:** Laguerre with $f(t)=1/(1+t)$ (strip $e^{-t}$ — it's the weight function):

```python
gaussian_quadrature(lambda t: 1/(1+t), 0, 1000, 5, method='laguerre')
```

**Output:** `0.5950840879689522` (error ≈ $1.3\times10^{-3}$ vs 0.59635).

#### C.10.4 Q4 — Laguerre exactness check: $\int_0^{\infty}e^{-x}x^4dx=4!=24$

**Question:** $n$-pt Gauss–Laguerre for $\int_0^{\infty}e^{-x}x^4dx$; compare with $4!=24$.

**How solved:** integrand minus weight is degree-4 polynomial → $N=3$ gives exactness through degree $2N-1=5\supset 4$:

```python
gaussian_quadrature(lambda t: t**4, 0, 1000, 3, method='laguerre')
```

**Output (`output4.txt`):**

```
Gaussian Laguerre Quadrature (N=3) = 24.0
Analytical value (4!) = 24
Difference = 0.0, Relative error = 0.0
```

**Takeaway (best viva demo of §6.8):** Gauss–Laguerre with 3 points is *exact* on $x^4$ — the $2n-1$ theorem in action.

---

---

### C.11 Mid-Semester Examination — `Midsem/`

#### C.11.1 Q1 — Microscopic Particle Diffusion (Ehrenfest Model)
- **Problem:** Simulate $N = 5000$ particles initially in the left compartment ($N_L = N, N_R = 0$) diffusing through a small hole over 30,000 steps.
- **Method:** LCG pseudo-random choice of particle $i \in [0, N-1]$ uniformly. State array toggles compartment.
- **Result:** $N_L(t)$ exponentially relaxes from 5000 towards macroscopic equilibrium ($N/2 = 2500$) with statistical fluctuations $\sim \sqrt{N} \approx 50$. Saved in `equilibrium.png`.

#### C.11.2 Q2 — Matrix Inversion via Doolittle LU Decomposition
- **Problem:** Find the inverse of a $5 \times 5$ matrix `midsem_mat1` using LU decomposition, verifying non-zero determinant.
- **Method:** $A = L U$ factorization, determinant $\det(A) = \prod U_{ii} = 68.717 \ne 0$. Column-by-column forward/backward solve $A \mathbf{x}_j = \mathbf{e}_j$.
- **Result:** Computed $A^{-1}$ rounded to 3 decimal places matching `output2.txt`.

#### C.11.3 Q3 — Wien's Displacement Law via Newton-Raphson
- **Problem:** Solve transcendental equation $(x - 5)e^x + 5 = 0$ for $x = hc / (\lambda_{\max} k T)$ to compute Wien's constant $b$.
- **Method:** 1D Newton-Raphson with $f(x) = (x - 5)e^x + 5$ and $f'(x) = (x - 4)e^x$.
- **Result:** Root converges to $x = 4.965114$. Wien's constant $b = hc / (x k) = 2.8978 \times 10^{-3} \text{ m}\cdot\text{K}$ matching CODATA value.

#### C.11.4 Q4 — Non-Uniform Beam Center of Mass via Midpoint Rule
- **Problem:** 2 m beam with density $\lambda(x) = x^2$. Find center of mass $x_{\text{cm}} = \int x \lambda dx / \int \lambda dx$ to 4 decimal places.
- **Method:** Composite midpoint rule with $N$ determined from error bound formulas ($N_1 = 116, N_2 = 200$).
- **Result:** Mass $M = 2.6667$, Moment $= 4.0000$, $x_{\text{cm}} = 1.5000\text{ m}$ (exact analytical value is $3L/4 = 1.5\text{ m}$).

#### C.11.5 Q5 — Ground State Hydrogen $\langle r \rangle$ via Gauss-Legendre Quadrature
- **Problem:** Compute expectation value $\langle r \rangle = \int_0^\infty 4 u^3 e^{-2u} du$ using radial wavefunction $R_{10}(r)$.
- **Method:** 6-point Gauss-Legendre quadrature with upper cutoff $b = 6.0 a_0$ where probability density decays below $10^{-4}$.
- **Result:** $\langle r \rangle = 1.5000 a_0$ (exact quantum mechanical value is $\frac{3}{2} a_0$).

---

### C.12 Differential Equations I (Euler & Predictor-Corrector) — `assign7a_fwbk/`

- **Problem:** Solve $dy/dx = y - x^2$ ($y(0)=0$) and $dy/dx = (x+y)^2$ ($y(0)=1$) using Forward Euler, Backward Euler, and Predictor-Corrector across step sizes $h \in \{0.02, 0.1, 0.3, 0.5\}$.
- **Method:**
  - Forward Euler: $y_{i+1} = y_i + h f(x_i, y_i)$
  - Backward Euler: Implicit step solved via Newton-Raphson at each increment.
  - Predictor-Corrector: Heun's slope averaging.
- **Findings:** Forward Euler oscillates and drifts significantly at $h = 0.5$; Backward Euler remains stable but introduces damping; Predictor-Corrector maintains high fidelity across all tested step sizes.

---

### C.13 Differential Equations II (Runge-Kutta 4th Order) — `assign7b_rk4/`

- **Q1 (Library):** Unified `rk4` implementation supporting both scalar ODEs and coupled first-order vector systems.
- **Q2:** Solution of $dy/dx = (x+y)^2$ compared against analytical $\tan(x + \pi/4) - x$ on $[0, \pi/5]$ with $h = 0.05, 0.15$.
- **Q3 (Damped SHO):** $\ddot{x} + \mu \dot{x} + \omega^2 x = 0$ ($k=1, m=1, \mu=0.15$) integrated over $t \in [0, 40]$. Generated phase portraits and energy decay curve $E(t)$.
- **Q4 (Chaotic Dynamics):** 3D integration of the Rössler attractor ($a=0.2, b=0.2, c=5.7$), producing the characteristic folded single-band strange attractor plot (`q4_rossler.png`).

---

### C.14 Boundary Value Problems & FEM — `assign7c_sht/`

- **Q1 (Beam Deflection):** Cantilever beam under load $q = 5000 \text{ N/m}$, Young's modulus $Y = 2 \times 10^{11} \text{ N/m}^2$. Solved via Shooting Method with RK4 + Lagrange interpolation to match $y(L) = 0$. Computed maximum depression $y_{\max} = -1.2500 \times 10^{-3} \text{ m}$ at midpoint $x = 1.0\text{ m}$.
- **Q2 (Heat Conduction):** Uninsulated rod with ambient heat dissipation $T'' + \alpha(T_a - T) = 0$, $T(0)=40^\circ\text{C}, T(10)=200^\circ\text{C}$. Shooting method determines location where $T = 100^\circ\text{C}$ ($x = 6.32\text{ m}$).
- **Q3 (1D FEM):** Steady-state Poisson $y''(x) = -2$, $y(0)=1, y(1)=3$ solved via Finite Element Method using piecewise linear hat functions and Cholesky decomposition. Matches analytical $y(x) = -x^2 + 3x + 1$ exactly at grid nodes.

---

### C.15 Partial Differential Equations — `assign7d_felem/`

- **Q1 (Explicit FTCS):** 1D Heat equation $u_t = u_{xx}$ with $u(x, 0) = 20|\sin(\pi x)|$ on a metal rod of length 2. Position grid $N_x = 20$, time grid $N_t = 2000$. Verified stability condition $\alpha = h_t / h_x^2 = 0.2 \le 0.5$. Visualized thermal decay profiles across time steps.
- **Q2 (Crank-Nicolson):** Convection-diffusion-reaction PDE $u_t = a u_{xx} + b u_x + c u$ ($a=2, b=-1, c=1$) with initial profile $u(x, 0) = 2 + x^2$ and boundary conditions $u(0, t)=2, u(1, t)=3$. Formulated matrix system $A \mathbf{u}^{j+1} = B \mathbf{u}^j + \mathbf{b}_{\text{bc}}$ and solved over 5000 time steps. Saved in `crank_nicolson_lines.png` and contour map.

---

### C.16 Eigenvalue Problems — `assign8_eval/`

- **Q1 (Library):** Power Iteration with Rayleigh quotient and automatic normalization at each iteration.
- **Q2:** Applied to $4 \times 4$ matrix `asgn8_mat1`:
  $$\text{Dominant eigenvalue } \lambda_1 = 8.000000, \quad \mathbf{v}_1 = [0.0, 0.0, 0.707107, 0.707107]^T$$
- **Q3:** 5x5 tridiagonal Toeplitz matrix `asgn8_mat2`:
  - First eigenvalue $\lambda_1 = 3.732051$, eigenvector $\mathbf{v}_1 = [0.288675, -0.5, 0.57735, -0.5, 0.288675]^T$.
  - Matrix deflation $A' = A - \lambda_1 \mathbf{u}_1 \mathbf{u}_1^T$ isolates second eigenvalue $\lambda_2 = 3.000000$, eigenvector $\mathbf{v}_2 = [-0.5, 0.5, 0.0, -0.5, 0.5]^T$.
  - Perfect agreement with analytical formula $\lambda_k = 2 + 2\cos(k\pi/6)$.

---

### C.17 Interpolation & Model Comparison — `assign9a_fit/`

- **Q1 (Library):** Lagrange interpolation returning both estimated value $y(x)$ and ascending polynomial coefficients $[c_0, c_1, \dots, c_n]$.
- **Q2 (Library):** Weighted linear regression returning slope, intercept, and Pearson correlation coefficient $r^2$.
- **Q3:** Lagrange interpolation on 5 data points estimating $y(6.7) = 33.500000$ with explicit coefficients printed.
- **Q4 (Model Race):** Experimental data fitted against Power Law ($y = a x^b$) vs. Exponential ($y = a e^{-bx}$):
  - Power Law: $a = 21.78, b = -0.53, R^2 = 0.9926$
  - Exponential: $a = 12.01, b = 0.05, R^2 = 0.9234$
  - Conclusion: Power Law provides significantly superior fit ($R^2 = 0.9926 > 0.9234$).

---

### C.18 Polynomial & Chebyshev Regression — `assign9b_fit/`

- **Q1 (Library):** Polynomial least squares of arbitrary order $k$ via normal equations, returning coefficients, Hankel matrix, and Frobenius condition number.
- **Q2 (Library):** Shifted Chebyshev basis regression on $[0, 1]$ using $\phi_m(x) = T_m(2x - 1)$.
- **Q3:** Cubic polynomial least squares on experimental dataset `assign9fit.txt`. Hankel matrix condition number $= 3.32 \times 10^4$.
- **Q4 (Stability Benchmark):** Comparative study:
  - Standard Monomial Hankel condition number: $\kappa = 33241.6$
  - Shifted Chebyshev matrix condition number: $\kappa = 14.8$
  - *Takeaway:* Orthogonal Chebyshev basis reduces ill-conditioning by over 3 orders of magnitude while yielding identical analytical curves.

---

### C.19 Complete Master Map: Repo Folders vs. Slide Sections

| Repo Folder | Topic / Routines | Slide Deck Reference |
|---|---|---|
| `Assgn0` | Loops, series sums, file I/O, `MyComplex` | `intro.pdf` |
| `assign1_rand` | LCG, $\pi$, radioactive decay, exponential sampling | `random.pdf` |
| `assign2_gj_lu` | Gauss-Jordan pivoting, inversion, LU decomposition | `gaussj.pdf`, `ludecomp.pdf` |
| `assign3_lu_chsk` | Forward/backward substitution, $6 \times 6$ LU, Cholesky | `ludecomp.pdf` |
| `assign4_jac_gs` | Jacobi, Gauss-Seidel, SOR relaxation ($\omega = 1.57$) | `itermat.pdf` |
| `assign5_bi_rf` | Root bracketing, Bisection, Regula Falsi | `root.pdf` |
| `assign5_fx_nr` | Fixed-point, 1D & Multivariate Newton-Raphson | `root.pdf` |
| `assign5_poly` | Horner's scheme, Laguerre roots, synthetic deflation | `root.pdf` |
| `assign6_midtrap` | Composite Midpoint and Trapezoidal quadrature | `integ.pdf` |
| `assign6b_simmc` | Simpson 1/3-rule, Monte Carlo variance $O(N^{-1/2})$ | `integ.pdf` |
| `assign6c_gq` | Gauss-Legendre & Gauss-Laguerre quadrature ($N \le 6$) | `integ.pdf` |
| `Midsem` | Comprehensive mid-semester examination suite | All slides 1–7 |
| `assign7a_fwbk` | Forward/Backward Euler, Heun predictor-corrector | `diffeq.pdf` (Part I) |
| `assign7b_rk4` | Runge-Kutta 4th Order, damped SHO, Rössler attractor | `diffeq.pdf` (Part I) |
| `assign7c_sht` | Shooting method (beam, rod), 1D Finite Element Method | `diffeq.pdf` (Part I) |
| `assign7d_felem` | Explicit FTCS heat equation, Crank-Nicolson PDE | `diffeq.pdf` (Part II) |
| `assign8_eval` | Power iteration, Rayleigh quotient, matrix deflation | `eigv.pdf` |
| `assign9a_fit` | Lagrange interpolation, weighted linear regression | `datafit.pdf` |
| `assign9b_fit` | Polynomial fit, Hankel conditioning, shifted Chebyshev | `datafit.pdf` |

---

*End of master notes — covers all 10 slide decks + full assignment repo walkthrough.*
