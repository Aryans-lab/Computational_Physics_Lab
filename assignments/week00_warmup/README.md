# Week 00 — Warmup Assignment

## Problem Set

This warmup covers Python fundamentals needed throughout the course.
All solutions follow the course programming conventions (file I/O, library
separation, no NumPy for core algorithms).

---

### Q1 – Loops & Arithmetic (`question1.py`)

**Problem:** Compute the sum of the first 20 even numbers and the factorial
of 8 using an explicit `for` loop.  
**Restriction:** No `math.factorial` or `sum()` for the computation itself.

| Quantity | Result |
|---|---|
| Sum of first 20 even numbers | 420 |
| 8! | 40 320 |

---

### Q2 – Series Sums: GP & HP (`question2.py`)

**Problem:** Sum N = 15 terms of:
- A **Geometric Progression** with first term t₀ = 1.25 and ratio r = 0.5
- A **Harmonic Progression** with underlying AP: first term t₀ = 1.25, common difference d = 1.5

**Restriction:** Closed-form formulae are **not** allowed; terms must be accumulated iteratively.

| Series | Sum |
|---|---|
| GP (15 terms) | 2.4999237… |
| HP (15 terms) | 2.4139570… |

> **Analysis:** The GP sum converges to t₀/(1−r) = 1.25/0.5 = **2.5** as N→∞.
> With 15 terms, the truncation error is approximately 7.6 × 10⁻⁵.

---

### Q3 – Matrix & Vector Operations from Files (`question3.py`)

**Problem:** Read matrices **A** (3×3), **B** (3×3) and column vectors **C**, **D** (3×1)
from ASCII data files in `data/`. Compute **AB**, **BC**, and **D·C**.

**Input data:**
```
A =                     B =
  2.0  -3.0   1.4         0.0  -1.0   1.0
  2.5   1.0  -2.0         1.5   0.5  -2.0
 -0.8   0.0   3.1         3.0   0.0  -2.0

C = [-2.0,  0.5,  1.5]^T
D = [ 1.0,  0.0, -1.0]^T
```

**Output** (see `output/question3_output.txt` for the exact file):
```
Matrix AB:
  [ -0.300000,  -3.500000,   5.200000]
  [ -4.500000,  -2.000000,   4.500000]
  [  9.300000,   0.800000,  -7.000000]

Matrix BC:
  [  1.000000]
  [ -5.750000]
  [ -9.000000]

Dot product D.C = -3.5
```

> **Floating-point note:** The (0,0) entry of **AB** has a residual ~10⁻¹⁶ error.
> This is a standard IEEE 754 artefact, not a coding bug.

---

### Q4 — (Not assigned this week)

---

### Q5 – Custom Complex Number Class (`question5.py`)

**Problem:** Implement `MyComplex` class supporting addition, subtraction,
multiplication, and modulus for c₁ = (1.3 − 2.2j) and c₂ = (−0.8 + 1.7j).

| Operation | Result |
|---|---|
| c₁ + c₂ | 0.5 − 0.5j |
| c₁ − c₂ | 2.1 − 3.9j |
| c₁ × c₂ | 2.7 + 3.97j |
| \|c₁\| | 2.5554 |
| \|c₂\| | 1.8788 |

> **Design note:** The class uses Python dunder methods (`__add__`, `__sub__`,
> `__mul__`) for natural operator syntax (`c1 + c2`), which is cleaner
> than the original `c1.add_complex(c1, c2)` pattern.

---

## Files in This Directory

```
week00_warmup/
├── README.md
├── question1.py
├── question2.py
├── question3.py
├── question5.py
├── data/
│   ├── asgn0_matA
│   ├── asgn0_matB
│   ├── asgn0_vecC
│   └── asgn0_vecD
└── output/
    └── question3_output.txt
```

## Running the Scripts

```bash
# from the repo root:
python assignments/week00_warmup/question1.py
python assignments/week00_warmup/question2.py
python assignments/week00_warmup/question3.py   # reads from data/ automatically
python assignments/week00_warmup/question5.py
```
