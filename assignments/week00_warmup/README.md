# Week 00 — Warmup Assignment

## Problem Set

This warmup covers Python fundamentals required throughout the course.
All solutions strictly follow the course programming conventions:
- Parameters/inputs read from **external files** via `sys.argv` (Rules 2 & 6)
- All output written to **separate output files** (Rule 3)
- Functions isolated in their own scripts (warmup only; from Week 01 onwards all routines go into `mylib.py`)
- No NumPy/SciPy for core computation (Rule 8)

---

### Q1 – Loops & Arithmetic (`question1.py`)

**Problem:** Compute the sum of the first N even numbers and the factorial of M,
where N and M are read from `data/q1_input.txt`.

**Input** (`data/q1_input.txt`): `20  8`

| Quantity | Result |
|---|---|
| Sum of first 20 even numbers | 420 |
| 8! | 40 320 |

**Output file:** `output/q1_output.txt`

---

### Q2 – Series Sums: GP & HP (`question2.py`)

**Problem:** Sum N terms of a GP and a HP.  
Parameters are read from `data/q2_input.txt`.  
Closed-form formulae are **not** used; terms are accumulated iteratively.

**Input** (`data/q2_input.txt`):  N=15, t₀=1.25, r=0.5, d=1.5

| Series | Sum |
|---|---|
| GP (15 terms) | 2.4999237… |
| HP (15 terms) | 2.4139570… |

> **Analysis:** The GP sum converges to t₀/(1−r) = **2.5** as N→∞.
> With 15 terms the truncation error is ≈77.6 × 10⁻⁶.

**Output file:** `output/q2_output.txt`

---

### Q3 – Matrix & Vector Operations from Files (`question3.py`)

**Problem:** Read matrices A, B and vectors C, D from `data/`. Compute AB, BC, D·C.

**Input data:**
```
A =                     B =
  2.0  -3.0   1.4         0.0  -1.0   1.0
  2.5   1.0  -2.0         1.5   0.5  -2.0
 -0.8   0.0   3.1         3.0   0.0  -2.0

C = [-2.0,  0.5,  1.5]^T
D = [ 1.0,  0.0, -1.0]^T
```

**Output** (`output/q3_output.txt`):
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

> **Floating-point note:** AB[0][0] has a residual ~10⁻¹⁶ error -- standard IEEE 754 artefact.

---

### Q4 — (Not assigned this week)

---

### Q5 – Custom Complex Number Class (`question5.py`)

**Problem:** Implement `MyComplex` class; compute operations for c₁ and c₂ read from file.

**Input** (`data/q5_input.txt`): `1.3  -2.2` and `-0.8  1.7`

| Operation | Result |
|---|---|
| c₁ + c₂ | 0.5 − 0.5j |
| c₁ − c₂ | 2.1 − 3.9j |
| c₁ × c₂ | 2.7 + 3.97j |
| \|c₁\| | 2.5553864678 |
| \|c₂\| | 1.8788294228 |

**Output file:** `output/q5_output.txt`

---

## Directory Layout

```
week00_warmup/
├── README.md
├── question1.py
├── question2.py
├── question3.py
├── question5.py
├── data/
│   ├── asgn0_matA       <- 3x3 matrix A
│   ├── asgn0_matB       <- 3x3 matrix B
│   ├── asgn0_vecC       <- 3x1 column vector C
│   ├── asgn0_vecD       <- 3x1 column vector D
│   ├── q1_input.txt     <- N, M for Q1
│   ├── q2_input.txt     <- N, t0, r, d for Q2
│   └── q5_input.txt     <- c1, c2 for Q5
└── output/
    ├── q1_output.txt
    ├── q2_output.txt
    ├── q3_output.txt
    └── q5_output.txt
```

## Running the Scripts

```bash
# from repo root:
python assignments/week00_warmup/question1.py \
    assignments/week00_warmup/data/q1_input.txt \
    assignments/week00_warmup/output/q1_output.txt

python assignments/week00_warmup/question2.py \
    assignments/week00_warmup/data/q2_input.txt \
    assignments/week00_warmup/output/q2_output.txt

python assignments/week00_warmup/question3.py \
    assignments/week00_warmup/output/q3_output.txt

python assignments/week00_warmup/question5.py \
    assignments/week00_warmup/data/q5_input.txt \
    assignments/week00_warmup/output/q5_output.txt
```
