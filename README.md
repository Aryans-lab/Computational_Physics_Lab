# PHY341/745 — Physics Computer Lab
### National Institute of Science Education and Research (NISER), Bhubaneswar

**Student:** Aryan Bandyopadhyay &nbsp;|&nbsp; **Roll No.:** 2411014  
**Course:** PHY341 (Int-MSc Sem-V) / PHY745 (Int-MSc-PhD Sem-III) · 4 Credits  
**Instructor:** Prof. Subhasish Basak (`sbasak@niser.ac.in`)  
**TAs:** Hemant Lohumi · Chinmoy Samanta  

---

## Overview

This repository contains all programming assignments, helper libraries, and data files for the **Physics Computer Lab** course at NISER. The course develops practical skills in scientific computing using **Python** (and optionally C/C++) to implement numerical methods from scratch — without relying on high-level library routines (NumPy, SciPy) for the core algorithms.

> **Key constraint from the course:** All numerical algorithms (matrix operations, root-finding, integration, ODEs, etc.) must be implemented manually in a custom library (`mylib.py`). System-built routines for major algorithms are **not** permitted.

---

## Repository Structure

```
Computational_Physics_Lab/
│
├── README.md                          ← This file
├── LICENSE
│
├── helper/                            ← Shared utility library (grows every week)
│   ├── mylib.py                       ← All numerical routines / algorithms
│   ├── main.py                        ← Generic front code (I/O + calls to mylib)
│   └── read_write.py                  ← Minimal file I/O snippet reference
│
└── assignments/
    └── week00_warmup/                 ← Warmup assignment (Python basics)
        ├── README.md                  ← Problem statements & notes
        ├── question1.py               ← Sum of 20 even numbers + factorial
        ├── question2.py               ← GP and HP series sums
        ├── question3.py               ← Matrix/vector operations from files
        ├── question5.py               ← Custom complex number class
        ├── data/
        │   ├── asgn0_matA             ← 3×3 matrix A (ASCII)
        │   ├── asgn0_matB             ← 3×3 matrix B (ASCII)
        │   ├── asgn0_vecC             ← 3×1 column vector C (ASCII)
        │   └── asgn0_vecD             ← 3×1 column vector D (ASCII)
        └── output/
            └── question3_output.txt   ← Verified output for Q3
```

---

## Course Syllabus

| # | Topic | Status |
|---|-------|--------|
| 0 | **Warmup** — Python basics, loops, file I/O, classes, matrix ops | ✅ Done |
| 1 | Random number generation | 🔜 |
| 2 | Gauss-Jordan elimination (linear systems, inverse, determinant) | 🔜 |
| 3 | LU decomposition | 🔜 |
| 4 | Jacobi / Gauss-Seidel iterative methods | 🔜 |
| 5 | Root finding — Regula Falsi, Newton-Raphson, Laguerre | 🔜 |
| 6 | Numerical integration — Midpoint, Trapezoidal, Simpson, Monte Carlo | 🔜 |
| 7 | ODEs — Euler, Predictor-Corrector, Runge-Kutta (IVP & BVP) | 🔜 |
| 8 | Least-square fitting | 🔜 |
| — | *If time permits:* PDEs & Eigenvalue problems | — |
| — | *If time permits:* ML regression & classification | — |

---

## Programming Conventions

All code in this repository strictly follows the course guidelines:

1. **Header comment** — Every file begins with the problem description, author name, and roll number.
2. **No hardwired input** — All input is read from external ASCII files via `sys.argv` or relative paths.
3. **Separate output files** — Results are written to output files; output is never just printed to console (except in warmup).
4. **Library / Front-code separation** — All numerical algorithms live in `mylib.py`. The `main.py` front code only handles I/O, file opening, and calls to library functions.
5. **No NumPy/SciPy for algorithms** — Core routines are implemented from scratch using pure Python.
6. **Loop-based initialization** — Matrices and series coefficients are always initialized inside loops (never manually assigned).
7. **Meaningful comments** — Each computing step is explained without being verbose.

---

## Warmup Assignment (Week 00)

**Problem set covering Python fundamentals required for the course.**

### Q1 — Loops & Arithmetic
Compute the **sum of the first 20 even numbers** and the **factorial of 8** using explicit loops (no `math.factorial`).

```
Sum of first 20 even numbers: 420
Factorial of 8: 40320
```

### Q2 — Series Sums (GP & HP)
Sum **N = 15 terms** of a Geometric Progression (r = 0.5, t₀ = 1.25) and a Harmonic Progression (d = 1.5, t₀ = 1.25) without using closed-form formulae.

```
Sum of GP series: 2.4999237060546875
Sum of HP series: 2.4139570733659186
```

### Q3 — Matrix & Vector Operations from Files
Read matrices **A** (3×3), **B** (3×3), column vectors **C** (3×1) and **D** (3×1) from ASCII data files. Compute:
- **AB** — matrix product
- **BC** — matrix-vector product
- **D·C** — dot product

```
Matrix AB:
  [-0.300000,  -3.500000,   5.200000]
  [-4.500000,  -2.000000,   4.500000]
  [ 9.300000,   0.800000,  -7.000000]

Matrix BC:
  [ 1.000000]
  [-5.750000]
  [-9.000000]

Dot product D·C: -3.5
```

> **Note:** The (0,0) entry of AB shows a small floating-point rounding artefact (~10⁻¹⁶), which is expected for IEEE 754 double precision.

### Q5 — Custom Complex Number Class
Implement a `MyComplex` class supporting **addition, subtraction, multiplication, and modulus** for `(1.3 − 2.2j)` and `(−0.8 + 1.7j)`.

```
Sum:        0.5 - 0.5j
Difference: 2.1 - 3.9j
Product:    2.7 + 3.97j
|c1|: 2.5554
|c2|: 1.8788
```

---

## Helper Library (`mylib.py`)

The shared library currently implements:

| Function | Description |
|----------|-------------|
| `read_tokens(filename)` | Read all whitespace-separated tokens from a file |
| `parse_numbers(tokens)` | Parse token list into floats (auto-detects leading N) |
| `read_second_column(filename)` | Parse two-column `index value` files |
| `sum_numbers(values)` | Summation loop |
| `read_matrix_add(tokens)` | Parse two n×m matrices for addition |
| `matrix_add(A, B, n, m)` | Element-wise matrix addition |
| `read_matrix_multiply(tokens)` | Parse two matrices for multiplication |
| `matrix_multiply(A, B, n, p, m)` | Standard O(n³) matrix multiplication |
| `matrix_to_string(C)` | Format a matrix as a printable string |

This library will be extended with new algorithms each week.

---

## How to Run

### Warmup scripts
```bash
# Question 1
python assignments/week00_warmup/question1.py

# Question 2
python assignments/week00_warmup/question2.py

# Question 3 (reads data files automatically from the data/ directory)
python assignments/week00_warmup/question3.py

# Question 5
python assignments/week00_warmup/question5.py
```

### Generic front code (helper/)
```bash
# Sum numbers in a file
python helper/main.py sum input.txt output.txt

# Matrix addition
python helper/main.py matadd input.txt output.txt

# Matrix multiplication
python helper/main.py matmul input.txt output.txt
```

---

## Marks Distribution

| Component | Marks |
|-----------|-------|
| Assignments | 35 |
| Mid-semester | 15 |
| DIY project | 20 |
| End-semester | 30 |
| **Total** | **100** |

---

## Academic Integrity

All code in this repository is original work by **Aryan Bandyopadhyay (Roll No. 2411014)**. As per course policy, suspected copying results in zero marks for all parties involved.

---

## License

This repository is released under the [MIT License](LICENSE) for academic reference purposes.
