# Week 02 — Gauss-Jordan Elimination & LU Decomposition

The first direct methods for $Ax = b$. Two philosophies: Gauss-Jordan pushes
the augmented matrix all the way to reduced row-echelon form; LU factorises
once ($PA = LU$) and then each solve is two triangular substitutions.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | GJ on a pivot-stress system $A = \begin{bmatrix}1&1&0\\0&0&1\\0&1&0\end{bmatrix}$, $b = (1,2,3)$ — column 2 has a zero diagonal, so a code without row pivoting dies here | $x = (-2, 3, 2)$ |
| 2 | `question2.py` | GJ on the 3-variable system $2y + 5z = 1$, $3x - y + 2z = -2$, $x - y + 3z = 3$ (augmented matrix in `data/asgn2_mat1`) | $x = (-2, -2, 1)$, residual $= 0$ to round-off |
| 3 | `question3.py` | Doolittle LU with partial pivoting on `data/asgn2_mat2`; the driver reassembles $L\,U$ and compares with $P\,A$ to *prove* the factorisation, then solves by forward-backward substitution | $P = (1,2,0)$ swap chain; $\|LU - PA\|_{\max} = 0$ to round-off |

## What's worth reading

- **Pivoting that actually pivots.** Partial pivoting in `lu_decomposition`
  swaps the rows of $L$ *below* the pivot (only the already-computed
  multipliers, columns $0 \dots i-1$) — swapping whole $L$ rows moves the
  unit diagonal and silently corrupts the factorisation.
- **The determinant comes for free.** $\det A = \prod U_{ii} \times
  \text{sign}(P)$, where the sign is computed from the permutation's cycle
  decomposition — an easy place for an off-by-one to flip ~40% of the signs.
- **Singularity is classified, not just raised.** GJ distinguishes
  inconsistent, dependent, and non-square inputs with specific
  `ValueError`s (all covered by unit tests).

## Running

```bash
python question1.py   # zero arguments: reads data/, writes output/q1_output.txt
python question2.py my_input.txt my_output.txt   # explicit paths still work
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
