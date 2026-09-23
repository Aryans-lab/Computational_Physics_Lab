# Week 03 — LU Forward-Backward & Cholesky

The payoff of factorisation: once $A$ is written as $LU$ (or $LL^T$ for
symmetric positive-definite $A$), a solve costs $O(n^2)$ instead of
$O(n^3)$, and the factor exposes the conditioning of the system.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Cholesky routine with an explicit symmetry check, verified on the SPD matrix `data/asgn3_mat2`; both failure modes (non-symmetric, non-positive-definite) are demonstrated to raise | $LL^T$ reproduces $A$ to $10^{-15}$; guards fire correctly |
| 2 | `question2.py` | The 6×6 system `data/asgn3_mat1`, `data/asgn3_vec1` solved by LU + forward-backward substitution | residual $< 10^{-12}$ |
| 3 | `question3.py` | Cholesky on `data/asgn3_mat2` with $b = (3,3,1,3)$; solve via $Ly = b$, $L^T x = y$ | $x = (0, 1, 1, 1)$ — the exact integer solution, recovered to machine precision |

## What's worth reading

- **Cholesky is LU that respects symmetry.** It costs about half the
  flops of a full LU ($\tfrac{1}{3}n^3$ vs $\tfrac{2}{3}n^3$) because only
  the triangle below the diagonal is ever touched — the driver's
  round-trip check $LL^T \approx A$ is the whole verification.
- **The symmetry check earns its keep.** A matrix that is SPD *up to a
  transcription error* fails the symmetry test before any square root is
  taken, turning a silent garbage answer into a clear `ValueError`.

## Running

```bash
python question1.py   # zero arguments: reads data/, writes output/q1_output.txt
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
