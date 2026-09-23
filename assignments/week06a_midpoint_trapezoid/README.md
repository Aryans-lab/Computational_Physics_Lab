# Week 06a — Closed Newton-Cotes: Midpoint & Trapezoidal

The two rules every quadrature course starts with: slice $[a, b]$ into $N$
pieces, and on each piece approximate $f$ by a constant (midpoint) or a
line (trapezoid). Both are second-order, $O(h^2)$, and both come with a
clean error bound that makes "how many subintervals do I need?" a one-line
calculation.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Validation on $\int_1^2 \frac{dx}{x} = \ln 2$: convergence-table with measured error ratios | ratios → 0.25 on doubling $N$ — the $O(N^{-2})$ law, measured not assumed |
| 2 | `question2.py` | Three lab integrals — $\int_1^2 dx/x$, $\int_0^{\pi/2} x\cos x\,dx$, $\int_0^1 x\arctan x\,dx$ — at $N = 4, 8, 15, 20$ vs the analytical values | all approaches track their $O(N^{-2})$ envelopes |

## Analytical values used

| Integral | Exact |
|---|---|
| $\displaystyle\int_1^2 \frac{dx}{x}$ | $\ln 2 = 0.6931471806$ |
| $\displaystyle\int_0^{\pi/2} x\cos x\,dx$ | $\frac{\pi}{2} - 1 = 0.5707963268$ |
| $\displaystyle\int_0^1 x\arctan x\,dx$ | $\frac{\pi}{4} - \frac12 = 0.2853981634$ |

## What's worth reading

- **The error bound is a *sizing tool*.**
  $\lVert E_{\text{mid}}\rVert \le \frac{(b-a)^3}{24N^2}\max|f''|$ and
  $\lVert E_{\text{trap}}\rVert \le \frac{(b-a)^3}{12N^2}\max|f''|$ let you
  solve for the minimal $N$ *before* integrating — the library ships
  `integration_error_bound_N` for exactly this (used by week06b).
- **Why the ratio column matters.** A rule claimed to be $O(h^p)$ must show
  its error dividing by $\approx 2^p$ when $N$ doubles. The table is the
  proof.

## Running

```bash
python question1.py   # zero arguments: writes output/q1_output.txt
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
