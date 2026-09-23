# Week 05b — Fixed-Point Iteration & Newton's Method

Two ways to convert $f(x) = 0$ into a fixed point $x = g(x)$: the
*designer's* choice (fixed-point) and the *geometric* choice (Newton).
Newton's price is a derivative; its reward is quadratic convergence — the
digit count *doubles* every step.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Library validation: fixed point on $x^2 - 2x - 3$; Newton 1-D on $x^2 - 2$ (quadratic-convergence check); multivariate Newton with a *numerical* Jacobian | roots at $3$, $\sqrt2$, and $(3, 4)$ with $\max\|F\| < 10^{-10}$ |
| 2 | `question2.py` | $f(x) = 3x + \sin x - e^x$ on $[-1.5, 1.5]$: bisection vs regula falsi vs Newton (analytic *and* numeric derivative), accuracy $10^{-6}$ | Newton reaches it in a handful of iterations where bisection needs ~22 |
| 3 | `question3.py` | Fixed-point for *both* roots of $x^2 - 2x - 3$: $x = \sqrt{2x+3}$ near $3$, $x = 3/(x-2)$ near $-1$ | both converge; the choice of rewrite $g$ is the whole art |
| 4 | `question4.py` | The coupled system $3x^2 - 2xy - 3 = 0$, $3y^2 - 4xy = 0$ from $(2.5, 3.1)$ | solution $(3, 4)$ in ~4 iterations |

## What's worth reading

- **The Newton step does not invert $J$.** The multivariate routine solves
  $J\,\delta x = -F$ with Gauss-Jordan — the same step as
  $-J^{-1}F$, without ever forming the inverse (which squares the
  condition number).
- **Convergence rates, measured.** Problem 2 puts all four methods next to
  each other so the table *is* the lesson: linear → superlinear → quadratic.
- **The fixed-point trap.** $g(x) = x^2 - 2x - 3 + x$ "works" algebraically
  but has $|g'| > 1$ at both roots and diverges from them; the two
  converging rewrites in problem 3 each satisfy $|g'| < 1$ at their target
  ($\tfrac13$ and $\tfrac19$).

## Running

```bash
python question4.py   # zero arguments: writes output/q4_output.txt
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
