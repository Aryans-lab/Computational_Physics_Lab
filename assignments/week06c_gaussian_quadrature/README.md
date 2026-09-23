# Week 06c — Gaussian Quadrature: Legendre & Laguerre

The end of the quadrature line: choose the *nodes* as free parameters, and
$n$ of them buy you exactness for every polynomial of degree $2n - 1$.
That is the best possible order for $n$ nodes (Riemann's theorem), and it
explains why 4 Gauss points beat 9 Simpson points.

Two families:

- **Gauss-Legendre** — finite interval $[a, b]$, mapped from the canonical
  $[-1, 1]$ by $x = \frac{b-a}{2}t + \frac{b+a}{2}$.
- **Gauss-Laguerre** — the semi-infinite interval $[0, \infty)$ with weight
  $e^{-x}$, which makes it the natural tool for the radial and
  Boltzmann-factor integrals of physics.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | The one routine, two families: prints the $n = 1 \dots 6$ node/weight tables and *proves* the order of precision ($\int x^{2n-1}$ exact in both families) | every entry at round-off ($\le 10^{-12}$) |
| 2 | `question2.py` | $\int_{-1}^{1} \frac{x^2}{1+x^4}\,dx$ by 4-point Gauss-Legendre, compared with Simpson | GQ: $0.4816354816$ vs true $0.4874954944$ — 4 evaluations where Simpson needs 9 |
| 3 | `question3.py` | $\int_0^\infty \frac{e^{-x}}{1+x}\,dx$ by 5-point Gauss-Laguerre (the weight *is* $e^{-x}$) | $0.595084\ldots$ vs $e\,E_1(1) = 0.596347\ldots$ — and no cutoff to choose |
| 4 | `question4.py` | $\int_0^\infty e^{-x} x^4\,dx$: a degree-4 polynomial, so the 3-node rule (exact to degree 5) must nail $4!$ exactly | $24.000000000000\ldots$, error at round-off |

## What's worth reading

- **The moment identities are the smoke test.** $\int_0^\infty x^k e^{-x}
  dx = k!$ and $\int_{-1}^1 x^k dx = \delta_{k\,\text{even}}\, 2/(k+1)$:
  any wrong node or weight in the table shows up immediately in problem 1's
  table.
- **One routine, two families.** `gaussian_quadrature(f, a, b, N, method)`
  dispatches on `method`; the Laguerre branch ignores $a, b$ by convention,
  which is the same interface shape a weighted-orthogonal-polynomial
  generalisation (Chebyshev, Hermite) would take.
- **Why problem 3 needs no truncation.** A naive quadrature of
  $\int_0^\infty$ has to pick a cutoff and check convergence; Laguerre's
  nodes already sample the tail according to $e^{-x}$.

## Running

```bash
python question1.py   # zero arguments: writes output/q1_output.txt
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
