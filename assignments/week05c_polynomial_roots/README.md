# Week 05c — Polynomial Roots: Laguerre's Method & Deflation

Newtons for polynomials: given only the coefficient list, find *all* the
real roots. Laguerre's method gives cubic convergence near a simple root,
and synthetic (Horner) deflation strips each root out, reducing the degree
by one for the next pass.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Library validation on $(x-1)(x-2)(x-3)$: every root checked by substitution, deflation chain traced | all three roots to $\sim 10^{-8}$ |
| 2 | `question2.py` | The lab's three polynomials — $x^4 - x^3 - 7x^2 + x + 6$; $x^4 - 5x^2 + 4$; $2x^5 - 19.5x^3 + 0.5x^2 + 13.5x - 4.5$ | P1, P2 exact to $\sim 10^{-10}$; P3's **double root** at $0.5$ splits into $0.49998 / 0.50002$ |

## What's worth reading

- **The double-root story is the point.** Multiplicity $m$ degrades
  Laguerre from cubic to *linear* convergence, so a double root found by
  deflation lands as two nearby values, not one exact one. The driver
  reports the split honestly and the driver's assertion absorbs it
  (tolerance $10^{-4}$) — a naive "must equal 0.5 exactly" test would
  fail a correct implementation.
- **Real roots only, by contract.** The library raises `ValueError` when
  the next Laguerre step has no real solution (the radicand is negative and
  stays negative) — a polynomial with complex roots is out of scope, and
  the failure is loud rather than a hallucinated real root.
- **Horner everywhere.** Evaluation and deflation both run through
  `polynomial_value` / `synthetic_division` in $O(n)$; the unit tests check
  both against `numpy.polyval` and known quotients.

## Running

```bash
python question2.py   # zero arguments: writes output/q2_output.txt
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
