# Week 06b — Simpson's 1/3 Rule & Monte Carlo Integration

Two philosophies of integration meet here. Simpson's rule is deterministic
and $O(h^4)$ — four digits of accuracy for four times the evaluations.
Monte Carlo is random and only $O(N^{-1/2})$ — *but* that rate does not
depend on the dimension, which is what makes it the only option for the
high-dimensional integrals of statistical physics.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Validation: Simpson's exactness on a cubic (its order of precision is 3), Monte-Carlo statistical check (estimate within a few $\sigma_I$ of exact) | Simpson on the cubic: $\sim 10^{-15}$; MC deviation $\sim O(1)\,\sigma$ |
| 2 | `question2.py` | $\int_1^2 dx/x$ and $\int_0^{\pi/2} x\cos x\,dx$ to 6 d.p., with $N$ **chosen from the error bounds**, not by trial | Simpson needs $N = 20$ and $22$ subintervals where the midpoint rule needs $289$ and $610$ |
| 3 | `question3.py` | Monte-Carlo $\int_{-1}^{1} \sin^2 x\,dx$ (exact $= 1 - \tfrac12\sin 2$), $N = 1100 \dots 50000$, plots of value and deviation vs $N$ | first $N$ with $|err| < 10^{-4}$ reported; deviation wanders in the $1/\sqrt{N}$ band |

## What's worth reading

- **The $N$ from the bound is the whole technique.**
  $N_{\text{mid}} \sim \epsilon^{-1/2}$ but $N_{\text{simp}} \sim
  \epsilon^{-1/4}$: for the same 6 d.p. target the gap is 14–28×. The
  driver prints the derivative bounds it used, so the sizing is auditable.
- **Reading the Monte-Carlo plots.** The estimate does *not* converge
  smoothly — it fluctuates around the truth with size
  $\sigma_I = (b-a)\sigma_f/\sqrt{N}$. That is not a bug in the method; it
  is the central limit theorem, visible.
- **Seeded, therefore reproducible.** The MC routine takes a seed; the
  golden-output test suite re-runs this driver and byte-compares the result,
  so "random" here still means exactly-reproducible.

## Running

```bash
python question3.py   # zero arguments: writes output/q3_output.txt + figures/q3_*.png
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
