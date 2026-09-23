# Week 01 — Pseudo-Random Numbers & Stochastic Simulations

Randomness is the tool that makes the next three weeks possible. This week
builds it from nothing: a full-period Linear Congruential Generator,
distribution sampling by inverse transform, and two physical applications —
a decay chain and an exponential distribution.

**Convention (course rule 8):** the core numerics use only the standard
library; Matplotlib is used solely for the figures.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Lag-correlation study: chaos map $x \mapsto 1 - \|c x (1 - x)\|$ vs the LCG | the LCG's $x_i$–$x_{i+1}$ points stay uncorrelated; the logistic map at $c = 1$ collapses onto a curve |
| 2 | `question2.py` | Full-period audit of the LCG (Hull-Dobell: $m = 2^{15}$, $a \equiv 1 \bmod 4$, $c$ odd) | all $32\,768$ states appear before any repeat |
| 3 | `question3.py` | Monte-Carlo $\pi$ with nested reseeding per $N = 20 \dots 5000$ | the residual plot tracks $O(N^{-1/2})$ |
| 4 | `question4.py` | Radioactive decay chain $A \to B \to C$, each nucleus deciding by $P = \lambda \Delta t$ | the classic rise-and-fall of $N_B$ emerges from pure counting |
| 5 | `question5.py` | Exponential deviates by inverse transform, $y = -\lambda^{-1}\ln u$ | sample mean of $50\,000$ draws within a couple of $\sigma$ of $1/\lambda$ |

## Running

Scripts are self-contained: they read from `data/`, write to
`output/qN_output.txt`, and save figures to `figures/` — all relative to the
folder, so they run with no arguments:

```bash
python question3.py          # or from anywhere: python <path>/question3.py
```

`mylib.py` in this folder is a verbatim bundle of the canonical library at
the repository root (`make bundle` regenerates it), which is what makes the
folder independently submittable.
