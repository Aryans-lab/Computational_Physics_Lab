# Week 05a — Bracketing: Bisection & Regula Falsi

Root-finding begins with the only method that is *guaranteed*: if a
continuous function changes sign on $[a, b]$, a root must be in there, and
bisection shrinks the interval by half each step. No derivatives, no
guessing — just the intermediate value theorem.

## Questions

| # | Driver | What it does | Headline result |
|---|---|---|---|
| 1 | `question1.py` | Library validation on $f(x) = x^2 - 2$: both methods vs $\sqrt 2$, plus the auto-bracketing demo starting from an interval with *no* sign change | both land within $10^{-10}$; `bracket_root` expands $[3,4]$ until it brackets the root |
| 2 | `question2.py` | $f(x) = \log(x/2) - \sin(5x/2)$ on $[1.0, 1.5]$, accuracy $10^{-6}$ | root ≈ $1.40193$ by both methods |
| 3 | `question3.py` | $f(x) = -x - \cos x$, *starting bracket $[2, 3]$ which contains no sign change* — the driver must first find the right interval | expansion lands on $[-1.08, 3.0]$; root $-0.7390851332\ldots$ |

## What's worth reading

- **The iteration count *is* the math.** Bisection needs
  $\lceil \log_2((b-a)/\epsilon) \rceil$ steps — the driver prints the
  measured count next to the $\log_2$ prediction so the linear rate is
  visible, not just claimed.
- **Regula falsi's quiet flaw.** It aims a secant instead of bisecting,
  usually converging faster — *unless* one endpoint sticks, in which case it
  creeps. Problem 3's wide bracket exposes that behaviour (100 iterations
  for $10^{-8}$), which is the honest reason "bisection is slow" is
  immediately followed by "Newton is not."
- **Bracket expansion.** When the given interval has no sign change, the
  routine pushes out the side with smaller $|f|$ (step growing 10% each
  push) until it does — a small piece of robustness that saves a lot of
  debugging.

## Running

```bash
python question2.py   # zero arguments: writes output/q2_output.txt
```

`mylib.py` here is a verbatim bundle of the canonical library
(`make bundle` regenerates it).
