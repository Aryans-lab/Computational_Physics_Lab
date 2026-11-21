# The talk, word for word

Thirteen minutes of speaking, two minutes of slack, five minutes of questions.
Read it out loud twice with a timer — it will feel slower than it reads.

Rules: never read the slide; the slide holds the picture, you hold the
sentence. Say the number, then say what it means, then move on.

---

## Slide 1 — title (0:30)

> "Sir, my project is on the one-dimensional Navier–Stokes equation — which
> turns out to be the Burgers equation — and on what Fourier analysis says
> about the shock that forms in it. Two of the ideas from your list, one
> equation, three figures."

*(Do not read the title aloud. Do not thank anyone yet.)*

---

## Slide 2 — one equation, two keywords (2:00)

> "The equation is `u_t + u u_x = ν u_xx`. In one dimension, incompressibility
> means `∂u/∂x = 0`, because the divergence of a one-component velocity is just
> `∂u/∂x`. If that is zero, both the advection term and the pressure gradient
> vanish. Put back a pressure gradient and what is left is this: the Burgers
> equation. So the one-dimensional Navier–Stokes equation *is* this equation."

> "It has two terms competing. The nonlinear term `u u_x` steepens the wave —
> points moving faster catch up with slower points ahead. The viscous term
> smooths it. A shock is the standoff between them."

> "The picture on the right is the whole project. A sine wave steepens, at
> `t = 1` it breaks, and after that it is a fixed front of finite thickness.
> The thickness is `4ν/U`, and that formula is one of four things I check
> against in the talk."

*(Point at the x = π line, then at the saturated profile at t = 3.)*

---

## Slide 3 — two solvers (2:00)

> "I wrote two solvers. The first is plain forward Euler with centred
> differences — fifteen lines. The problem is stability: for the diffusion
> term you need `Δt < Δx²/2ν`. At my resolution that means sixty thousand
> steps, which is why the second solver exists."

> "The second treats diffusion implicitly, with Crank–Nicolson. That turns each
> step into a tridiagonal linear system, which I solve with the Thomas
> algorithm from `mylib.py` — `O(N)` per step. It removes the `Δt` restriction
> entirely; I tested it at fifty times the explicit limit and it stays finite."

> "Advection is done explicitly with a predictor–corrector step, and the whole
> step is split symmetrically — half diffusion, full advection, half
> diffusion — so that the scheme stays second order. That matters, and I will
> show you the measurement in a moment."

*(Do not explain the Thomas algorithm unless asked. If asked: forward sweep,
back substitution, O(N), diagonally dominant so no pivoting.)*

---

## Slide 4 — figure 1: is it right? (2:30)

> "Three checks before any physics."

> "First, I switch the nonlinear term off. Then the equation is the heat
> equation, which has the exact solution `exp(-νt) sin x`. The numerical
> answer matches it to `8 × 10⁻⁶` — that check validates the library
> diffusion routine and my use of it."

> "Second, with the nonlinear term back on, I compare runs at N = 64 to 512
> against a reference at N = 4096, run with a quarter of the time step. The
> errors fall with measured orders 2.02, 2.00 and 2.02 — second order, as the
> splitting argument predicts. When I used a single unsymmetrical split
> instead, the same measurement gave 1.84, 1.62, 1.43, collapsing towards
> first order — so the symmetric split is not decoration, it is the difference
> between first and second order."

> "Third, the stability boundary: at 0.9 of the von Neumann limit the solution
> decays smoothly; at 1.1 it explodes within two hundred steps. The boundary
> is where theory puts it."

---

## Slide 5 — figure 2: the shock (2:30)

> "Now the physics. Top left: profiles at five times. The steepening is the
> characteristic equation — each point moves at its own speed, so the profile
> folds over itself. It breaks at `t_s = -1/min u₀'`, which for `sin x` is
> exactly 1. Measured: 1.00."

> "Top right: the maximum slope. Before breaking it follows `1/(t_s - t)`;
> after, viscosity caps it at `U²/2ν`. At the end I measure 28.4 against a
> prediction of 27.1."

> "Bottom left: the profile at t = 3 against the exact steady solution,
> `-U tanh[U(x-π)/2ν]`. They agree to two per cent of the amplitude. That
> tanh is an exact solution — you can verify it by substitution in three
> lines."

> "Bottom right: thickness against viscosity, four runs. The prediction is
> `4ν/U`, and I measure `δ ∝ ν^1.000`. The points sit a few per cent below
> the line — bottom right of the panel — and the deficit shrinks as the shock
> gets better resolved. I will come back to that in the limitations."

---

## Slide 6 — figure 3: the spectrum (2:30)

> "This is the Fourier part. A shock is a jump, and a jump has a long tail in
> wavenumber. Formally: integrate `c_k` by parts once and the jump survives as
> a boundary term, so `c_k ∝ 1/k`, and the energy `E(k) = |c_k|²` goes as
> `k⁻²`."

> "Left: measured spectra for two viscosities, fitted over `4 ≤ k ≤ 12`, which
> is below the dissipative range. Slopes `-2.06` and `-2.03`. Both are straight
> lines and they lie on top of each other over the whole inertial range."

> "Right: the same data multiplied by `k²`. If the law is exact this is flat —
> and it is flat, over more than a decade, at about 0.09. Then it falls. The
> knee is where viscosity has eaten the harmonics, at `k ≈ 1/δ`. It sits at 31
> for `ν = 0.01` and at 63 for `ν = 0.005` — halving the viscosity doubles the
> cut-off.

> "So the sentence is: viscosity does not change the slope of the spectrum, it
> changes where the spectrum ends."

---

## Slide 7 — what it means, and limitations (2:00)

> "Three things I take away. One: two schemes with the same order can differ
> enormously in cost, and the difference is entirely the stability limit. Two:
> every number in this project was predicted before it was measured —
> breaking time, thickness, scaling, spectral slope — which is what makes it a
> test rather than a picture. Three: Fourier analysis separates two things
> that look merged in real space; the slope is the shape of the shock, the
> cut-off is the viscosity."

> "Limitations, before you ask. One equation, one shock — this is a test case,
> not a simulation of turbulence. The thickness is 5 per cent low because the
> shock is still relaxing, not exactly steady, and the deficit shrinks as ν
> falls. The spectrum covers one decade, enough to tell `k⁻²` from `k⁻³` and
> nothing finer."

---

## Slide 8 — references (1:00)

> "Three references. Burgers 1948 gave the equation and the motivation.
> Whitham's chapters 2 to 4 gave the breaking time and the tanh shock with its
> `4ν/U` thickness — that is what figures 2(c) and 2(d) check. Bec and Khanin's
> review gave the `k⁻²` prediction for the spectrum; I used the first section
> and deliberately left out the multifractal theory, which is beyond this
> course."

> "Code, figures and every number are in the repository, under
> `burglence_2.0`. Thank you — I am happy to take questions."

---

## Timing discipline

| if you are at | you should be on |
|---|---|
| 2 min | slide 2 finished |
| 4 min | slide 3 finished |
| 6.5 min | slide 4 finished |
| 9 min | slide 5 finished |
| 11.5 min | slide 6 finished |
| 13.5 min | slide 7 finished |
| 14.5 min | slide 8 finished |

If you are more than a minute behind at the 9-minute mark, drop the third
bullet of slide 7 and go straight to limitations.

## The three questions you will certainly get

**"Why is your thickness 5 % off?"** → §7.1 of `EXPLAINER.md`. Quasi-steady,
not steady; ratio improves as ν falls; scaling is exact.

**"Why not just use the explicit scheme?"** → sixty thousand steps versus two
thousand, same accuracy. Show figure 1(c).

**"What is k⁻²?"** → A jump. Integrate by parts once, the jump survives, the
coefficients go as 1/k, the energy as 1/k².
