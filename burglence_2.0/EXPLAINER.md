# The complete explainer

Everything in this project, from the first idea to the last number, written so
that you can defend any line of it without notes.

**How to use this document.** Read Parts 1–2 once for the story, then Parts 3–4
with pen and paper (do the derivations yourself — that is the only way they
become yours), then Parts 5–7 with the code open beside you. Part 10 is the
viva. If you learn one page, learn §3.3–3.4 (the two formulas) and §7 (why
your numbers are a few per cent off).

Contents
1. How the project came about
2. The literature, and what was taken from it
3. The physics, derived from nothing
4. The numerics, derived from nothing
5. The code, function by function
6. The three experiments, number by number
7. Error analysis: why 5 % and not 0.5 %
8. What is yours and what is textbook
9. The talk
10. Thirty viva questions
11. Traps

---

# 1. How the project came about

## 1.1 The brief

The DIY component (20 marks) asked for a problem that uses the ideas of the
course: self-avoiding walks, adaptive Runge–Kutta, SVD, chaos, Fourier
analysis, and the one-dimensional Navier–Stokes equation, singly or combined.
Grading is on coding, plotting and interpretation — not on physics breadth.

## 1.2 The eliminations, in order

**Chaos was out immediately.** About 40 % of the batch is doing the logistic
map or Lorenz. More importantly, a chaos project has no exact answer to check
against: you plot two trajectories that separate and the plot *is* the result.
There is nothing to be wrong about, and nothing to be right about either.

**SVD was out** for a similar reason: it is a linear-algebra technique looking
for a data set, and no natural data set in this course needs it.

**Self-avoiding walks** are interesting but are pure Monte Carlo — no PDE, no
spectra, and the interesting part (the exponent ν ≈ 0.588) needs enormous
statistics to pin down.

**Adaptive Runge–Kutta** is a tool, not a problem. On its own it is a
glorified ODE solver.

That leaves two: **Fourier analysis** and **1-D Navier–Stokes** — and they fit
together, which is the whole design.

## 1.3 The realisation that made the project

The 1-D Navier–Stokes equation looks impossible. Here is the algebra that
makes it trivial, and this is the single most important thing to be able to say
out loud:

> In one dimension, incompressibility reads ∂u/∂x = 0, because the divergence
> of a one-component velocity field is just ∂u/∂x. If ∂u/∂x = 0 then the
> advection term u·∇u = u ∂u/∂x vanishes — and so does the pressure gradient,
> which in 1-D is the only thing that could still drive the flow.

So the full 1-D momentum equation collapses to almost nothing. Reinstate one
term — a pressure gradient, or equivalently a body force — and you get

  u_t + u u_x = ν u_xx

which is the **viscous Burgers equation**. It is the only non-trivial
one-dimensional reduction of Navier–Stokes, it is a single PDE, and — this is
the point — **it keeps the two effects that make fluid dynamics hard**:

* `u u_x` is nonlinearity. It steepens a smooth wave until it breaks.
* `ν u_xx` is viscosity. It refuses to let the wave break.

A shock is exactly the compromise between them. So one equation gives you a
genuinely fluid-dynamical phenomenon, at a cost of one PDE and 250 lines of
code.

## 1.4 Why Fourier analysis is the right second keyword

A shock is a structure with a *long tail in wavenumber*: it is smooth
everywhere except at one point, so its Fourier coefficients decay slowly, like
1/k. Smooth functions decay exponentially; shocks decay algebraically. That
gives a prediction — **E(k) ∝ k⁻²** — which is

* a second, completely independent check on the same simulation,
* a one-line NumPy job (`np.fft.rfft`), so it costs no extra code,
* and the standard language of turbulence, so the result connects to the
  literature rather than being a curiosity.

The two keywords are therefore not bolted together: the PDE produces the
structure, and the Fourier transform measures it.

## 1.5 The five design decisions nobody will ask about but everybody should

| decision | choice | why |
|---|---|---|
| initial condition | `u₀ = sin x` | its breaking time is analytically known (`t_s = 1`), and after breaking it becomes a sawtooth — which is exactly what produces a `k⁻²` spectrum. One initial condition serves both experiments. |
| boundary conditions | `u = 0` at both ends | `sin x` is odd about `x = 0`; Burgers preserves oddness; so the exact solution vanishes there. The boundary condition is *exact*, not a convenience. (§3.8) |
| viscosity | `ν = 0.01` | compromise: shock thickness `δ = 4ν/U ≈ 0.044`, which is 14 grid cells at `N = 2048` — resolved — while the run still takes seconds. |
| resolution | `N = 2048` | 14 points across the shock. Below ~5 the thickness measurement degrades (measured, §6.2). |
| time step | `dt = 0.6 Δx / max|u|` | CFL with a safety factor. For the Crank–Nicolson solver diffusion imposes no limit, so only advection matters. |

The one thing that is *not* a free choice: `dt < Δx²/2ν` for the explicit
solver. At `N = 2048` that is `dt < 4.7 × 10⁻⁵`, i.e. 60 000 steps — which is
exactly why the implicit solver exists. (§4.3–4.4)

---

# 2. The literature, and what was taken from it

Three references. Learn the *shape* of each, not the page numbers.

## 2.1 Burgers (1948), *A mathematical model illustrating the theory of turbulence*

**What it is.** A seven-page paper introducing `u_t + u u_x = ν u_xx` as the
simplest equation with both a nonlinear transfer of energy and a dissipative
term. Burgers' motivation was turbulence: real turbulence has an energy
cascade, and he wanted the smallest model where energy is fed in at large
scales by nonlinearity and removed at small scales by viscosity.

**What we take.** The equation, and the framing sentence: it is the smallest
model that contains both effects. That is the sentence to open the talk with.

**What we do not take.** Everything statistical — the paper's later sections
are about turbulence spectra and shock statistics. Out of scope.

## 2.2 Whitham (1974), *Linear and Nonlinear Waves*, chapters 2–4

**What it is.** The standard textbook. Chapter 2 is the method of
characteristics; chapter 3 is wave breaking (including the formula for the
breaking time); chapter 4 is the effect of viscosity — the tanh front and its
thickness.

**What we take — and this is what the project "replicates".**

* the breaking time `t_s = -1/min u₀'`, derived in four lines (§3.2), and
  tested to `1.00` against the measured strain divergence;
* the exact steady shock `u = -U tanh[U(x-x₀)/2ν]` — *verified by substitution*
  in §3.4, then compared with the numerical profile to 2 % of `U`;
* its thickness `δ = 4ν/U`, measured at four viscosities, scaling as `ν^1.00`.

Figure 2(c) and 2(d) of the report are effectively Whitham's chapter 4,
computed from scratch and checked against the formula.

## 2.3 Bec & Khanin (2007), *Burgers turbulence*, Physics Reports 447

**What it is.** A 66-page review. Only the first section matters here: it sets
up the problem, defines the energy spectrum, and states that for random smooth
initial data the spectrum develops a `k⁻²` range once shocks form.

**What we take.** The prediction `E(k) ∝ k⁻²`, and the vocabulary (inertial
range, dissipative cut-off) that lets you describe figure 3 in the language a
physicist expects.

**What we do not take — and be explicit about this if asked.** The review's
real content is multifractal structure functions, `ζ_p = min(p,1)`, shock
statistics, Lagrangian limits. That is genuinely graduate-level and it is
deliberately absent here. Saying "I read section 2 and used the spectrum
prediction; the rest is beyond the course" is a strong answer, not a weak one.

## 2.4 How to describe the reading honestly

> "Burgers' 1948 paper gave me the equation and the motivation. Whitham's
> chapters 2–4 gave me the two analytic results I test against — the breaking
> time and the tanh shock with its `4ν/U` thickness. Bec and Khanin's review
> gave me the spectral prediction `E(k) ∝ k⁻²`. I deliberately did not use the
> multifractal part of the review; it needs machinery well beyond this course."

Three sentences, each true, each tied to a figure.

---

# 3. The physics, derived from nothing

## 3.1 Reading the equation

```
u_t  +  u u_x  =  ν u_xx
```

Take each term separately.

* `u_t` — how fast the velocity at a fixed point changes.
* `u u_x` — transport: the fluid carries its own velocity profile along with
  it. Where `u` is positive, the profile is swept to the right; where the
  profile is decreasing, faster fluid catches up with slower fluid ahead and
  the profile *steepens*. This is the term that makes shocks.
* `ν u_xx` — diffusion: velocity smooths out, curving the profile towards a
  straight line. This is the term that prevents shocks.

Set `ν = 0` and nothing stops the steepening. Keep `ν > 0` and the steepening
is arrested at the scale where the two terms balance.

## 3.2 Characteristics and the breaking time

For `ν = 0`, the equation is a first-order PDE and can be solved exactly by
following the fluid.

**Step 1 — the characteristic equation.** Along a curve `x(t)`,

```
du/dt = u_t + (dx/dt) u_x .
```

Comparing with `u_t + u u_x = 0`, the two agree if `dx/dt = u`. Along such a
curve `du/dt = 0`, so **`u` is constant along characteristics** and each
characteristic is a straight line moving at its own fluid velocity:

```
x(t) = x₀ + u₀(x₀) t .
```

**Step 2 — the slope.** Differentiate the PDE with respect to `x`:

```
(u_x)_t + u (u_x)_x + u_x² = 0 .
```

That is again a transport equation along the *same* characteristics, but with
a source term:

```
d(u_x)/dt = -u_x² .
```

**Step 3 — solve it.** Separating variables,

```
u_x(t) = u₀'(x₀) / [1 + t u₀'(x₀)] .
```

Two things follow.

* Where `u₀' > 0` the slope decays — the profile flattens.
* Where `u₀' < 0` the denominator shrinks and the slope **blows up** at
  `t = -1/u₀'(x₀)`.

**Step 4 — the breaking time.** The first point to break is the one with the
most negative slope:

```
t_s = -1 / min_x u₀'(x) .
```

For `u₀ = sin x`, `min u₀' = -1` at `x = π`, so `t_s = 1`.

This is the whole of the first experiment's theory. It is four lines, it is
exact, and the simulation reproduces it to two decimal places.

## 3.3 What actually happens at `t_s`

At `t_s` the slope is infinite: the solution of the PDE stops being a function
(it becomes multivalued). Physically the fluid on the left overtakes the fluid
on the right and, in a real gas, a **shock** forms — a thin region where the
velocity jumps.

In the inviscid equation the shock is a genuine discontinuity. A discontinuity
is not a classical solution, so one extends the idea: a **weak solution**, a
function that satisfies the equation in an integrated (distributional) sense
everywhere except on the jump, where it must satisfy the **Rankine–Hugoniot**
condition `s = (u_L + u_R)/2` — the shock speed is the average of the
velocities on the two sides. (For the symmetric shock in our problem,
`u_L = +U`, `u_R = -U`, so the shock sits still: `s = 0`, which is exactly what
figure 2(a) shows.)

Weak solutions are not unique; the physical one is picked out by the **entropy
condition** — characteristics must run *into* the shock, not out of it.

With `ν > 0` none of this is needed: the solution stays smooth for all time,
the discontinuity is replaced by a thin but finite front, and the answer is
unique. That is why the project keeps `ν` in the equation even though it is
small. The limit `ν → 0` of the viscous solution selects the physically
correct weak solution — this is the "vanishing viscosity" idea.

## 3.4 The exact shock profile (derive it, do not quote it)

Claim: for any `U > 0` and any `x₀`,

```
u(x) = -U tanh[ U (x - x₀) / (2ν) ]
```

solves the equation with `u_t = 0`.

**Verify by substitution.** With `a = U/(2ν)`:

```
u      = -U tanh(a(x-x₀))
u_x    = -U a sech²(a(x-x₀)) = -(U²/2ν) sech²(...)
u_xx   = -(U²/2ν) · (-2a sech² tanh) = (U³/2ν²) sech² tanh
```

Left-hand side: `u u_x = (-U tanh)·(-(U²/2ν) sech²) = (U³/2ν) sech² tanh`.

Right-hand side: `ν u_xx = ν · (U³/2ν²) sech² tanh = (U³/2ν) sech² tanh`.

They are identical, for every `U`, `x₀`, `ν`. So the tanh is an exact steady
solution — not an approximation — and `u_t = 0` as assumed.

**Its thickness.** The maximum slope is `U²/2ν` (set `sech² = 1` at the
centre). The jump across the shock is `2U`. The natural definition of the
thickness is

```
δ = (jump) / (maximum slope) = 2U / (U²/2ν) = 4ν/U .
```

That factor of 4 is not magic; it is just "jump divided by steepest slope" for
a tanh. If you are asked where 4 comes from, that is the answer.

**What it means physically.** Thicker fluid ⇒ thicker shock. Faster fluid
(bigger `U`) ⇒ thinner shock. Since `ν` has units of `length²/time` and `U` of
`length/time`, `ν/U` is a length — the only length in the problem. Any
thickness must be a constant times it, and for a tanh the constant is 4.

## 3.5 Energy balance (why the flow decays, and why that matters)

Multiply the equation by `u` and integrate over the domain:

```
d/dt ∫ ½u² dx = -ν ∫ u_x² dx + [boundary terms]
```

The nonlinear term integrates to `-∫ u³/3 dx` evaluated at the ends — zero for
periodic boundaries, and zero here because `u = 0` at both ends. So:

> **Energy is never created; it is only ever destroyed by viscosity, at a rate
> proportional to `∫u_x²`.**

Two consequences you can use in the viva:

1. A shock, where `u_x` is huge, is a *very* efficient sink of energy. The
   fluid's kinetic energy drains away fastest exactly where the shock is.
2. With our boundaries, energy also leaves through the ends, so the amplitude
   `U` slowly decays. That is why figure 2(b) plots `U(t)²/2ν` rather than a
   constant, and why `U` at `t = 3` is 0.736 rather than 1.

## 3.6 The late-time shape: a sawtooth

After breaking, the profile is a straight line on either side of the shock
(with our boundaries it relaxes towards zero at the ends). The reason is
geometric: in regions where `u_x` is small compared with `1/t`, the
characteristics are nearly parallel, so the profile is carried along almost
unchanged — a straight line stays a straight line. Only near the shock, where
`u` changes abruptly, does anything interesting keep happening.

So the asymptotic solution is a **sawtooth**: linear segments joined by a
sharp front. Hold on to that, because the whole of experiment 3 follows from
it.

## 3.7 Why a sawtooth has a `k⁻²` spectrum

Write the Fourier series `u(x) = Σ c_k e^{ikx}`. The coefficients are

```
c_k = (1/2π) ∫ u(x) e^{-ikx} dx .
```

Integrate by parts:

```
c_k = (1/2π) [ u e^{-ikx}/(-ik) ]_ends  +  (1/(2πik)) ∫ u'(x) e^{-ikx} dx .
```

The first term vanishes (periodicity / zero ends). The second says:
**`c_k` is `1/k` times the Fourier coefficient of `u'`**. If `u'` is a
well-behaved function (which it is, except at the shock, where `u'` is a spike),
then `c_k ∝ 1/k`, and therefore

```
E(k) = |c_k|² ∝ k⁻² .
```

For a *smooth* function, you can integrate by parts over and over and the
coefficients fall faster than any power. **The algebraic `k⁻²` tail is the
fingerprint of a jump.** That single sentence is the interpretation of
figure 3.

Where does it stop? Where viscosity erases the jump, i.e. at wavenumbers
`k ≈ 1/δ`. Halving `ν` halves `δ` and therefore doubles the cut-off — measured
as 31 → 63 in figure 3(b).

## 3.8 Why the boundary conditions are exact

`u₀ = sin x` is odd about `x = 0`: `u₀(-x) = -u₀(x)`.

If `u` is odd, then `u_x` is even, `u u_x` is odd, and `u_xx` is even. So the
whole right-hand side `-u u_x + ν u_xx` is **odd**, and the equation preserves
oddness exactly. An odd function is zero at `x = 0`. By periodicity it is also
zero at `x = 2π`.

So `u(0,t) = u(2π,t) = 0` is not an approximation or a convenience: it is what
the exact solution does. (The test `test_mean_stays_zero` verifies the
consequence: the mean of `u` stays at zero to machine precision.)

---

# 4. The numerics, derived from nothing

## 4.1 The grid

`N` intervals, `N+1` points, `x_j = j Δx` with `Δx = 2π/N`. The two end points
are never updated — they are the boundary condition. Everything below happens
on the interior points `1 … N-1`.

## 4.2 Difference formulas (Taylor, not memory)

Expand about `x_j`:

```
u_{j+1} = u_j + Δx u' + Δx²/2 u'' + Δx³/6 u''' + …
u_{j-1} = u_j - Δx u' + Δx²/2 u'' - Δx³/6 u''' + …
```

Subtract: `u_{j+1} - u_{j-1} = 2Δx u' + O(Δx³)`, so

```
u' ≈ (u_{j+1} - u_{j-1}) / (2Δx)      error O(Δx²)
```

Add: `u_{j+1} + u_{j-1} = 2u_j + Δx² u'' + O(Δx⁴)`, so

```
u'' ≈ (u_{j+1} - 2u_j + u_{j-1}) / Δx²   error O(Δx²)
```

Both are second order. Those two lines are `advection()` and `diffusion()` in
`burgers.py`.

## 4.3 The explicit scheme, and where `Δx²/2ν` comes from

```
u_j^{n+1} = u_j^n + Δt [ -u_j^n (u_{j+1}^n - u_{j-1}^n)/(2Δx)
                        + ν (u_{j+1}^n - 2u_j^n + u_{j-1}^n)/Δx² ]
```

One line of code. Now the stability, by von Neumann: substitute a Fourier mode
`u_j^n = g^n e^{iθj}` (θ is the dimensionless wavenumber, `0 ≤ θ ≤ π`). Every
term just multiplies, and

```
g(θ) = 1 - i C sin θ - 2r(1 - cos θ),
C = |u| Δt/Δx   (CFL number),     r = ν Δt/Δx² .
```

Stability means `|g| ≤ 1` for every θ. Take the worst case for the diffusion
part, `θ = π` (the grid-scale mode, `cos θ = -1`, `sin θ = 0`):

```
g = 1 - 4r    ⇒    |1 - 4r| ≤ 1    ⇒    0 ≤ r ≤ ½ .
```

That is the condition

```
Δt ≤ Δx² / (2ν) .
```

**Why it hurts.** Halve `Δx` and `Δt` must fall by a factor **four**; total
work scales as `1/Δx³`. At `N = 2048`, `ν = 10⁻²` the limit is
`Δt < 4.7 × 10⁻⁵`, i.e. ~60 000 steps for `t = 3`. (The centred advection term
adds a CFL-type restriction `Δt ≲ Δx/|u|`, but at these parameters the
diffusion limit is much the stricter one.)

Experiment 1(c) tests this directly: at `0.9 ×` the limit the solution decays
smoothly, at `1.1 ×` it reaches 10⁶ within 200 steps. The measured boundary is
where the theory puts it.

## 4.4 Crank–Nicolson: removing the restriction

The explicit scheme evaluates diffusion at the *old* time. The implicit idea
is to evaluate it at the *new* time; the symmetric idea, which is second order
in time, is to average:

```
u^{n+1} - u^n = (Δt ν / 2) (D2 u^{n+1} + D2 u^n)
```

where `D2` is the second-difference operator. Rearranged:

```
(I - r/2 D2) u^{n+1} = (I + r/2 D2) u^n ,      r = νΔt/Δx² .
```

This is a **linear system** — everything on the left is known except `u^{n+1}`.
But look at the matrix: each row involves only `u_{j-1}, u_j, u_{j+1}`. It is
**tridiagonal**:

```
⎡ 1+r   -r/2    0    …  ⎤ ⎡u₁⎤     ⎡ rhs₁ ⎤
⎢ -r/2  1+r   -r/2   …  ⎥ ⎢u₂⎥  =  ⎢ rhs₂ ⎥
⎢  0   -r/2    1+r   …  ⎥ ⎢u₃⎥     ⎢ rhs₃ ⎥
⎣  …                    ⎦ ⎣ … ⎦     ⎣  …   ⎦
```

Its von Neumann amplification factor is

```
g(θ) = [1 - r(1-cos θ)] / [1 + r(1-cos θ)] ,
```

whose modulus is ≤ 1 **for every `r`**. So there is no diffusion restriction
at all — and the test `test_cn_survives_a_huge_time_step` runs at `50 ×` the
explicit limit and stays finite.

The price is that each step now requires solving a linear system. Which is
where the Thomas algorithm comes in.

## 4.5 The Thomas algorithm

Solving a full `N × N` system costs `O(N³)` (Gauss) or `O(N²)` (with LU). A
tridiagonal system can be solved in `O(N)` by simple elimination.

Write the system as `a_i u_{i-1} + b_i u_i + c_i u_{i+1} = d_i`.

**Forward sweep** — eliminate the sub-diagonal:

```
for i = 1 … n-1:
    m       = a_i / b_{i-1}
    b_i    -= m c_{i-1}
    d_i    -= m d_{i-1}
```

**Back substitution:**

```
u_n = d_n / b_n
for i = n-1 … 1:  u_i = (d_i - c_i u_{i+1}) / b_i
```

Roughly 8 arithmetic operations per unknown — `O(N)` in both time and memory.
That is `mylib.thomas_algorithm`, called twice per step in `step_cn`.

**Is it safe?** The algorithm fails only if a pivot vanishes. Here
`|b| = 1 + r > 2·|r/2| = |a| + |c|`, so the matrix is *strictly diagonally
dominant* and no pivot can vanish. Good answer if asked about stability of the
linear solve.

## 4.6 Splitting, and why the halves matter

The full step has two pieces — advection `A` and diffusion `D` — which do not
commute (`AD ≠ DA`). You cannot apply them simultaneously with two separate
black boxes, so you split.

**Lie splitting** (diffusion then advection, once each):

```
u^{n+1} = e^{Δt A} e^{Δt D} u^n      local error O(Δt²)  ⇒  global O(Δt)
```

**Strang splitting** (half, full, half):

```
u^{n+1} = e^{Δt/2 D} e^{Δt A} e^{Δt/2 D} u^n   local error O(Δt³) ⇒ global O(Δt²)
```

The symmetry is what does it: the error of the first half-step is cancelled by
the second to one order higher.

**I measured this.** With a plain single split, the convergence orders came
out `1.84, 1.62, 1.43` — degrading towards first order exactly as the theory
predicts. With the symmetric split they are `2.02, 2.00, 2.02`. That is a
genuinely nice thing to be able to say in a viva, and it is what `step_cn`
implements:

```python
u = _implicit_diffusion(u, 0.5 * r)      # first half of the diffusion
k1 = advection(u, dx)                    # advection, predictor-corrector
k2 = advection(u + dt * k1, dx)
u = u + 0.5 * dt * (k1 + k2)
return _implicit_diffusion(u, 0.5 * r)   # second half
```

The advection step is *predictor–corrector* (Heun): take a full Euler step,
evaluate the slope there, then average the two slopes. That is what makes it
second order as well.

## 4.7 How you measure the order of a scheme

You need an answer to compare with. There is no exact solution after breaking,
so we use **self-convergence**: compute the solution on a fine grid
(`N = 4096`) and treat it as exact, then compute on `64, 128, 256, 512` and
measure

```
E(N) = rms( u_N - u_ref ) .
```

Second order means `E ∝ N⁻²`, i.e.

```
order = log₂[ E(N)/E(2N) ] .
```

Two subtleties worth knowing:

* the reference must be *more* accurate than the runs being tested, otherwise
  the measured order saturates. We ran the reference with `Δt/4` as well as
  `Δx/2`, so its error is far below the smallest `E` measured;
* the coarse grid points are a subset of the fine ones (`N` divides 4096), so
  no interpolation is needed — just take every `4096/N`-th point.

Result: `2.02, 2.00, 2.02`.

## 4.8 The three diagnostics

**`max_slope`** — `max_j |u_{j+1}-u_j|/Δx`. This is the quantity that diverges
at `t_s`; it is what figure 2(b) plots.

**`shock_width`** — `(max u − min u) / max_slope`. For a tanh this equals
`4ν/U` exactly (§3.4). It assumes a *single* shock, which is true here.

**`energy_spectrum`** — take `u`, subtract the mean, `c = np.fft.rfft(u)/n`,
`E = |c|²`, `k = 2π · fftfreq`. The absolute normalisation is irrelevant
because only the *slope* is quoted.

**`fit_power_law`** — fit `log E = a + b log k` over a stated window with
`mylib.poly_fit`, and quote `b`. The window `4 ≤ k ≤ 12` is chosen *below* the
dissipative cut-off (31), never at the point where the curve bends.

---

# 5. The code, function by function

`burgers.py` (257 lines) is the whole physics. Nothing else contains any.

| function | what it does | one thing to say about it |
|---|---|---|
| `make_grid(N)` | `N+1` points on `[0, 2π]` | endpoint included, `Δx = 2π/N` |
| `initial_sine(x)` | `sin x` | the only initial condition used |
| `advection(u, dx)` | `-u u_x`, centred | second order, zero on boundaries |
| `diffusion(u, dx, nu)` | `ν u_xx`, centred | second order |
| `step_explicit` | one Euler step | stable only if `Δt < Δx²/2ν` |
| `_implicit_diffusion(u, r)` | one CN step, one Thomas solve | where `mylib.thomas_algorithm` is called |
| `step_cn` | half-D, full-A, half-D | second order, no diffusion limit |
| `stable_dt` | the `Δt` rule | CFL, plus `Δx²/2ν` only for explicit |
| `solve` | the time loop | stores profiles at requested times |
| `max_slope` | diagnostic | diverges at `t_s` |
| `shock_width` | diagnostic | `jump / max slope` |
| `energy_spectrum` | diagnostic | one FFT |
| `fit_power_law` | diagnostic | `mylib.poly_fit` on logs |
| `l2_difference` | error norm | used by the convergence test |

The core, verbatim, is about twenty lines:

```python
def advection(u, dx):
    f = np.zeros_like(u)
    f[1:-1] = -u[1:-1] * (u[2:] - u[:-2]) / (2.0 * dx)
    return f

def diffusion(u, dx, nu):
    f = np.zeros_like(u)
    f[1:-1] = nu * (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dx ** 2
    return f

def step_explicit(u, dx, dt, nu):
    return u + dt * (advection(u, dx) + diffusion(u, dx, nu))
```

`u[1:-1]` is the interior; `f[0] = f[-1] = 0` keeps the boundaries fixed at
zero. That's it — that is the entire explicit solver.

---

# 6. The three experiments, number by number

## 6.1 Experiment 1 — is the code right? (`fig01_sanity`)

**(a) Diffusion only.** Set the nonlinear term aside: `u_t = νu_xx` has the
exact solution `exp(-νt) sin x`. The diffusion half-step is the library's
Crank–Nicolson, so this validates mylib and our use of it.
*Result: max error `7.95 × 10⁻⁶`.*

**(b) Convergence.** Reference `N = 4096` with `Δt/4`; test runs at
`N = 64…512`; `t = 0.5` (before breaking, so the solution is smooth).
*Result: orders `2.02, 2.00, 2.02`.*

**(c) Stability.** `N = 512`, `ν = 0.01` ⇒ `Δx²/2ν = 7.53 × 10⁻³`. Run at
`0.9×` and `1.1×`. *Result: `0.9×` decays to `0.788`; `1.1×` blows up to
infinity within a few hundred steps.*

**What to say:** "Three independent checks — an exact solution, a measured
order, and a stability boundary — all agree with theory before I show a single
result."

## 6.2 Experiment 2 — the shock (`fig02_shock`)

Setup: `ν = 0.01`, `N = 2048`, `Δt = 0.6Δx`, `T = 3`, `u₀ = sin x`.

**(a) Profiles at t = 0, 0.5, 1, 2, 3.** Steepening, a vertical front at
`x = π` appearing at `t ≈ 1`, then a silhouette that stops changing.

**(b) `max|u_x|` against time.** Before `t = 1` it tracks `1/(t_s − t)`;
afterwards it follows `U(t)²/2ν`. At the final time: measured **28.4** against
`U²/2ν = 27.1` with `U = 0.736`.

**(c) The profile at `t = 3` against `−U tanh[U(x−π)/2ν]`.** The two curves
are on top of each other; the largest deviation over the plotted window is
`2.0 × 10⁻²`, i.e. **2 % of `U`**.

**(d) Thickness against viscosity.** Four runs, resolution scaled so the shock
always spans ≥ 7 cells:

| ν | N | δ measured | 4ν/U | points in shock |
|---|---|---|---|---|
| 0.040 | 512 | 0.1690 | 0.1851 | 15.1 |
| 0.020 | 512 | 0.0843 | 0.0885 | 7.2 |
| 0.010 | 1024 | 0.0422 | 0.0433 | 7.1 |
| 0.005 | 2048 | 0.0211 | 0.0214 | 7.0 |

Fitted slope: **1.000**. The ratios are `0.91, 0.95, 0.97, 0.99` — they
approach 1 from below, which is the honest sign of a resolution effect (§7).

## 6.3 Experiment 3 — the spectrum (`fig03_spectrum`)

Setup: `ν = 0.01` and `0.005`, `N = 2048`, `T = 2`, same initial condition.

**(a) `E(k)` on log–log.** Straight lines. Fitted over `4 ≤ k ≤ 12`:
**−2.060** and **−2.027**.

**(b) `E(k)k²`.** Flat where the power law holds (a plateau at ≈ 0.09), then a
knee. Cut-off (where it has fallen 30 % below the plateau): **k ≈ 31** at
`ν = 0.01`, **k ≈ 63** at `ν = 0.005`. Halving `ν` doubles the cut-off, as
`1/δ` predicts.

**What to say:** "Viscosity does not change the slope of the spectrum; it
changes where the spectrum ends."

---

# 7. Error analysis: why 5 % and not 0.5 %

This is the part that separates a careful student from a careless one. Know
these four numbers and their causes.

### 7.1 Thickness is 5 % low (`0.0518` vs `0.0544`)

Two effects, in opposite directions.

*Discretisation pushes it the wrong way.* For a tanh sampled at spacing `h`,
the steepest discrete slope is `2U tanh(a)/h` with `a = Uh/4ν`, and the ratio
to the true slope is `tanh(a)/a < 1`. Since `δ = jump/slope`, discretisation
makes `δ` slightly *too large* — by only 0.2 % here (`a ≈ 0.07`), so this is
not the explanation.

*The profile is not yet exactly steady.* The measured maximum slope (28.4) is
5 % **higher** than `U²/2ν` (27.1) for the instantaneous amplitude, and `δ` is
low by the same 5 %. The front is still relaxing: with our boundaries the
amplitude is decaying, so the solution is quasi-steady rather than steady. The
tanh comparison in panel (c) confirms it — the shapes agree to 2 % of `U`, not
to 0.02 %.

**Viva answer:** "Five per cent low, and it is not a bug: the analytic
`4ν/U` is for an exactly steady shock, while ours is still relaxing because
the amplitude decays at the boundaries. The ratio improves as `ν` falls —
0.91, 0.95, 0.97, 0.99 — and the *scaling* is `ν^1.000`, which is the robust
result."

### 7.2 Resolution criteria

Two different requirements, often confused:

* to *represent* a shock in a plot, ~3–4 points across `δ` is enough;
* to *measure* its thickness from a maximum slope, you need more, because the
  maximum slope is a two-point quantity. The sweep table shows the bias at
  7 points per shock is ~3 %.

Rule of thumb used here: **at least 7 points across `δ`**, i.e.
`Δx ≤ δ/7 = 4ν/(7U)`.

### 7.3 The spectral fit window

`4 ≤ k ≤ 12` is not arbitrary: `k = 12` is well below the cut-off (31), and
`k = 4` is far above the box scale. Fitting beyond the knee would drag the
slope down — that is why the same code gives `−2.3` if you fit to `k = 40` at
`ν = 0.01`. Quote the window, always.

### 7.4 What is *not* an error

The strain **overshoot** at `t ≈ 1.5` (the dots rise above the plateau before
settling): that is real physics — right after breaking the front is thinner
than the asymptotic one because the amplitude is still near its maximum, and
`U²/2ν` is larger then. The red curve uses the instantaneous `U(t)`, so the
dots sit on it.

---

# 8. What is yours and what is textbook

Be precise about this; vagueness is what gets punished.

**Textbook / from the references (say so):**
the equation (Burgers 1948); the breaking-time derivation; the tanh profile and
`δ = 4ν/U` (Whitham ch. 4); `E(k) ∝ k⁻²` for a sawtooth (Bec–Khanin §2);
Crank–Nicolson, the Thomas algorithm, von Neumann analysis, Strang splitting.

**Yours (this is a real list):**
* the choice of problem, and the reduction argument that makes 1-D
  Navier–Stokes *be* Burgers;
* the design log — `sin x`, Dirichlet-via-oddness, `ν = 0.01`, `N = 2048`;
* the implementation: two solvers, the diagnostics, the whole pipeline in one
  file and three scripts;
* the *measurements*: every number in the tables, produced by your runs;
* the error analysis (§7) — including noticing that the plain split degrades
  to first order and fixing it with a symmetric split, and running both to
  show it;
* the comparison that shows the tanh profile and the `4ν/U` law hold to a few
  per cent, and the spectral cut-off moving as `1/δ`.

That is a normal, respectable undergraduate project: textbook physics,
independently implemented and independently *checked*. Nobody expects a new
equation; they expect you to know which parts are which.

---

# 9. The talk

Full spoken script: `TALK_SCRIPT.md`. Structure and timing:

| slide | min | message |
|---|---|---|
| 1 | 0.5 | title |
| 2 | 2 | one equation, two keywords; what a shock is |
| 3 | 2 | two solvers; the `Δx²` limit; what each is checked against |
| 4 | 2.5 | figure 1: exact solution, order 2, stability boundary |
| 5 | 2.5 | figure 2: `t_s = 1`, tanh, `4ν/U`, `ν^1.00` |
| 6 | 2.5 | figure 3: `k⁻²`, cut-off moves as `1/δ` |
| 7 | 2 | takeaway + limitations (say them before you are asked) |
| 8 | 1 | references, questions |

The only sentences that must be word-perfect are the three "money lines":

1. "In one dimension incompressibility forces `u_x = 0`, so advection and the
   pressure gradient both vanish — what survives is Burgers."
2. "The nonlinear term steepens until `t_s = 1`; viscosity replaces the
   discontinuity with a `tanh` front of thickness `4ν/U`."
3. "A shock is a jump, a jump gives `c_k ∝ 1/k`, so the spectrum is `k⁻²`;
   viscosity decides only where that spectrum ends."

---

# 10. Thirty viva questions

**On the problem**

1. *Why this topic?* — One PDE that contains a real fluid phenomenon and has
   analytic answers for every quantity I measure; and it pairs naturally with
   Fourier analysis, which is the second keyword.
2. *Where is Navier–Stokes?* — 1-D incompressibility kills advection and
   pressure; reinstating a pressure gradient gives Burgers.
3. *Why not 2-D/3-D?* — No exact solutions, orders of magnitude more compute,
   and the marks are for coding and interpretation, not for resolution.
4. *Is this turbulence?* — No. It is the standard model of a compressible 1-D
   flow; the `k⁻²` spectrum is the one piece of turbulence vocabulary that
   survives in 1-D.
5. *Is it chaos?* — No: no sensitive dependence, no positive Lyapunov
   exponent; the flow relaxes to a steady shock.

**On the physics**

6. *Derive `t_s`.* — §3.2, four lines.
7. *Why `t_s = 1` for `sin x`?* — `min u₀' = -1` at `x = π`.
8. *What happens at `t_s`?* — The classical solution becomes multivalued; a
   shock forms; with `ν > 0` it is smoothed to a finite front.
9. *Show the tanh is a solution.* — §3.4, by substitution; both sides equal
   `(U³/2ν) sech² tanh`.
10. *Where does the 4 in `4ν/U` come from?* — jump `2U` divided by peak slope
    `U²/2ν`.
11. *Why does the amplitude decay?* — Viscous dissipation, plus leakage through
    the boundaries; the energy balance is in §3.5.
12. *Why is the profile a sawtooth?* — Characteristics are nearly parallel
    away from the shock, so straight segments are transported unchanged.
13. *Why `k⁻²`?* — One integration by parts turns `c_k` into `1/k` times the
    transform of `u'`, which is a spike. Algebra in §3.7.
14. *Why does the spectrum end?* — Viscosity erases structure below `δ`, i.e.
    above `k ≈ 1/δ`.

**On the numerics**

15. *What order is your scheme, and how do you know?* — Two, measured
    `2.02/2.00/2.02` against a fine reference.
16. *Why does the explicit scheme need `Δx²/2ν`?* — Von Neumann with
    `θ = π`: `g = 1 − 4r`, so `r ≤ ½`. Derivation in §4.3.
17. *Why Crank–Nicolson?* — Second order in time and unconditionally stable
    for diffusion; the `Δx²` limit disappears.
18. *What does the Thomas algorithm do, and how fast?* — `O(N)` tridiagonal
    solve; the matrix is diagonally dominant so no pivoting is needed.
19. *You split the operators — doesn't that cost accuracy?* — A plain split
    does; I measured it (`1.84, 1.62, 1.43`). The symmetric Strang split
    restores second order (`2.02, 2.00, 2.02`).
20. *Why is the advection step predictor–corrector?* — So it is second order
    too; plain Euler would be first order and would drag the whole scheme down.
21. *How is the boundary handled?* — `u = 0` at both ends is exact by oddness
    of `sin x`; the interior indices `1:-1` never touch the ends.
22. *Why does the mean stay zero?* — The exact solution is odd, and the scheme
    preserves that; the test checks it to machine precision.
23. *What is `max_slope`, and what is it for?* — `max|Δu/Δx|`; it diverges at
    `t_s` and saturates at `U²/2ν`.
24. *How do you measure shock thickness?* — jump over peak slope; for a tanh
    that is exactly `4ν/U`.
25. *How did you fit the spectrum?* — `log E` against `log k` over `4 ≤ k ≤ 12`,
    least squares with `mylib.poly_fit`, window below the cut-off.

**On the results**

26. *Is `−2.06` good enough?* — Yes: the prediction is `−2` and the deviations
    are at the few-per-cent level, with the window stated. The cut-off
    behaviour (31 → 63) is a sharper test than the slope itself.
27. *Your thickness is 5 % low — is the code wrong?* — No; §7.1. Quasi-steady
    rather than steady, ratio improving as `ν` falls, scaling exactly 1.
28. *What is your biggest error bar?* — Resolution of the front when measuring
    thickness; quantified in the sweep table as 3 % at 7 points per shock.
29. *What would you do with another month?* — Force the equation so the
    amplitude is stationary; go to `ν = 0.002` at `N = 4096` for a decade more
    of `k⁻²`; add an upwind scheme to show the diffusive alternative.
30. *What did you learn?* — That a numerical answer is only worth something
    if you can say how wrong it is; every number here came with a prediction
    and an error estimate.

---

# 11. Traps

1. **Do not call it turbulence.** It is Burgers turbulence at most; the word
   invites the multifractal questions you did not do.
2. **Do not claim the tanh is exact for your data.** It is exact for the steady
   problem; yours is quasi-steady, 2 % of `U`.
3. **Do not quote a slope without a window.** `−2.06` for `4 ≤ k ≤ 12`.
4. **Do not claim the spectrum is resolved.** One decade, enough to separate
   `k⁻²` from `k⁻³`, nothing subtle.
5. **Do not pretend the splitting is free.** Say it costs a doubling of the
   linear solves and buys second order.
6. **Do not say "the explicit scheme is inaccurate".** It is accurate and
   second order; it is *expensive* because of the `Δx²` restriction.
7. **If you do not know something, say so** and say what you would do to find
   out. That is worth more than a bluff.
