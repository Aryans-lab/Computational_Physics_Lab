# Three references, and what each one gives this project

The course asks for two or three papers. These are the three; each one is
responsible for a specific claim in the report, so nothing is cited
decoratively.

---

### 1. J. M. Burgers, *A mathematical model illustrating the theory of turbulence*, Advances in Applied Mechanics **1**, 171–199 (1948).

The origin of the equation. Burgers introduces

```
u_t + u u_x = ν u_xx
```

as the simplest model that has both the nonlinear transfer of energy (the
`u u_x` term) and its dissipation (the `ν u_xx` term), and points out that in
one dimension the pressure term drops out completely.

**What we take from it:** the equation itself, and the motivation for studying
it — it is the smallest model that shows steepening and dissipation competing.
**What we do not take:** the statistical/turbulence discussion, which needs
machinery far beyond this project.

---

### 2. G. B. Whitham, *Linear and Nonlinear Waves*, Wiley (1974), chapters 2–4.

Chapter 2 is characteristics; chapter 3 is the breaking of a wave (the
`t_s = -1/min u₀'` result); chapter 4 covers the viscous smoothing of the
discontinuity and the tanh profile.

**What we take from it:**
* the breaking time (derived in the report in four lines, not quoted);
* the exact travelling/steady shock
  `u = -U tanh[U(x-x₀)/2ν]` with thickness `δ = 4ν/U` — this is the formula
  Experiment 2 checks against;
* the statement that after breaking the solution is a sawtooth, which is the
  starting point for the Fourier analysis.

**This is the paper being "replicated":** figs. 2(c) and 2(d) are the tanh
profile and the `4ν/U` scaling, computed from scratch and compared with the
formula.

---

### 3. J. Bec and K. Khanin, *Burgers turbulence*, Physics Reports **447**, 1–66 (2007).

A long review; only the first ~15 pages are used. Section 2 sets up the
problem and defines the energy spectrum; the early part of section 3 states
that for a randomly chosen smooth initial condition the spectrum develops a
`k⁻²` range once shocks have formed.

**What we take from it:** the spectral prediction `E(k) ∝ k⁻²` for a
shock-dominated solution, and the language used to describe it (inertial range,
dissipative cut-off).
**What we do not take:** the rigorous multifractal / structure-function theory
(`ζ_p = min(p,1)`, Lagrangian limits, shock statistics). That was the previous
version of this project and it is genuinely graduate-level; it is deliberately
absent here.

---

## Optional fourth, if a numerical-methods reference is wanted

R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge
University Press (2002), chapters 4 and 8 — for the von Neumann stability
condition `dt < dx²/2ν` and for the notion of numerical diffusion. Not needed
for the talk; useful if sir asks "where does your stability limit come from?".

---

## How to cite them in the report

* Introduction, for the equation: **[1]**.
* For `t_s`, the tanh profile and `δ = 4ν/U`: **[2]**, chapters 2–4.
* For `E(k) ∝ k⁻²`: **[3]**, section 2 and the start of section 3.

Three citations in three places, each doing real work.
