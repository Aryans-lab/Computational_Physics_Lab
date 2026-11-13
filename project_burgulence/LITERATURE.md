# Literature review — annotated bibliography

**Project:** shock waves in a one-dimensional fluid (viscous Burgers equation)
**Student:** Aryan Bandyopadhyay (2411014)

Read in this order. The three marked ★ are the ones to cite in the four-page
report; the rest are the background that makes the three defensible in a viva.

---

## 1. What the equation is, and where it comes from

### ★ J. M. Burgers, *A mathematical model illustrating the theory of turbulence*,
Adv. Appl. Mech. **1**, 171–199 (1948).

The original proposal: keep the nonlinear advection of the Navier–Stokes
equation, throw away the pressure gradient, and see what survives. What
survives is a one-dimensional model in which the energy cascade of real
turbulence becomes a cascade of *shocks* — steepening, colliding, merging.
**What to take:** the motivation for calling this the 1-D Navier–Stokes
equation, and the physical picture of a "shock train".

### G. B. Whitham, *Linear and Nonlinear Waves* (Wiley, 1974), ch. 2–4.

Chapters 2 and 3 give the method of characteristics and the exact breaking
time for a simple wave; chapter 4 gives the viscous shock profile as a
travelling wave of the Burgers equation, which is exactly the
$u=-U\tanh[U(x-x_0)/2\nu]$ used as an analytic benchmark in experiment 2.
**What to take:** the derivation of $t_s=-1/\min u_0'$ and of $\delta=4\nu/U$.

---

## 2. The exact solution

### J. D. Cole, *On a quasi-linear parabolic equation occurring in aerodynamics*,
Quart. Appl. Math. **9**, 225–236 (1951).

### E. Hopf, *The partial differential equation $u_t+uu_x=\mu u_{xx}$*,
Comm. Pure Appl. Math. **3**, 201–230 (1950).

Independently discovered: the substitution $u=-2\nu\,\partial_x\ln\phi$ maps the
nonlinear Burgers equation onto the linear heat equation. This is the project's
reference solution. **What to take:** the transformation and the periodic heat
kernel used in `src/exact.py`.

> Practical lesson learned the hard way: the Cole–Hopf integral must be
> evaluated with the initial data interpolated *band-limited* onto the
> quadrature grid. Linear interpolation of a coarse grid injects an
> $O(\Delta x^2)$ error into the benchmark, which then masquerades as the
> truncation error of the solver being tested.

---

## 3. Burgulence — the statistics of many shocks

### ★ J. Bec and K. Khanin, *Burgers turbulence*, Phys. Rep. **447**, 1–66 (2007).

The definitive review. Sections 1–3 are enough for this project: the shock
adiabatic invariant, the $E(k)\propto k^{-2}$ spectrum of a shock train, the
asymptotic shock density, and the bifractal structure functions
$\zeta_p=\min(p,1)$. **What to take:** every prediction that experiment 4 tests,
and the warning that these are *asymptotic* results — finite $\nu$ and a finite
box both bend them, which is exactly what the fit window has to handle.

### U. Frisch, *Turbulence* (Cambridge University Press, 1995), ch. 5–6.

Kolmogorov's 1941 theory and the $\zeta_p=p/3$ prediction, included here only as
the foil: the measured exponents in experiment 4 are a factor of three away from
it at $p=4$, which is what makes the measurement meaningful.

---

## 4. The numerical methods

### ★ P. A. Lagerstrom, *Matched Asymptotic Expansions* (Springer, 1988).

Inner/outer structure of the viscous shock: the outer solution is the
inviscid sawtooth, the inner solution is the tanh profile, and the shock
thickness is where the two meet. **What to take:** the language to explain why
a shock is a *singular perturbation* and why the grid has to resolve the inner
layer — the physical statement behind experiment 5.

### J. P. Boyd, *Chebyshev and Fourier Spectral Methods*, 2nd ed. (Dover, 2001).

Chapters 2 and 4 cover the Fourier pseudo-spectral method, aliasing and the
2/3 rule; chapter 10 covers integrating-factor and split-step time stepping.
**What to take:** the de-aliasing mask and the Strang splitting in
`src/solvers.py`.

### S. A. Orszag, *On the elimination of aliasing in finite-difference schemes by
filtering high-wavenumber components*, J. Atmos. Sci. **28**, 1074 (1971).

The origin of the 2/3 rule.

### G. Strang, *On the construction and comparison of difference schemes*,
SIAM J. Numer. Anal. **5**, 506–517 (1968).

Why half-step/full-step/half-step is second order. The spectral solver uses
this: exact diffusion in Fourier space, RK4 for the advection.

### R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*
(Cambridge University Press, 2002), ch. 4, 6, 8.

Upwind, Lax–Wendroff and the modified-equation analysis that produces
$\nu_{\rm num}=\tfrac12U\Delta x(1-C)$ — the formula that experiment 5 verifies
to 1%. Chapter 8 also explains why a centred scheme rings where an upwind
scheme smears, which is the substance of `fig01`–`fig02`.

### J. von Neumann and R. D. Richtmyer, *A method for the numerical calculation of
hydrodynamic shocks*, J. Appl. Phys. **21**, 232–237 (1950).

The origin of the stability analysis used in experiment 1: the amplification
factor of the grid-scale mode. The measured growth rates agree with the
prediction to 0.3%.

---

## 5. Where the writing conventions come from

### N. J. Higham, *Handbook of Writing for the Mathematical Sciences*,
2nd ed. (SIAM, 1998).

For the structure of a short technical report: one claim per figure, figures
that can be read without the text, and every reported number traceable to a
run. The `results/` directory is that principle applied literally.

---

## The three to cite

If the report may only carry three references:

1. **Bec & Khanin (2007)** — the physics of Burgulence, and the predictions.
2. **Burgers (1948)** — the model itself.
3. **Whitham (1974) ch. 4** — the analytic shock profile and the breaking time.

Add LeVeque (2002) if a numerical-methods citation is wanted for the grid
viscosity formula, and Cole (1951) if the exact solution needs attribution
rather than a derivation.

---

## Reading plan (4 hours total)

| Time | What | Why now |
|---|---|---|
| 45 min | Whitham ch. 2–4 | you must be able to derive $t_s$ and the tanh profile on a whiteboard |
| 60 min | Bec & Khanin §1–3 | the review's figures are the standard the project's figures are judged against |
| 45 min | LeVeque ch. 4, 6 | modified-equation argument for $\nu_{\rm num}$ |
| 40 min | Boyd ch. 2 | why the spectral method wins, and what the 2/3 rule does |
| 30 min | Cole (1951) | the transformation, skim only |

*Every entry above should be checked against the DOI in `refs.bib` before it is
copied into the report.*
