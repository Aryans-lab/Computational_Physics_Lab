# Report outline — what each paragraph must do

**Target:** exactly four pages, six figures, six references.
**Rule from the course:** the report must be in your own words. `main.tex` is a
working scaffold with the equations, figure placement and *correct numbers* —
rewrite every sentence, and check that you can defend each claim aloud before
you submit it.

---

## Page budget

| Section | Space | What it is for |
|---|---|---|
| Abstract | 6 lines | the question, the method, the three numbers |
| 1 Introduction | ½ page | why this equation, what $t_s$ and $\delta$ are |
| 2 Methods | ¾ page | five solvers, the exact reference, three diagnostics |
| 3 Results | 2 pages | five subsections, one figure each |
| 4 Discussion | ¼ page | what generalises, what does not |
| 5 Conclusion + refs | ¼ page | one paragraph and six lines |

If it overflows: drop Fig. 3 (`fig03`) — its content is a single sentence —
and shorten §3.5. Do **not** drop Figs. 1, 4 or 5; they are the argument.

---

## Paragraph by paragraph

### Abstract

One sentence of setup (`in 1-D the incompressibility constraint kills advection
and pressure, leaving Burgers`), one sentence of question (`when a shock is a
few grid points wide, what is physics and what is the mesh?`), one sentence of
method (five solvers, Cole–Hopf reference, three analytic predictions), one
sentence of the three headline numbers, one sentence of conclusion.

### 1. Introduction

- Paragraph 1: the equation; where it comes from physically (pressureless /
  compressible 1-D flow); state that the report is about numerics, not fluid
  dynamics.
- Paragraph 2: derive $t_s=-1/\min u_0'$ **in three lines** via characteristics.
  Say what happens after $t_s$ and why viscosity replaces it by a $\tanh$ front
  of thickness $4\nu/U$. *These two formulas are the entire theory the report
  needs.*
- Paragraph 3: the plan, in one sentence per experiment.

### 2. Methods

- Bullet-list the five schemes with their formal order. One line each. Do not
  describe the algorithms — the course already did.
- The Cole–Hopf reference deserves a short paragraph, including why the
  interpolation must be band-limited (this is a genuine methodological point
  and a likely viva question).
- Diagnostics: define $\delta$, $E(k)$, $S_p(r)$ compactly.

### 3. Results — one paragraph per figure, claim first

Each subsection opens with the **number**, then explains the figure:

1. **Verification.** "0.94, 2.02, 2.00 and exponential." Then the stability
   agreement (0.3%) and the cost ratio (1639 : 80 : 40).
2. **Breaking.** "$t_s=1$; strain $1.249$ vs $1.250$; saturation $247.5$ vs
   $250$; $\delta\propto\nu^{0.985}$."
3. **Spectra.** "Slope $-2.030\pm0.001$; cut-off $\propto\nu^{-0.971}$." The
   sentence that matters: viscosity does not change the slope, only the range.
4. **Burgulence.** "$\zeta_p = 0.52, 0.92, 1.00, 0.94, 0.92$ versus
   $\min(p,1)$ and $p/3$; the two predictions are a factor 3 apart at $p=4$."
   Be honest about the 8% deficit at $p=1$ and why.
5. **Grid viscosity.** "$\nu_{\rm num}=\frac12U\Delta x(1-C)$, verified to 1%
   over five decades; at $N=256$ the grid is 3.7 times more viscous than the
   fluid." End with the criterion $\Delta x\lesssim0.2\delta$.

### 4. Discussion

Exactly three paragraphs: (i) convergence order does not predict shock quality;
(ii) integral quantities are robust, local quantities are not; (iii) numerical
dissipation should be computed, not hoped away.
Then two limitations: the 1.5-decade inertial range, and the decaying
(non-stationary) statistics.

### 5. Conclusion

One paragraph, no new claims. Reuse the phrase "which features of a numerical
solution are physics and which are the mesh" — the report should answer the
question it asked.

---

## The numbers, in one place

| Quantity | Theory | Measured |
|---|---|---|
| upwind order | 1 | 0.94 |
| Lax–Wendroff order | 2 | 2.02 |
| Crank–Nicolson order | 2 | 2.00 |
| spectral error, $N=16\to64$ | exponential | $7.5\times10^{-4}\to10^{-9}$ |
| growth rate, FTCS $0.9\Delta t_{\rm crit}$ | $-82.3$ | $-82.1$ |
| growth rate, FTCS $1.1\Delta t_{\rm crit}$ | $+55.0$ | $+55.1$ |
| growth rate, CN at $20\Delta t_{\rm crit}$ | $-1.66$ | $-1.66$ |
| steps at $N=256$ | — | 1639 / 80 / 40 |
| breaking time | 1 | 1.00 |
| strain at $t=0.2$, $0.5$ | 1.250, 2.000 | 1.249, 1.993 |
| saturated strain | 250 | 247.5 |
| $\delta(\nu)$ | $4\nu/U$ | $\propto\nu^{0.985}$, prefactor 1.05–1.10 |
| spectrum | $k^{-2}$ | $k^{-2.030\pm0.001}$ |
| cut-off | $U/4\nu$ | $\propto\nu^{-0.971}$ |
| $\zeta_{0.5}, \zeta_1, \zeta_2, \zeta_3, \zeta_4$ | 0.5, 1, 1, 1, 1 | 0.52, 0.92, 1.00, 0.94, 0.92 |
| $\nu_{\rm num}$ | $\frac12U\Delta x(1-C)$ | ratio 1.00, slope 1.003 |
| $\nu_{\rm num}/\nu$ at $N=256$ | — | 3.7 |

Every one of these is reproduced by `python run_all.py` and stored in
`results/*.json`.

---

## Ten questions you will be asked

1. **"Why is this the 1-D Navier–Stokes equation?"**
   1-D incompressibility forces $\partial_xu=0$, so advection and pressure both
   vanish; Burgers is the surviving reduction once a pressure gradient is
   dropped (pressureless/compressible 1-D flow). Cite Burgers (1948).

2. **"How do you know the exact solution is exact?"**
   Two independent checks: it reproduces the initial condition as $t\to0$, and
   its finite-difference residual under the PDE is $4\times10^{-9}$
   (`tests/test_core.py::test_cole_hopf_satisfies_pde`). Also it is
   grid-independent to $10^{-12}$ when the quadrature grid is refined.

3. **"Why does the spectral scheme converge exponentially and then stop?"**
   Space: the solution is analytic between shocks, so Fourier coefficients
   decay faster than any power. Time: Strang splitting is second order, so
   below about $10^{-9}$ the time-stepping error dominates — that is the floor
   in Fig. 1a.

4. **"Why not use an adaptive step size?"**
   The diffusion term is stiff ($\nu k_{\max}^2$), so an explicit stability
   limit, not an accuracy requirement, sets $\Delta t$. Strang splitting and
   Crank–Nicolson remove that limit directly; adaptivity would fix the wrong
   problem.

5. **"Is your shock width measurement biased?"**
   Yes, by 5–10%, and that is stated: $\delta$ is defined as jump over maximum
   slope, which for a discrete profile underestimates the slope when the shock
   spans few points. The *scaling* with $\nu$ is unaffected.

6. **"Why is $\zeta_1$ 0.92 instead of 1?"**
   Finite box: at $r$ comparable to the integral scale, $S_1$ saturates and the
   local slope falls. Compare the local-slope table in
   `results/e4_structure_functions.txt`.

7. **"Could the $\zeta_p$ agreement be a coincidence?"**
   No: at $p=4$ the two predictions differ by 0.33 and the error bar is 0.05.

8. **"What is aliasing and did you deal with it?"**
   The product $u^2$ has bandwidth $2k_{\max}$; modes above $k_{\max}$ fold
   back. The 2/3 rule zeroes the top third before forming the product. Without
   it the high-$k$ spectrum is contaminated — check `dealias=False`.

9. **"Why is the numerical viscosity formula $\frac12U\Delta x(1-C)$?"**
   Truncation error of the upwind difference: the leading correction term is a
   second derivative with exactly that coefficient (modified-equation analysis,
   LeVeque ch. 8).

10. **"What would you do next?"**
    Force the equation to reach a statistically stationary state; add a
    WENO/TVD scheme and repeat the $\zeta_p$ measurement; extend to the 2-D
    case where vortex stretching (absent in 1-D) changes the physics.

---

## Submission checklist

- [ ] Report rewritten in your own words, 4 pages, all six figures referenced.
- [ ] `python run_all.py` finishes cleanly on a fresh checkout.
- [ ] `python tests/test_core.py` prints 15/15.
- [ ] Figures are the PDFs from `figures/`, not screenshots.
- [ ] Three references cited in the text, not only listed.
- [ ] Code headers carry name, roll number, problem statement.
- [ ] Zip the folder (or link the GitHub path) and upload before the deadline.
