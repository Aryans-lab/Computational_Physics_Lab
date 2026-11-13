# DIY Project — Brainstorm, Literature Map & Selection

**Course:** PHY341 / PHY745 — Physics Computer Lab (DIY = 20 marks)
**Student:** Aryan Bandyopadhyay (2411014), NISER
**Deadline shape:** codes + ≤4-page report before Endsem; 15 + 5 min presentation.

> This file is the *decision document*, not the project. It fixes the topic, the
> quantitative questions, the figure list, the paper trail, and the build plan.

---

## 0. What the brief actually rewards (read off the course material)

From `Slides/intro.pdf`, `COURSE_CONTEXT.md` and the master notes:

| Constraint | Consequence for topic choice |
|---|---|
| "Focus on **coding, plotting and interpretation of data** in the context of the problem" | The deliverable is *quantitative diagnostics*: convergence plots, error bars, scalings. A topic that produces only a pretty picture is a weak topic. |
| "**Not expected to discuss the physics** of the problem" | Deep theory is a liability, not an asset. Choose a system whose equations are 2 lines long but whose *numerics* are rich. |
| "You can use **more than one** idea from the list" | The best projects sit at an intersection — it reads as ambition without extra risk. |
| "Find physics problems where these ideas are **applied**" | Every method must be attached to a physical system, not a textbook function. |
| Report ≤ 4 pages, **own words** | Small number of sharp results, not breadth. 3 results + 1 artifact is ideal. |
| Marks: **DIY 20**, code audited like assignments | Code quality and reproducibility matter as much as the figures. |
| Course ethos: **everything from scratch** (`mylib.py`) | Reusing/extending `mylib.py` is a *feature*. It also means "no black boxes" is part of the story we can tell in the report. |

**What "publishable" realistically means here.** Not a research paper. The honest
ladder is:

1. **Zenodo DOI on the GitHub repo** (guaranteed, 10 minutes, citable, looks great on a CV).
2. **arXiv preprint** in `physics.comp-ph` / `nlin.CD` — needs an endorser, but a
   tidy 8–10 page numerical-experiment note is exactly what such preprints are.
3. **AJP / European Journal of Physics / American Journal of Physics "computational
   physics" note** — refereed, pedagogy-oriented, and *very* receptive to "here is a
   numerical artifact students should know about" papers.

Any topic below that scores high on "artifact" rather than "simulation" has a
realistic shot at (2) and (3). Pure reproductions (SAW exponents, Lorenz attractor)
do not.

---

## 1. Scoring criteria

Six axes, 1–5. Weights reflect *your* stated priorities: finishable, good on GitHub,
interesting, plausible preprint, useful for the summer-2027 internship story.

- **Feasible** — can it be finished and debugged in the remaining semester?
- **Visual** — quality/quantity of publication-grade figures.
- **GitHub** — does the repo look like a real project (library + CLI + tests + docs)?
- **Novel** — is there a result a referee/recruiter has not seen 500 times?
- **Course fit** — how directly it satisfies the DIY rubric (diagnostics, interpretation).
- **CV synergy** — does it reinforce the N-body / instrumentation / detector narrative?

---

## 2. The keyword list → candidate physics problems

| Keyword | Candidate physics problem | Numerics you would actually build | Verdict |
|---|---|---|---|
| **Self-avoiding random walk** | Polymer scaling; SAW on 2D/3D lattices vs. simple random walk | Pivot algorithm, autocorrelation/blocking, log–log fits with χ², Flory theory comparison | Easy, pretty, *completely reproduced* since 1988 |
| **Adaptive Runge–Kutta** | Eccentric orbits, close encounters, stiff kinetics, chaotic flows | Embedded RK pairs (RKF45 / Dormand–Prince 5(4)), PI step control, Richardson cross-check | Strong — but needs a *system* to be about |
| **Singular Value Decomposition** | Abel inversion of line-integrated emission; low-rank denoising of 2D fields; PCA of spectra | One-sided Jacobi SVD, truncated SVD + Tikhonov, L-curve, condition-number analysis | Excellent "experimental physics" flavour; less "physics" per se |
| **Chaos** | Lorenz-63 convection, Rössler, driven pendulum, three-body scattering | Lyapunov spectra (Benettin/QR), Poincaré sections, bifurcation sweeps, invariant measures | Saturated on GitHub unless paired with a real question |
| **Fourier analysis** | Energy spectra of a turbulent field; lock-in recovery of a buried signal; spectral leakage | FFT from scratch (radix-2), windowing, Welch PSD, non-uniform sampling | Good support act; weak as a solo act for 20 marks |
| **1-dim Navier–Stokes** | Burgers equation: shock formation, Burgers turbulence, dissipation range | Fourier pseudospectral + ETDRK4, dealiasing, upwind vs spectral comparison, structure functions | Most *physics* per line of code; hardest to finish cleanly |

The three intersections worth taking seriously:

- **A. Adaptive RK × Chaos** — "which numbers from a chaotic simulation can you trust?"
- **B. Adaptive RK × Chaos × gravitation** — same question in the three-body problem.
- **C. 1D Navier–Stokes × Fourier** — "Burgulence": shocks, spectra, and what numerical
  dissipation fakes.

---

## 3. Finalist A (RECOMMENDED) — *Trusting the Butterfly*

### Adaptive step control, Lyapunov horizons, and which numbers from a chaotic simulation you can actually believe

**Physics system.** Lorenz-63, the three-mode truncation of Rayleigh–Bénard
convection (σ = 10, ρ = 28, β = 8/3). Three equations, no theory needed, and it is
*the* canonical chaotic flow with a published λ₁ ≈ 0.906 — a number we can
reproduce, which doubles as a validation of the whole pipeline.

**Why it is not "just another Lorenz plot".** Everyone plots the butterfly. Nobody in
a lab course asks the question that matters: *the trajectory is garbage after ~30
Lyapunov times — so why do published statistical results from chaotic simulations
still work?* That answer is **shadowing**, it is a real research topic (Hammel–Yorke–
Grebogi; Sauer–Grebogi–Yorke), and it can be measured with a laptop.

### Three experiments (each = one result, one figure, one paragraph of interpretation)

**E1 — Build the integrator and prove it is correct.**
Dormand–Prince 5(4) embedded pair (7 stages, FSAL), local error from the 5th/4th
difference, PI controller `h ← h·(tol/err)^{1/5}` with safety 0.9, rejection counter,
`h_min/h_max` guards. Cross-check the error estimate by Richardson step-doubling.
Validate order on an *analytic* problem: eccentric Kepler orbit (e = 0.9), where the
global error must fall as h⁵ and where fixed-step RK4 wastes ~10–50× the function
evaluations near pericentre. **Figure:** error vs. cost, adaptive vs. fixed-step, plus
a histogram of accepted step sizes through one orbit.

**E2 — Measure the horizon.**
Integrate the same IC at tolerance τ and at a reference τ_ref = 10⁻¹³; track
d(t) = ‖y_τ − y_ref‖. Prediction: d(t) ≈ δ₀ e^{λ₁ t}, so a plot of ln d vs t has
slope λ₁ = 0.906 and the *forecast horizon* obeys

    T_valid(τ) ≈ (1/λ₁) ln(Δ_valid/δ₀).

**This is the money plot:** T_valid λ₁ collapses onto a single straight line in
log₁₀ τ over 9 decades of tolerance, with slope ≈ 1/λ₁. You are measuring a
dynamical invariant by watching your own numerical error grow.

**E3 — Statistical fidelity survives trajectory death.**
For τ ∈ [10⁻², 10⁻¹²], run N = 200-member ensembles for T = 500 (≈ 450 Lyapunov
times) and compute four observables:

- λ₁ itself (Benettin two-trajectory + periodic renormalisation),
- ⟨x⟩, Var(z) (invariant measure moments),
- the power spectral density of x(t), and its low-frequency slope,
- the mean return time to a Poincaré section z = 27.

**Punchline:** the trajectory is decorrelated from the reference after T ≈ 20–40, but
the moments converge to <1 % at τ ≈ 10⁻⁶, i.e. **six to eight orders of magnitude
sooner**. One log–log figure, two families of curves, and the sentence that carries
the report: *chaotic simulations are not predictions, they are samplers — and they are
good samplers long before they are good predictors.*

### E4 — The artifact (this is what makes it preprint-grade)

An adaptive solver accepts **non-uniform** time steps. Two consequences that bite
real published work, and that we can demonstrate:

1. **Sampling bias.** Histogramming state at accepted steps weights each interval
   equally, not each interval *by its duration*. Regions where h is large are
   systematically under-sampled: bias O(⟨h⟩ · d ln h/dt). Show the measured bias in
   ⟨x⟩ and in the PSD against the uniform-grid (correct) answer, and fix it by
   weighting with h or interpolating to a uniform grid. **This effect is invisible to
   the E2 test — the trajectory is fine, the statistics are subtly wrong.**
2. **Spectral distortion.** An FFT of an adaptively-stepped series with non-uniform
   sample times produces spurious high-frequency power. Show the naive FFT vs the
   resampled FFT against the PSD from a fixed-step high-resolution reference.

(E4 is the "no one told me this" result. It is exactly the kind of thing a referee at
AJP accepts, because it is *practically* useful and *reproducibly* demonstrated.)

**Bonus (optional, 1 figure):** the same machinery on the gravitational three-body
problem, where a close encounter forces h → 0 and the adaptive controller either saves
you or silently destroys energy conservation. Ties the project to your N-body work.

### What it hits from the keyword list
Adaptive Runge–Kutta (core) + Chaos (core) + Fourier analysis (E4 spectra + E3 PSD).
Three of six, with the first two genuinely central.

### References (core three, all verified)
1. E. N. Lorenz, *Deterministic nonperiodic flow*, J. Atmos. Sci. **20**, 130 (1963).
2. J. R. Dormand & P. J. Prince, *A family of embedded Runge–Kutta formulae*,
   J. Comput. Appl. Math. **6**, 19 (1980).
3. S. M. Hammel, J. A. Yorke & C. Grebogi, *Do numerical orbits of chaotic dynamical
   processes represent true orbits?*, J. Complexity **3**, 136 (1987).

Supporting (pick 2): T. Sauer, C. Grebogi & J. A. Yorke, *How long do numerical
chaotic solutions remain valid?*, Phys. Rev. Lett. **79**, 59 (1997); E. N. Lorenz,
*Computational periodicity as observed in a simple system*, Tellus A **58**, 549 (2006);
E. Hairer, S. P. Nørsett & G. Wanner, *Solving Ordinary Differential Equations I*
(Springer, 1993), Ch. II.

### Risks
- Convergence of ensemble statistics needs care (block averaging, not naive σ/√N —
  autocorrelated samples!). **Mitigation:** this is itself a nice methods paragraph,
  and the course already taught χ²/error scaling.
- Python speed: N = 200 × 10⁶ steps is a lot. **Mitigation:** vectorise the ensemble
  in NumPy (200 trajectories advance in lockstep with a *shared* step size is wrong —
  instead use per-trajectory step sizes with masked updates, or run 200 independent
  integrations at 10⁶ steps each, ~minutes each). Feasible; measure it early.
- E4 must be demonstrated cleanly or it looks like a bug. **Mitigation:** prove it on
  a system with an exact invariant measure (e.g. harmonic oscillator) before Lorenz.

**Effort:** 2–3 focused weeks. **Risk of not finishing:** low.

---

## 4. Finalist B — *Chaos in the gravitational three-body problem*

**Question.** After how many Lyapunov times does a three-body simulation stop being
reversible, and do the *outcome statistics* (escape speed distribution, binary
semimajor-axis distribution, collision fraction) survive when individual orbits do not?

**Numerics.** Adaptive RK45 with a close-encounter-aware controller; optional
regularisation (KS or the "logarithmic Hamiltonian" trick); fixed-step leapfrog and a
symplectic integrator for comparison; exponential-divergence measurement; ensembles of
10⁴–10⁵ random initial triples.

**Figures.** Divergence of forward/backward integrations; energy error vs. tolerance;
outcome-distribution stability across 4 decades of tolerance; the classic "which body
escapes" sensitivity map.

**Why it is attractive.** It is *literally* the subject of published papers
(Boekholt & Portegies Zwart), it plugs directly into your N-body project with
Dr. Khandai, and the pictures are spectacular. A student note here has a plausible
arXiv home in `astro-ph.IM` / `physics.comp-ph`.

**Why I rank it second.** Close encounters + regularisation + 10⁴ orbits is a
4–6 week project with several failure modes (energy drift, step-size pathology at
pericentre, performance). The statistical machinery is the same as A but with a
harder dynamical engine. **Best version of B = A with the three-body as the
application instead of Lorenz.**

**References.** T. C. N. Boekholt & S. F. Portegies Zwart, *On the reliability of
N-body simulations*, Comput. Astrophys. Cosmol. **2**, 2 (2015); T. C. N. Boekholt
et al., MNRAS **452**, 3941 (2015); Boekholt & Portegies Zwart, MNRAS **476**, 496 (2018).

---

## 5. Finalist C — *Burgulence: 1D Navier–Stokes, shocks and spectra*

**Question.** How much of the "turbulence" you see in a 1-D Navier–Stokes
(Burgers) simulation is physics and how much is your grid?

**Numerics.** u_t + u u_x = ν u_xx (+ stochastic forcing). Fourier pseudospectral
solver with 2/3-rule dealiasing and ETDRK4 exponential time differencing; comparison
against a first-order upwind finite-difference solver at the same nominal resolution.

**Diagnostics.** Shock formation time from smooth IC (t_s = 1 for u₀ = sin x);
shock thickness δ ≈ ν/U; energy spectrum E(k) ~ k⁻² for a shock train; structure
functions S_p(r) ~ r^{ζ_p} with the bifractal prediction ζ_p = min(p, 1); PDF of
velocity gradients with its |s|^{-3/2} tail.

**The punchline.** The measured ζ_p *changes* when the grid is too coarse — numerical
dissipation mimics viscosity, and the "measured" scaling exponent is a property of your
mesh, not of the fluid. A resolution-convergence study of an exponent is a genuinely
nice piece of data interpretation.

**Assessment.** The most *physical* of the three, the most beautiful figures
(space–time heat maps of coalescing shocks), and the hardest to get right: spectral
methods, stiff time stepping, long runs. It also has zero overlap with anything already
in your repo, which cuts both ways. **Effort:** 4 weeks, medium–high risk.

**References.** J. Bec & K. Khanin, *Burgers turbulence*, Phys. Rep. **447**, 1
(2007); J. M. Burgers, Adv. Appl. Math. Mech. **1**, 171 (1948); Basdevant et al.,
*Spectral and finite difference solutions of the Burgers equation*, Comput. Fluids
**14**, 23 (1986) *(verify volume before citing)*.

---

## 6. The rest of the list (why they are not the pick)

**Self-avoiding random walk.** Build the pivot algorithm, measure ν = 0.5876 (3D) and
3/4 (2D exact), do the Flory comparison, block-average the error bars. It is a *lovely*
2-week project and the clearest possible demonstration of "Monte Carlo + scaling +
fitting", and it reuses your LCG. But: the answer is known to 7 digits (Clisby 2010),
the pictures are walks, and there is no original question left to ask at this level.
Good safety net, weak CV artifact.

**Singular value decomposition.** The best physics application is **Abel inversion**:
you measure a line-integrated emission profile (plasma column, flame, discharge — very
much your experimental world), and must invert an ill-conditioned triangular operator
to get radial emissivity. Truncated SVD vs Tikhonov, L-curve, noise amplification
~1/σ_min. Add PCA of a 2-D field video as a second act. Excellent "experimental data
analysis" signal, distinctly different from your current repo, but it is signal
processing more than physics, and the "discovery" is small.

**Fourier analysis.** Strongest solo version is a **lock-in amplifier from scratch**:
bury a 1 µV sinusoid in 1/f + white noise, recover it by phase-sensitive detection,
and measure SNR gain vs. integration time and equivalent noise bandwidth; then beat it
with an FFT + Welch PSD, and note when each wins. It is *directly* connected to your
ESR/lock-in lab work and would read beautifully in an internship interview — but as a
20-mark DIY it is one method and one system.

**1-dim Navier–Stokes (literal reading: channel flow).** Plane Poiseuille / Womersley
flow, Crank–Nicolson in time, startup transients. Correct, quick, and dull; no
interesting question survives.

**Rejected outright:** double pendulum / Rössler galleries / Mandelbrot / "Lorenz
attractor in 3-D". Saturated, and they are plots, not experiments.

---

## 7. Scorecard

| Project | Feasible | Visual | GitHub | Novel | Course fit | CV synergy | **Total /30** |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **A. Adaptive RK × Chaos (Trusting the Butterfly)** | 5 | 4 | 5 | 4 | 5 | 4 | **27** |
| **B. Three-body reliability** | 3 | 5 | 5 | 5 | 4 | 5 | **27** |
| **C. Burgulence (1D NS × Fourier)** | 3 | 5 | 5 | 4 | 4 | 2 | **23** |
| SVD / Abel inversion | 4 | 4 | 4 | 3 | 3 | 5 | **23** |
| SAW critical exponents | 5 | 4 | 4 | 2 | 4 | 2 | **21** |
| Lock-in + FFT | 5 | 3 | 4 | 3 | 3 | 5 | **23** |

A and B tie on points; **A wins on risk** (identical intellectual content, much lower
chance of a half-finished project in week 12) and B wins on ceiling. The clean
compromise: **build A, and use the three-body problem as A's final application
section** — same code, one extra figure, both stories told.

---

## 8. Recommendation

> **Take A, with the three-body as the closing application.**

One-sentence pitch for the report abstract:

> *Chaotic simulations are routinely run to thousands of Lyapunov times, long after
> individual trajectories have decorrelated from any true orbit. Using a
> from-scratch adaptive Runge–Kutta 5(4) integrator, we quantify the gap between
> trajectory fidelity and statistical fidelity in the Lorenz system: moments and
> spectra converge six to eight orders of magnitude earlier in tolerance than
> trajectories do, while naive sampling of adaptively-stepped output introduces a
> tolerance-dependent bias in those same statistics.*

Why this one, specifically for you:

1. **It cannot fail.** If the ensembles are slow, shrink N; if E4 misbehaves, E2 and E3
   already carry a 4-page report.
2. **It is a *question*, not a demonstration.** Every recruiter has seen a Lorenz plot;
   nobody has seen "here is the tolerance at which my statistics stop moving".
3. **It builds on what you already have** — `mylib.rk4` becomes `mylib.rk45`; the
   Rössler assignment becomes a validation case; the Fourier material becomes E4.
   Your repo gains a coherent arc instead of a folder of homework.
4. **It has a natural internship sentence**: "I built an adaptive integrator and
   quantified when numerical error destroys trajectory vs. statistical predictions, and
   found a sampling bias that silently corrupts spectra." That is a sentence an
   experimental group building detectors and a theory group running simulations both
   understand.
5. **It is the honest path to a preprint.** "Practical artifact + measured convergence
   law" is publishable pedagogy; "I reproduced the Lorenz attractor" is not.

---

## 9. If you pick A — build plan

### Repo layout (consistent with the existing conventions)

```
diy_project/
├── README.md                 # the project front page: question, 4 figures, results table
├── BRAINSTORM.md             # this file
├── LITERATURE.md             # annotated bibliography (10 papers, 3-4 sentences each)
├── refs.bib                  # BibTeX, DOI-verified
├── src/
│   ├── integrators.py        # rk45_dopri (FSAL, PI control), rk4 fixed, richardson_check
│   ├── systems.py            # lorenz63, rossler, kepler, three_body (r^-2 softened)
│   ├── diagnostics.py        # lyapunov_benettin, psd_welch, invariant moments, block averaging
│   └── experiments/
│       ├── e1_integrator_validation.py
│       ├── e2_lyapunov_horizon.py
│       ├── e3_statistical_fidelity.py
│       └── e4_sampling_bias.py
├── figures/                  # publication-grade PDF + PNG, generated, never hand-edited
├── data/                     # small CSV/JSON of converged observables (committed)
├── results/                  # summary tables (LaTeX + CSV)
└── tests/                    # analytic-solution tests for the integrator
```

### Figure list (target 5 figures for the report, 8 for slides)

| # | Figure | Type | Purpose |
|---|---|---|---|
| F1 | Kepler e = 0.9: global error vs. f-evaluations, adaptive vs fixed RK4 | log–log | proves the integrator is 5th order and cheaper |
| F2 | Lorenz: accepted step size vs. time, coloured by local error | time series | shows the controller "feels" the dynamics |
| F3 | ln d(t) vs t for 5 tolerances, with fitted slopes = λ₁ | semi-log | the measured horizon law |
| F4 | Trajectory error **and** 4 statistical errors vs tolerance | log–log, one panel | **the punchline** |
| F5 | Sampling-bias artifact: biased vs corrected ⟨x⟩ and PSD | 2 panels | the artifact |
| F6 | Three-body: energy error vs tolerance; divergence vs. Lyapunov time | 2 panels | application / outreach figure |
| F7 | Poincaré section + invariant density | scatter/heatmap | visual anchor for "sampler not predictor" |
| F8 | Validation: PSD of Lorenz x(t) vs published broadband shape | spectrum | Fourier cross-check |

### Six-week schedule (part-time, alongside coursework)

| Week | Deliverable |
|---|---|
| 1 | `integrators.py` + `systems.py`; E1 passes (order 5 on Kepler, and on y' = y with exact solution) |
| 2 | Lyapunov machinery (Benettin + QR); reproduce λ₁ = 0.906 ± 0.005; E2 horizon plot |
| 3 | Ensemble driver + block-averaged error bars; E3 running overnight |
| 4 | E4 sampling bias (prove on harmonic oscillator first, then Lorenz); Fourier pipeline |
| 5 | Three-body section, README, `refs.bib`, Zenodo DOI, figures finalised |
| 6 | Slides + talk rehearsal; report skeleton handed over for your own writing |

### 4-page report skeleton (you write the prose)

1. **Motivation** (½ page) — the reliability question; what shadowing does and does not
   guarantee.
2. **Method** (1 page) — DOPRI5, error control, Benettin exponents, ensemble/block
   averaging, PSD. Equations only, no derivations.
3. **Results** (1¾ pages) — F3, F4, F5 with the numbers.
4. **Limitations & outlook** (½ page) — double precision floor, stiff vs. non-stiff,
   the three-body case.
5. **References** (¼ page) — 5–6 entries.

### 15 + 5 minute talk arc

1. The butterfly (30 s, no explanation needed) → "this trajectory is wrong after t ≈ 30".
2. Question: then why do we trust the statistics? (1 min)
3. The integrator, and proof it is right (2 min)
4. Measuring the horizon (3 min)
5. Statistic vs. trajectory: the one figure (4 min)
6. The artifact: your spectrometer/PSD is lying to you (3 min)
7. Three-body coda (1.5 min)
8. What I would do next / questions (1 min) — leave 5 min for the Q&A the professor
   will spend on *"how do you know the reference is right?"* (answer: tolerance
   stacking + Richardson + published λ₁).

---

## 10. Literature to read before week 1 (in this order)

1. Lorenz 1963 — skim §1–2 only, the equations and the "nonperiodic" claim. 20 min.
2. Dormand & Prince 1980 — the table of coefficients and the error-control paragraph.
   30 min. (Or Hairer–Nørsett–Wanner Ch. II.4 for the readable version.)
3. Hammel, Yorke & Grebogi 1987 — the shadowing argument, the refinement procedure. 1 h.
4. Sauer, Grebogi & Yorke 1997 — shadowing time vs. Lyapunov time; the scaling. 45 min.
5. Lorenz 2006 (Tellus A) — computational periodicity: what goes wrong with big steps.
   30 min.
6. Boekholt & Portegies Zwart 2015 — how the same logic plays out in N-body. 1 h.
7. (If C is ever reconsidered) Bec & Khanin 2007 — review, read §1 and the
   structure-function section only.

---

*All DOIs to be collected in `refs.bib` and checked at build time. Citations marked
"verify" above must not be copied into the report until the volume/page is confirmed.*
