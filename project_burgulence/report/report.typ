// ---------------------------------------------------------------------------
//  DIY Project report -- PHY341/745 Physics Computer Lab, NISER
//  Aryan Bandyopadhyay (2411014)
//
//  Typeset with Typst (compile: see ../make_pdfs.js, or `typst compile`).
//  Figures: ../figures/*.png, produced by ../run_all.py
//
//  NOTE: this is a scaffold with the correct numbers. Rewrite the prose in
//  your own words before submitting.
// ---------------------------------------------------------------------------

#set document(title: "Shock waves in a one-dimensional fluid", author: "Aryan Bandyopadhyay")
#set page(paper: "a4", margin: (x: 1.75cm, top: 1.5cm, bottom: 1.5cm),
  numbering: "1 / 4")
#set text(font: "STIX Two Text", size: 9.6pt, lang: "en")
#show math.equation: set text(font: "STIX Two Math", size: 10pt)
#set math.equation(numbering: "(1)")
#set par(justify: true, leading: 0.62em, spacing: 0.62em)
#set heading(numbering: "1.")
#show heading: it => block(above: 1.0em, below: 0.55em, {
  set text(size: 10.8pt, weight: 700, fill: rgb("#1f4e79"))
  it
})
#show heading.where(level: 2): it => block(above: 0.75em, below: 0.35em, {
  set text(size: 9.9pt, weight: 700, fill: rgb("#1f4e79"))
  it
})
#set figure.caption(separator: [ ])
#show figure.caption: set text(size: 8.3pt)
#show figure: set block(breakable: false, above: 0.55em, below: 0.55em)
#show ref: it => text(fill: rgb("#1f4e79"), it)

#align(center)[
  #text(size: 15pt, weight: 700, fill: rgb("#1f4e79"))[
    Shock waves in a one-dimensional fluid]
  #v(-0.4em)
  #text(size: 11pt)[How much of what you see is the grid? \ A numerical study of the viscous Burgers equation]
  #v(-0.2em)
  #text(size: 8.8pt)[Aryan Bandyopadhyay \ (2411014) \ — PHY341/745 Physics Computer Lab, NISER Bhubaneswar]
]
#v(-0.6em)
#line(length: 100%, stroke: 0.6pt + rgb("#1f4e79"))
#v(-0.2em)

#pad(left: 1.2em, right: 1.2em)[
  #text(size: 8.9pt)[
    *Abstract.* In one spatial dimension the incompressibility constraint removes
    both advection and the pressure gradient from the Navier--Stokes equation,
    leaving the viscous Burgers equation $u_t + u u_x = nu u_(x x)$ as the only
    non-trivial one-dimensional reduction. It contains two competing effects ---
    nonlinear steepening, which breaks a smooth wave in finite time, and
    viscosity, which prevents the discontinuity --- and therefore isolates, in the
    simplest possible setting, the question this report is about: when a shock is
    only a few grid points wide, which features of a numerical solution are
    physics and which are properties of the mesh? Using five independent solvers,
    an exact reference from the Cole--Hopf transformation, and three analytic
    predictions (shock thickness, $k^(-2)$ spectrum, bifractal scaling
    exponents), we measure where each scheme stops being trustworthy, and show
    that the numerical viscosity of a first-order scheme can exceed the physical
    viscosity by a factor of seven on a grid that would otherwise look perfectly
    adequate.
  ]]
#v(-0.3em)

= Introduction

The Burgers equation
$ u_t + u u_x = nu u_(x x), quad x in [0, 2 pi), quad u(0,t) = u(2 pi,t) $ <eq:burgers>
is the standard one-dimensional model of a compressible, dissipative flow. Each
point of the fluid travels at its own speed $u$, so a smooth profile steepens: a
gradient $u_x < 0$ grows as $1\/ (1 + t u_x)$ and reaches infinity at
$ t_s = -1 / op("min")_x u_0'(x). $ <eq:ts>
Viscosity prevents the singularity and replaces it by a stationary shock whose
thickness is set by the only length scale available,
$ delta = 4 nu \/ U, $ whose profile is the exact travelling wave
$u = -U tanh[U (x - x_0) \/ 2 nu]$.

Two facts make this an ideal test bed for numerical work. First, every quantity
of interest has an analytic value, so numerical error can never be mistaken for
physics. Second, the interesting structures --- shocks --- are precisely the
structures that discretisations handle worst. This report therefore does not try
to _explain_ the physics of shocks; it uses them as a ruler.

= Methods

== Five solvers

@eq:burgers is integrated with (i) first-order upwind, (ii)
Lax--Wendroff, (iii) forward-time centred-space (FTCS), (iv) Crank--Nicolson in
the diffusion with a Heun step in the advection, and (v) a Fourier
pseudo-spectral scheme with $2\/3$-rule de-aliasing and Strang splitting (exact
diffusion, fourth-order Runge--Kutta advection). The periodic tridiagonal systems
in (iv) are solved with the Thomas algorithm plus a Sherman--Morrison correction.
All routines were written for this project; only array arithmetic uses NumPy.

#align(center)[
  #table(
    columns: (auto, auto, auto, auto),
    align: (left, left, left, left),
    stroke: none,
    inset: (x: 6pt, y: 2.2pt),
    table.hline(stroke: 0.7pt),
    [*scheme*], [*order*], [*time step at $N = 256$*], [*comment*],
    table.hline(stroke: 0.4pt),
    [upwind], [1], [185], [diffusive],
    [Lax--Wendroff], [2], [93], [dispersive],
    [FTCS], [2], [—], [unconditionally unstable for diffusion],
    [Crank--Nicolson + Heun], [2], [40], [no diffusion limit],
    [spectral + Strang], [exp.], [80], [no advection limit],
    table.hline(stroke: 0.7pt),
  )
]

== Exact reference

The Cole--Hopf substitution $u = -2 nu partial_x ln phi$ linearises
@eq:burgers into the heat equation, giving
$ phi(x,t) = integral_0^(2 pi) K(x - y, t) thin phi_0(y) dif y, quad
  phi_0 = exp[-1/(2 nu) integral_0^y u_0] $ <eq:colehopf>
with $K$ the periodic heat kernel. The integral is evaluated on an 8192-point
grid, into which $u_0$ is interpolated _band-limited_: linear interpolation of a
coarse grid would inject an $O(Delta x^2)$ error into the benchmark itself. The
resulting reference is grid-independent to $10^(-12)$, and its residual under
@eq:burgers is $4 times 10^(-9)$.

== Diagnostics

Three measurements carry the argument: the shock thickness
$delta = ("jump in " u) \/ max abs(u_x)$, which equals $4 nu \/ U$ for the exact
profile; the energy spectrum $E(k)$, obtained from a radix-2 FFT written from
scratch and verified against a naive DFT to $10^(-15)$; and the structure
functions $S_p (r) = chevron.l abs(u(x+r) - u(x))^p chevron.r tilde.op r^(zeta_p)$.

= Results

== Verification: three orders, and a stability boundary

@fig:validation (a) compares each scheme with the exact solution at $T = 0.5$
(before breaking, $nu = 0.1$) with the time step held fixed, so that only the
spatial error is visible. The measured orders are $0.94$ (upwind), $2.02$
(Lax--Wendroff) and $2.00$ (Crank--Nicolson); the spectral scheme drops from
$7.5 times 10^(-4)$ to $10^(-9)$ between $N = 16$ and $N = 64$ and then flattens
at the Strang-splitting floor. Panel (b) checks the stability boundary directly:
a grid-scale mode is seeded and its growth measured, giving $-82.1$, $+55.1$ and
$-1.66$ per unit time against von Neumann predictions of $-82.3$, $+55.0$ and
$-1.66$. Panel (c) shows the price: for the same problem at $N = 256$, explicit
diffusion needs 1639 steps, Strang splitting 80, Crank--Nicolson 40.

#figure(
  image("../figures/fig01_validation.png", width: 93%),
  caption: [(a) Spatial convergence against the exact Cole--Hopf solution.
  (b) Amplification of a seeded grid-scale mode near the explicit stability
  limit: the measured rates match von Neumann theory to $0.3%$.
  (c) Number of time steps needed at $N = 256$.],
) <fig:validation>

== Breaking, and the thickness of what is left

@fig:shock (a) is the space--time diagram for $nu = 2 times 10^(-3)$. The wave
steepens, breaks at $t_s = 1$ and thereafter consists of a shock of fixed
thickness. The maximum strain (@fig:scaling, a) follows $1\/(t_s - t)$ before
breaking --- measured values $1.249$ and $1.993$ against $1.250$ and $2.000$ at
$t = 0.2$ and $0.5$ --- and saturates afterwards at $U^2 \/ 2 nu = 250$; the
measured value is $247.5$. The asymptotic profile is indistinguishable from the
exact $tanh$ solution (@fig:shock, b). Fitting the thickness over five
viscosities gives $delta prop nu^(0.985)$ with a prefactor $5$--$10%$ above
$4 nu \/ U$, consistent with the finite amplitude of the decaying sawtooth.

#figure(
  image("../figures/fig02_shock_formation.png", width: 93%),
  caption: [(a) Space--time diagram; characteristics cross at $t_s = 1$.
  (b) The asymptotic profile at $t = 2.2$ against the exact steady shock.],
) <fig:shock>

#figure(
  image("../figures/fig03_shock_scaling.png", width: 93%),
  caption: [(a) Strain rate: the inviscid divergence $1\/(t_s - t)$ and the
  viscous saturation $U^2 \/ 2 nu$. (b) Shock thickness against viscosity; the
  dashed line is $4 nu \/ U$ with no free parameter.],
) <fig:scaling>

== The spectrum: a power law with a movable end

A periodic shock train is a sawtooth, and the Fourier coefficients of a sawtooth
decay as $k^(-2)$. @fig:spectra (a) confirms this over the whole inertial range,
with a fitted slope of $-2.030 plus.minus 0.001$, and shows that the three
viscosities are indistinguishable there: _viscosity does not change the shape of
the spectrum, only where it ends_. The compensated spectrum (@fig:spectra, b)
makes the end visible as a knee, and (@fig:spectra, c) shows that the cut-off
scales as $k_"cut" prop nu^(-0.971)$, that is, as $1\/ delta$.

#figure(
  image("../figures/fig04_spectra.png", width: 93%),
  caption: [(a) Energy spectra. (b) Compensated spectrum $E(k) k^2$: flat means
  $k^(-2)$, and the knee is the dissipative cut-off. (c) Cut-off against
  viscosity.],
) <fig:spectra>

== Burgulence: many shocks, one exponent

With random band-limited initial data the flow becomes a train of shocks of
random strength. The structure functions (@fig:zeta, a) are then power laws over
the range between the shock thickness and the box, and the measured exponents are
$ zeta_p = 0.52, thin 0.92, thin 1.00, thin 0.94, thin 0.92 quad
  "for" quad p = 0.5, thin 1, thin 2, thin 3, thin 4 , $
against the bifractal prediction $zeta_p = op("min")(p, 1)$ and the Kolmogorov
prediction $zeta_p = p\/3$. The two are a factor of three apart at $p = 4$,
where the measurement has a $5%$ uncertainty, so the data discriminate between
them decisively. The $8%$ deficit at $p = 1$ is the known finite-box effect: at
separations comparable to the integral scale, $S_1$ saturates and the local slope
falls.

#figure(
  image("../figures/fig05_structure_functions.png", width: 93%),
  caption: [(a) Structure functions for 24 realisations; the shaded band is the
  fitting window. (b) Measured exponents against the two competing predictions.],
) <fig:zeta>

== How much of the viscosity is the grid?

The truncation error of first-order upwind is a diffusion term,
$ nu_"num" = U thin Delta x thin (1 - C) \/ 2, quad C = U Delta t \/ Delta x . $ <eq:nunum>
Measuring it directly, from the damping of a single Fourier mode in a pure
advection problem, reproduces @eq:nunum to $1%$ over five decades in $Delta x$,
with a fitted exponent of $Delta x^(1.003)$ (@fig:grid, a). For a physical
$nu = 2 times 10^(-3)$ this means that at $N = 256$ the scheme is $3.7$ times
more viscous than the fluid, and at $N = 4096$, where five points fit inside
$4 nu \/ U$, it is still $23%$ too viscous (@fig:grid, b). The practical
criterion is therefore not a resolution in points per wavelength but a
resolution in _points per shock_: $Delta x lt.eq 0.2 thin delta$.

#figure(
  image("../figures/fig06_grid_viscosity.png", width: 93%),
  caption: [(a) Numerical viscosity of upwind, measured from mode damping,
  against @eq:nunum. (b) The same quantity compared with a physical
  $nu = 2 times 10^(-3)$: below $N approx 1024$ the grid is the more viscous
  fluid.],
) <fig:grid>

= Discussion

Three conclusions survive the numbers. First, the convergence order of a scheme
says very little about its behaviour at a shock: all four schemes have a known
order, and yet they disagree qualitatively about what a shock looks like
(smeared by upwind, ringed by the two centred schemes). Second, integral
quantities are much more robust than local ones: the $k^(-2)$ slope and the
bifractal exponents are insensitive to viscosity over three orders of magnitude,
while the shock thickness itself is not. Third, numerical dissipation is not a
small correction that can be neglected at high resolution --- it is a _measurable_
quantity with a known formula, and the honest practice is to compute it rather
than hope it is small.

The main limitation is the range of scales. With $N = 4096$ the inertial range
spans about one and a half decades, enough to separate $k^(-2)$ from $k^(-3)$ but
not enough to resolve the sub-leading corrections to $zeta_p$; and the decaying
nature of the flow means that the statistics are taken at a single time rather
than in a statistically stationary state. Forcing the equation would fix the
second issue at the cost of a second control parameter.

= Conclusion

The one-dimensional Navier--Stokes equation is small enough to be solved five
ways and rich enough that the five answers differ, in a controlled and
quantifiable manner. By measuring the schemes against exact solutions at every
step --- convergence order, breaking time, shock thickness, spectrum and
structure functions --- the project separates what the fluid does from what the
grid does, and finds the boundary between them at roughly five grid points per
shock.

= Code and reproducibility

All the numerical work lives in `project_burgulence/`, and every number quoted
above is written to `results/` as plain text, JSON or CSV at run time, so the
report can be checked against the code without re-running anything. The solver
core is four modules: `fftlib.py` (radix-2 FFT and its inverse), `solvers.py`
(the five schemes), `exact.py` (Cole--Hopf reference, characteristics, steady
shock) and `diagnostics.py` (shock thickness, spectra, structure functions and
their fit windows). The five experiments correspond one-to-one with the sections
above, `run_all.py` runs them in order, and `tests/test_core.py` checks fifteen
invariants --- symmetry, conservation, the residual of the exact solution, the
FFT against a naive DFT, the observed order of each scheme --- in under three
seconds.

Only array arithmetic uses NumPy. The FFT, the tridiagonal solver, the time
integrators and the spectral derivative are all written for this project,
following the `mylib.py` convention of the course.

#align(center)[
  #table(
    columns: (auto, auto, auto),
    align: (left, left, left),
    stroke: none,
    inset: (x: 6pt, y: 2.2pt),
    table.hline(stroke: 0.7pt),
    [*script*], [*what it measures*], [*figure*],
    table.hline(stroke: 0.4pt),
    [`e1_validation.py`], [convergence order, von Neumann growth, cost], [1],
    [`e2_shock_formation.py`], [space--time diagram, strain, shock thickness], [2, 3],
    [`e3_spectra.py`], [energy spectra, dissipative cut-off], [4],
    [`e4_structure_functions.py`], [$S_p(r)$ and the exponents $zeta_p$], [5],
    [`e5_grid_viscosity.py`], [numerical viscosity from mode damping], [6],
    table.hline(stroke: 0.7pt),
  )
]

#v(0.4em)
#line(length: 100%, stroke: 0.4pt)
#v(-0.3em)
#text(size: 8.6pt)[
  *References.*
  #par(hanging-indent: 1.1em)[
  [1] J. M. Burgers, "A mathematical model illustrating the theory of
  turbulence", _Adv. Appl. Mech._ *1*, 171 (1948). \
  [2] J. Bec and K. Khanin, "Burgers turbulence", _Phys. Rep._ *447*, 1 (2007). \
  [3] G. B. Whitham, _Linear and Nonlinear Waves_ (Wiley, 1974), ch. 2--4. \
  [4] R. J. LeVeque, _Finite Volume Methods for Hyperbolic Problems_ (Cambridge
  University Press, 2002), ch. 4, 6. \
  [5] J. P. Boyd, _Chebyshev and Fourier Spectral Methods_, 2nd ed. (Dover,
  2001), ch. 2. \
  [6] J. D. Cole, "On a quasi-linear parabolic equation occurring in
  aerodynamics", _Quart. Appl. Math._ *9*, 225 (1951). \
  ]
]
